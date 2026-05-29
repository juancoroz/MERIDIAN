//  gui_subprocess.cpp -- Subprocess implementation
//
//  POSIX-only.  Uses fork()/execvp()/waitpid() and a pipe for stdout
//  capture.  The worker thread reads the pipe line by line, appends
//  each line to the log under a mutex, and on EOF reaps the child.
//
//  Why not popen()?  popen() doesn't give us the child's pid (so we
//  can't signal it for cancellation), and it merges stderr only via
//  shell redirection.  Going one level lower with fork/exec is the
//  same amount of code and gives us the pid for free.
#include "gui_subprocess.h"

#include <cerrno>
#include <chrono>
#include <csignal>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <fcntl.h>
#include <sys/types.h>
#include <sys/wait.h>
#include <unistd.h>

namespace gui {

namespace {
// Cap the log at this many lines so a verbose subprocess can't
// exhaust memory.  We keep the most recent lines.
constexpr size_t LOG_CAP = 10000;
}  // anon

Subprocess::Subprocess() = default;

Subprocess::~Subprocess() {
    wait();
}

void Subprocess::wait() {
    if (worker_ && worker_->joinable()) {
        worker_->join();
    }
    worker_.reset();
}

std::string Subprocess::last_line() const {
    std::lock_guard<std::mutex> lock(log_mu_);
    return last_line_;
}

std::vector<std::string> Subprocess::log_lines() const {
    std::lock_guard<std::mutex> lock(log_mu_);
    return log_;  // copy
}

void Subprocess::signal(int sig) {
    int pid = child_pid_.load();
    if (pid > 0) {
        ::kill(pid, sig);
    }
}

bool Subprocess::start(const std::vector<std::string>& argv) {
    if (running_.load()) return false;
    if (argv.empty()) return false;

    // Reset state.  Caller may be re-running on a previously-used
    // Subprocess instance.
    if (worker_ && worker_->joinable()) worker_->join();
    worker_.reset();
    {
        std::lock_guard<std::mutex> lock(log_mu_);
        log_.clear();
        last_line_.clear();
    }
    exit_code_.store(0);
    elapsed_seconds_.store(0.0);
    done_.store(false);

    // Make a copy of argv to own inside the worker (C-string pointers
    // need a backing store that outlives execvp).
    auto argv_copy = std::make_shared<std::vector<std::string>>(argv);

    running_.store(true);

    worker_.reset(new std::thread([this, argv_copy]() {
        const auto t0 = std::chrono::steady_clock::now();

        // Build C-string array for execvp.
        std::vector<char*> cargs;
        cargs.reserve(argv_copy->size() + 1);
        for (auto& a : *argv_copy) cargs.push_back(const_cast<char*>(a.c_str()));
        cargs.push_back(nullptr);

        // Pipe for child's stdout+stderr (combined into one stream).
        int pipefd[2];
        if (pipe(pipefd) != 0) {
            std::lock_guard<std::mutex> lock(log_mu_);
            last_line_ = std::string("pipe() failed: ") + std::strerror(errno);
            log_.push_back(last_line_);
            exit_code_.store(-1);
            running_.store(false);
            done_.store(true);
            return;
        }

        pid_t pid = ::fork();
        if (pid < 0) {
            ::close(pipefd[0]);
            ::close(pipefd[1]);
            std::lock_guard<std::mutex> lock(log_mu_);
            last_line_ = std::string("fork() failed: ") + std::strerror(errno);
            log_.push_back(last_line_);
            exit_code_.store(-1);
            running_.store(false);
            done_.store(true);
            return;
        }

        if (pid == 0) {
            // Child: redirect stdout and stderr to the pipe write end,
            // then exec.  No std::cerr here -- we're in a forked
            // process and need to be careful about libc state.
            ::close(pipefd[0]);
            ::dup2(pipefd[1], 1);  // stdout
            ::dup2(pipefd[1], 2);  // stderr
            ::close(pipefd[1]);
            ::execvp(cargs[0], cargs.data());
            // exec failed -- write something useful and exit.
            ::fprintf(stderr, "execvp(%s) failed: %s\n",
                      cargs[0], std::strerror(errno));
            ::_exit(127);
        }

        // Parent: close write end, save pid, read from read end.
        ::close(pipefd[1]);
        child_pid_.store(pid);

        // Read line-by-line.  Use a small FILE* wrapper for getline().
        FILE* f = ::fdopen(pipefd[0], "r");
        if (!f) {
            ::close(pipefd[0]);
            std::lock_guard<std::mutex> lock(log_mu_);
            last_line_ = "fdopen() failed";
            log_.push_back(last_line_);
            exit_code_.store(-1);
            running_.store(false);
            done_.store(true);
            return;
        }

        // Drive a small timer alongside the read so the GUI sees an
        // updating elapsed_seconds value.  We poll with a short
        // timeout to keep the elapsed counter ticking even when the
        // child is quiet.
        int fd = ::fileno(f);
        // Make the fd non-blocking so reads don't stall the elapsed
        // counter for many seconds.
        int flags = ::fcntl(fd, F_GETFL, 0);
        if (flags != -1) ::fcntl(fd, F_SETFL, flags | O_NONBLOCK);

        std::string partial;
        char buf[4096];
        while (true) {
            // Update elapsed time
            const auto now = std::chrono::steady_clock::now();
            const double secs =
                std::chrono::duration<double>(now - t0).count();
            elapsed_seconds_.store(secs);

            ssize_t n = ::read(fd, buf, sizeof(buf));
            if (n > 0) {
                partial.append(buf, n);
                // Split on '\n' into log lines.
                size_t pos = 0;
                while (true) {
                    size_t nl = partial.find('\n', pos);
                    if (nl == std::string::npos) break;
                    std::string line = partial.substr(pos, nl - pos);
                    pos = nl + 1;
                    // Strip trailing \r
                    if (!line.empty() && line.back() == '\r') line.pop_back();
                    std::lock_guard<std::mutex> lock(log_mu_);
                    if (log_.size() >= LOG_CAP) {
                        // Drop oldest 10% to amortize the shift cost
                        log_.erase(log_.begin(),
                                   log_.begin() + (LOG_CAP / 10));
                    }
                    log_.push_back(line);
                    if (!line.empty()) last_line_ = line;
                }
                partial.erase(0, pos);
            } else if (n == 0) {
                // EOF -- child closed its end of the pipe.
                break;
            } else {
                // n < 0 -- EAGAIN means no data ready yet; sleep a bit.
                if (errno == EAGAIN || errno == EWOULDBLOCK) {
                    std::this_thread::sleep_for(
                        std::chrono::milliseconds(50));
                    continue;
                }
                if (errno == EINTR) continue;
                break;
            }
        }
        // Flush any trailing partial line that lacked a newline.
        if (!partial.empty()) {
            std::lock_guard<std::mutex> lock(log_mu_);
            log_.push_back(partial);
            last_line_ = partial;
        }
        ::fclose(f);

        // Reap the child.
        int status = 0;
        while (::waitpid(pid, &status, 0) == -1 && errno == EINTR) {
            // retry
        }
        int rc;
        if (WIFEXITED(status)) {
            rc = WEXITSTATUS(status);
        } else if (WIFSIGNALED(status)) {
            // Killed by signal -- report a conventional code.
            rc = 128 + WTERMSIG(status);
        } else {
            rc = -1;
        }
        exit_code_.store(rc);

        // Final elapsed update
        const auto t1 = std::chrono::steady_clock::now();
        elapsed_seconds_.store(
            std::chrono::duration<double>(t1 - t0).count());

        child_pid_.store(-1);
        running_.store(false);
        done_.store(true);
    }));

    return true;
}

}  // namespace gui
