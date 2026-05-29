//  gui_subprocess.h -- launch external commands, capture stdout
//
//  Used by the Sweeps tab to launch mission_mc / mission_sobol /
//  build_mc_viewer / xdg-open without blocking the GUI thread.
//
//  Each Subprocess instance runs one command on a worker thread.
//  Stdout (merged with stderr) is captured line-by-line into a shared
//  ring buffer the GUI can poll for the "last status line" display.
//  Process exit code and elapsed wall time are also published.
//
//  Lifecycle:
//      Subprocess sp;
//      sp.start({"./mission_mc", "monte_carlo.json"});
//      // ... each frame ...
//      if (sp.is_running()) {
//          show_spinner();
//          show_text(sp.last_line());
//      } else if (sp.done() && sp.exit_code() == 0) {
//          show_result();
//      }
//
//  After a process exits, the Subprocess is finished but not reset;
//  start() can be called again to launch a fresh subprocess on the
//  same instance.
#ifndef ROCKET6DOF_GUI_SUBPROCESS_H
#define ROCKET6DOF_GUI_SUBPROCESS_H

#include <atomic>
#include <chrono>
#include <memory>
#include <mutex>
#include <string>
#include <thread>
#include <vector>

namespace gui {

class Subprocess {
public:
    Subprocess();
    ~Subprocess();

    // Non-copyable, non-movable (owns a thread + atomics).
    Subprocess(const Subprocess&)            = delete;
    Subprocess& operator=(const Subprocess&) = delete;

    // Launch the command.  argv[0] is the executable; argv[1..] are
    // arguments.  Returns false if a process is already running on
    // this instance.  The command runs on a worker thread; the GUI
    // thread should not block here.
    bool start(const std::vector<std::string>& argv);

    // True while the subprocess is still executing.  Cleared as soon
    // as the worker thread observes process termination.
    bool is_running() const { return running_.load(); }

    // True if a previous start() has been called and the worker thread
    // has finished (whether the process succeeded or failed).  Note:
    // is_running() and done() can both be false simultaneously, e.g.
    // before start() is ever called.
    bool done() const { return done_.load(); }

    // Process exit code.  0 on success, nonzero on failure.  -1 if
    // the process couldn't be launched.  Only valid after done() is
    // true.
    int  exit_code() const { return exit_code_.load(); }

    // Wall clock elapsed seconds.  Updated continuously while running.
    double elapsed_seconds() const { return elapsed_seconds_.load(); }

    // The most recent line written to stdout/stderr (combined).
    // Truncated to ~256 chars in practice.  Returns empty string if
    // no output yet.
    std::string last_line() const;

    // Full captured output, up to a cap (~10000 lines).  Used by the
    // GUI to show the run's log scrollback in a child window.
    std::vector<std::string> log_lines() const;

    // If the subprocess is running, signal it (SIGTERM by default;
    // pass SIGKILL=9 for forceful termination).  Caller can poll
    // is_running() to confirm termination.
    void signal(int sig = 15);

    // Wait for the worker thread to finish (joins it if joinable).
    // Used by the destructor and during clean shutdown.  Safe to call
    // multiple times.
    void wait();

private:
    std::atomic<bool>     running_         {false};
    std::atomic<bool>     done_            {false};
    std::atomic<int>      exit_code_       {0};
    std::atomic<int>      child_pid_       {-1};
    std::atomic<double>   elapsed_seconds_ {0.0};

    mutable std::mutex                       log_mu_;
    std::vector<std::string>                 log_;       // accumulated, capped
    std::string                              last_line_; // most recent line

    std::unique_ptr<std::thread> worker_;
};

}  // namespace gui

#endif
