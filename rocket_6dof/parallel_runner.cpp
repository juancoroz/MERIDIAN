//  parallel_runner.cpp -- Subprocess-based parallel execution
//
//  Spawns N mission_worker subprocesses via fork()/execvp(), each
//  consuming a chunk of the input sample matrix.  Aggregates the
//  per-row Y values into a single ParallelResult.
//
//  Worker IPC: temp CSV files in tmp_dir.  Each worker writes its
//  output CSV to a unique path; the orchestrator concatenates after
//  all workers finish.  Worker stderr is left attached to the
//  parent's stderr so failures are visible.

#include "parallel_runner.h"

#include <algorithm>
#include <cerrno>
#include <chrono>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <fcntl.h>
#include <fstream>
#include <thread>
#include <signal.h>
#include <sstream>
#include <string>
#include <sys/types.h>
#include <sys/wait.h>
#include <unistd.h>
#include <vector>

namespace rocket6dof {

namespace {

// Build the chunk CSV file for a worker, returning its path.
// Format: row_id,path1,path2,...  (one row per sample in the chunk)
std::string write_chunk_csv(const std::string& tmp_dir,
                            int worker_idx,
                            pid_t parent_pid,
                            const std::vector<std::string>& paths,
                            const std::vector<std::vector<double>>& sample_matrix,
                            int row_lo, int row_hi)
{
    std::ostringstream oss;
    oss << tmp_dir << "/r6dof_par_" << parent_pid << "_w" << worker_idx << "_in.csv";
    std::string path = oss.str();

    std::FILE* f = std::fopen(path.c_str(), "w");
    if (!f) {
        std::fprintf(stderr, "parallel: cannot open %s for writing: %s\n",
                     path.c_str(), std::strerror(errno));
        return "";
    }
    std::fprintf(f, "row_id");
    for (const auto& p : paths) std::fprintf(f, ",%s", p.c_str());
    std::fputc('\n', f);
    for (int r = row_lo; r < row_hi; ++r) {
        std::fprintf(f, "%d", r);
        for (double v : sample_matrix[r]) std::fprintf(f, ",%.17g", v);
        std::fputc('\n', f);
    }
    std::fclose(f);
    return path;
}

std::string output_chunk_path(const std::string& tmp_dir,
                              int worker_idx, pid_t parent_pid)
{
    std::ostringstream oss;
    oss << tmp_dir << "/r6dof_par_" << parent_pid << "_w" << worker_idx << "_out.csv";
    return oss.str();
}

// Split N rows into n_workers contiguous chunks of ~equal size.
// Returns vector of (lo, hi) row index pairs; chunks empty when
// n_workers > N (extra workers get empty ranges and are skipped).
std::vector<std::pair<int,int>> split_chunks(int n_rows, int n_workers) {
    std::vector<std::pair<int,int>> out;
    int per_chunk = (n_rows + n_workers - 1) / n_workers;  // ceil
    for (int w = 0; w < n_workers; ++w) {
        int lo = w * per_chunk;
        int hi = std::min(lo + per_chunk, n_rows);
        if (lo >= n_rows) break;
        out.emplace_back(lo, hi);
    }
    return out;
}

// Read a worker's output CSV into a map row_id -> output vector.
// Expected header: row_id,<output_names>,nan_flag
bool read_worker_output(const std::string& path,
                        const std::vector<std::string>& output_names,
                        std::vector<std::vector<double>>& row_outputs,
                        std::vector<bool>& nan_flags)
{
    std::ifstream f(path);
    if (!f) {
        std::fprintf(stderr, "parallel: cannot open %s for reading\n", path.c_str());
        return false;
    }
    std::string line;
    if (!std::getline(f, line)) {
        std::fprintf(stderr, "parallel: empty output file %s\n", path.c_str());
        return false;
    }
    // Parse header
    std::vector<std::string> header;
    {
        std::string cur;
        for (char c : line) {
            if (c == ',') { header.push_back(cur); cur.clear(); }
            else if (c == '\r') {}
            else cur.push_back(c);
        }
        header.push_back(cur);
    }
    if (header.empty() || header[0] != "row_id") {
        std::fprintf(stderr, "parallel: bad header in %s (expected row_id first)\n", path.c_str());
        return false;
    }
    // Map output_name -> column index in worker CSV
    std::vector<int> out_col_idx(output_names.size(), -1);
    int nan_col = -1;
    for (size_t i = 1; i < header.size(); ++i) {
        if (header[i] == "nan_flag") { nan_col = static_cast<int>(i); continue; }
        for (size_t o = 0; o < output_names.size(); ++o) {
            if (header[i] == output_names[o]) {
                out_col_idx[o] = static_cast<int>(i);
                break;
            }
        }
    }
    for (size_t o = 0; o < output_names.size(); ++o) {
        if (out_col_idx[o] < 0) {
            std::fprintf(stderr, "parallel: output '%s' missing from worker CSV\n",
                         output_names[o].c_str());
            return false;
        }
    }

    while (std::getline(f, line)) {
        if (line.empty()) continue;
        // Parse cells
        std::vector<std::string> cells;
        std::string cur;
        for (char c : line) {
            if (c == ',') { cells.push_back(cur); cur.clear(); }
            else if (c == '\r') {}
            else cur.push_back(c);
        }
        cells.push_back(cur);
        if (cells.size() != header.size()) {
            std::fprintf(stderr, "parallel: row column-count mismatch in %s (got %zu, expected %zu)\n",
                         path.c_str(), cells.size(), header.size());
            return false;
        }
        int rid = std::atoi(cells[0].c_str());
        if (rid < 0 || rid >= static_cast<int>(row_outputs.size())) {
            std::fprintf(stderr, "parallel: row_id %d out of range\n", rid);
            return false;
        }
        for (size_t o = 0; o < output_names.size(); ++o) {
            row_outputs[rid][o] = std::strtod(cells[out_col_idx[o]].c_str(), nullptr);
        }
        if (nan_col >= 0) {
            nan_flags[rid] = (std::atoi(cells[nan_col].c_str()) != 0);
        }
    }
    return true;
}

} // anon

ParallelResult run_chunked_parallel(
    const std::string& config_path,
    const std::vector<std::string>& input_paths,
    const std::vector<std::vector<double>>& sample_matrix,
    const std::vector<std::string>& output_names,
    int n_workers,
    const std::string& worker_bin,
    const std::string& tmp_dir,
    const std::string& progress_label)
{
    const int n_rows = static_cast<int>(sample_matrix.size());
    const int n_out  = static_cast<int>(output_names.size());
    ParallelResult result;
    result.ok = false;
    result.elapsed_s = 0.0;
    result.row_outputs.assign(n_rows, std::vector<double>(n_out, std::nan("")));
    result.nan_flags.assign(n_rows, false);

    if (n_workers < 1) n_workers = 1;
    auto chunks = split_chunks(n_rows, n_workers);
    if (chunks.empty()) {
        // n_rows == 0: trivially done
        result.ok = true;
        return result;
    }

    pid_t parent_pid = ::getpid();

    // Write per-worker input chunks first; record paths for cleanup
    struct WorkerJob {
        int idx;
        int row_lo, row_hi;
        std::string in_path;
        std::string out_path;
        pid_t pid;
    };
    std::vector<WorkerJob> jobs;
    jobs.reserve(chunks.size());
    for (size_t i = 0; i < chunks.size(); ++i) {
        WorkerJob j;
        j.idx     = static_cast<int>(i);
        j.row_lo  = chunks[i].first;
        j.row_hi  = chunks[i].second;
        j.in_path = write_chunk_csv(tmp_dir, j.idx, parent_pid,
                                    input_paths, sample_matrix,
                                    j.row_lo, j.row_hi);
        j.out_path = output_chunk_path(tmp_dir, j.idx, parent_pid);
        j.pid = -1;
        if (j.in_path.empty()) {
            std::fprintf(stderr, "parallel: aborting (failed to write chunk %d)\n", j.idx);
            return result;
        }
        jobs.push_back(j);
    }

    std::printf("[%s] %d rows split across %d workers (~%d rows/worker)\n",
                progress_label.c_str(), n_rows, static_cast<int>(jobs.size()),
                (n_rows + static_cast<int>(jobs.size()) - 1) / static_cast<int>(jobs.size()));
    // Flush so any pipe consumer (notably the GUI's subprocess wrapper)
    // sees this line before the workers start emitting their PROGRESS
    // heartbeats.  Without this, block-buffered stdout would withhold
    // the line until the parent does its next big write or exits.
    std::fflush(stdout);

    auto t0 = std::chrono::steady_clock::now();

    // Fork all workers
    for (auto& j : jobs) {
        pid_t pid = ::fork();
        if (pid < 0) {
            std::fprintf(stderr, "parallel: fork() failed: %s\n", std::strerror(errno));
            // Try to clean up: kill any spawned children, unlink temp files
            for (const auto& jj : jobs) {
                if (jj.pid > 0) ::kill(jj.pid, SIGTERM);
                if (!jj.in_path.empty()) ::unlink(jj.in_path.c_str());
                ::unlink(jj.out_path.c_str());
            }
            return result;
        }
        if (pid == 0) {
            // Child: exec the worker.  argv: <bin> <config> <chunk_in> <chunk_out>
            // Suppress worker stdout so the progress lines from many workers
            // don't interleave; stderr stays attached for diagnostics.
            int devnull = ::open("/dev/null", O_WRONLY);
            if (devnull >= 0) { ::dup2(devnull, 1); ::close(devnull); }

            std::vector<char*> argv;
            argv.push_back(const_cast<char*>(worker_bin.c_str()));
            argv.push_back(const_cast<char*>(config_path.c_str()));
            argv.push_back(const_cast<char*>(j.in_path.c_str()));
            argv.push_back(const_cast<char*>(j.out_path.c_str()));
            argv.push_back(nullptr);
            ::execvp(worker_bin.c_str(), argv.data());
            // If exec fails, the child must exit non-zero
            std::fprintf(stderr, "parallel: execvp(%s) failed: %s\n",
                         worker_bin.c_str(), std::strerror(errno));
            std::_Exit(127);
        }
        j.pid = pid;
    }

    // Wait for all workers
    int n_failed = 0;
    for (auto& j : jobs) {
        int status = 0;
        pid_t got = ::waitpid(j.pid, &status, 0);
        if (got != j.pid) {
            std::fprintf(stderr, "parallel: waitpid(%d) failed: %s\n",
                         (int)j.pid, std::strerror(errno));
            n_failed++;
            continue;
        }
        if (!WIFEXITED(status) || WEXITSTATUS(status) != 0) {
            std::fprintf(stderr, "parallel: worker %d (pid %d) exited with %d\n",
                         j.idx, (int)j.pid,
                         WIFEXITED(status) ? WEXITSTATUS(status) : -1);
            n_failed++;
        }
    }

    auto t1 = std::chrono::steady_clock::now();
    result.elapsed_s = std::chrono::duration<double>(t1 - t0).count();

    if (n_failed > 0) {
        std::fprintf(stderr, "parallel: %d worker(s) failed; aborting\n", n_failed);
        for (const auto& j : jobs) {
            if (!j.in_path.empty()) ::unlink(j.in_path.c_str());
            ::unlink(j.out_path.c_str());
        }
        return result;
    }

    // Read worker outputs
    for (const auto& j : jobs) {
        if (!read_worker_output(j.out_path, output_names,
                                result.row_outputs, result.nan_flags)) {
            std::fprintf(stderr, "parallel: failed to read worker %d output\n", j.idx);
            // Best-effort cleanup before returning
            for (const auto& jj : jobs) {
                ::unlink(jj.in_path.c_str());
                ::unlink(jj.out_path.c_str());
            }
            return result;
        }
    }

    // Cleanup temp files
    for (const auto& j : jobs) {
        ::unlink(j.in_path.c_str());
        ::unlink(j.out_path.c_str());
    }

    std::printf("[%s] done in %.2fs (%.1f ms/run effective; %.1f ms/run per-worker)\n",
                progress_label.c_str(),
                result.elapsed_s,
                1000.0 * result.elapsed_s / std::max(1, n_rows),
                1000.0 * result.elapsed_s * jobs.size() / std::max(1, n_rows));

    result.ok = true;
    return result;
}

ResolvedWorkers resolve_n_workers(int requested, int auto_cap) {
    ResolvedWorkers r;
    if (requested >= 1) {
        // Explicit user choice -- respect it without applying the auto cap.
        // The user took the trouble to type a number; they know their env.
        r.n = requested;
        char buf[64];
        std::snprintf(buf, sizeof(buf), "user-specified: %d", requested);
        r.source = buf;
        return r;
    }
    if (requested < 0) {
        // Negative is nonsense; clamp to 1 and tell the user.
        r.n = 1;
        r.source = "invalid (clamped from negative)";
        return r;
    }
    // requested == 0 -> auto-detect.
    unsigned hw = std::thread::hardware_concurrency();
    if (hw == 0) {
        // Some containers and odd platforms report 0.  Fall back to 1.
        r.n = 1;
        r.source = "auto: hardware_concurrency() returned 0, using 1";
        return r;
    }
    int detected = static_cast<int>(hw);
    char buf[128];
    if (detected > auto_cap) {
        r.n = auto_cap;
        std::snprintf(buf, sizeof(buf),
                      "auto: %d logical cores detected, capped at %d "
                      "(set n_workers explicitly to override)",
                      detected, auto_cap);
    } else {
        r.n = detected;
        std::snprintf(buf, sizeof(buf), "auto: %d logical cores detected", detected);
    }
    r.source = buf;
    return r;
}

} // namespace rocket6dof
