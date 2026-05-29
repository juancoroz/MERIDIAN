//  mission_progress.h -- shared progress publication
//
//  Used by mission_runner.cpp's Logger to publish current mission time
//  and a completion flag, and by mission_gui.cpp's worker-thread
//  machinery to read those values for the progress bar.
//
//  This is a single-writer, single-reader pattern: one running mission
//  writes, one GUI thread reads.  std::atomic<double> with default
//  memory ordering is sufficient.
//
//  Default state (not running): current_t = 0, total_t = 0,
//  is_running = false.  The caller should set total_t to the expected
//  t_end BEFORE starting the mission, so the GUI can compute a fraction.
//
//  Why a global rather than a callback parameter?  The Logger lives
//  inside mission_runner.cpp's anonymous namespace and is created
//  fresh inside run_mission().  Threading a progress callback through
//  would require changing run_mission()'s signature and the Logger's
//  constructor, both of which would ripple through CLI usage too.
//  A global atomic is the smallest possible coupling: the CLI ignores
//  it, the GUI reads it, the simulation always writes it.
#ifndef ROCKET6DOF_MISSION_PROGRESS_H
#define ROCKET6DOF_MISSION_PROGRESS_H

#include <atomic>

namespace rocket6dof {
namespace progress {

// Current mission time, updated each time the Logger's rpt() runs.
// Reset to 0 at the start of every run_mission() call.
extern std::atomic<double> current_t;

// Expected total mission duration (t_end from the config).  Set by
// run_mission() before the integration loop starts.  GUI uses
// current_t / total_t to draw a progress bar.
extern std::atomic<double> total_t;

// True while a mission is in flight.  Set true at the start of
// run_mission(), cleared at the end.  The GUI thread polls this to
// know when to stop showing the progress bar and start showing
// post-run plots.
extern std::atomic<bool> is_running;

// Set by run_mission() at completion.  0 = success, nonzero = failure.
extern std::atomic<int> last_rc;

}  // namespace progress
}  // namespace rocket6dof

#endif
