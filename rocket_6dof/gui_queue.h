//  gui_queue.h -- Queue tab: run multiple missions sequentially
//
//  State:
//    - vector of QueueItem (path + status)
//    - worker thread that runs them in sequence
//
//  UI:
//    - List of queued configs with per-item status
//    - "Add config..." button (file picker)
//    - "Remove" button per item
//    - "Run queue" / "Cancel queue" buttons
//
//  The worker thread reuses run_mission_ext (the same in-process entry
//  point the single-shot Run button uses).  Each queue entry runs
//  serially; mission_progress::current_t/total_t/is_running publish
//  for the existing header progress bar, so the user gets per-mission
//  progress without extra UI.
#ifndef ROCKET6DOF_GUI_QUEUE_H
#define ROCKET6DOF_GUI_QUEUE_H

#include "gui_file_picker.h"

#include <atomic>
#include <memory>
#include <string>
#include <thread>
#include <vector>

namespace gui {

struct QueueItem {
    std::string path;
    enum Status { Queued, Running, Done, Failed, Cancelled };
    Status      status      = Queued;
    double      elapsed_sec = 0.0;
    int         rc          = 0;
    // Final altitude parsed from the summary, if available.  Set on
    // success only.  -1 means "not parsed".
    double      final_alt_km = -1.0;
};

struct QueueTabState {
    std::vector<QueueItem> items;

    // File picker for the "Add config..." button.
    FilePicker             picker;

    // Worker thread + control flags.  The thread iterates over items,
    // updating their status as it goes.  cancel_requested is set by
    // the Cancel button; the worker checks it between items (cannot
    // cancel mid-mission because run_mission_ext is synchronous).
    std::unique_ptr<std::thread> worker;
    std::atomic<bool>            worker_done       {false};
    std::atomic<bool>            cancel_requested  {false};

    // Display text describing the queue's last terminal state, e.g.
    // "Queue done: 3 ok, 1 failed".
    std::string  summary_text;
};

// Render the Queue tab.  Returns true if the user clicked "Run queue"
// this frame; caller is responsible for spawning the worker thread
// (the threading machinery lives in mission_gui.cpp's main loop).
bool render_queue_tab(QueueTabState& q);

}  // namespace gui

#endif
