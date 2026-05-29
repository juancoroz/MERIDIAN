//  gui_queue.cpp -- Queue tab implementation
//
//  The actual thread spawning happens in mission_gui.cpp because the
//  worker needs access to the AppState's config (for setting up the
//  temp file path, etc.) and the existing mission progress atomics.
//  This file just renders the UI and reports user actions.
#include "gui_queue.h"
#include "imgui.h"

#include <cstdio>

namespace gui {

namespace {

const char* status_text(QueueItem::Status s) {
    switch (s) {
        case QueueItem::Queued:    return "queued";
        case QueueItem::Running:   return "running";
        case QueueItem::Done:      return "done";
        case QueueItem::Failed:    return "failed";
        case QueueItem::Cancelled: return "cancelled";
    }
    return "?";
}

ImU32 status_color(QueueItem::Status s) {
    switch (s) {
        case QueueItem::Queued:    return IM_COL32(150, 150, 150, 255);
        case QueueItem::Running:   return IM_COL32(255, 220,  80, 255);
        case QueueItem::Done:      return IM_COL32( 80, 200, 120, 255);
        case QueueItem::Failed:    return IM_COL32(220,  80,  80, 255);
        case QueueItem::Cancelled: return IM_COL32(150, 100, 100, 255);
    }
    return IM_COL32(200, 200, 200, 255);
}

}  // anon

bool render_queue_tab(QueueTabState& q) {
    bool run_clicked = false;

    ImGui::TextDisabled(
        "Queue multiple mission configs to run sequentially.  Each "
        "config runs through the same in-process pipeline as the "
        "Basic-tab Run button; progress for the active mission appears "
        "in the header progress bar above.");
    ImGui::Spacing();
    ImGui::Separator();
    ImGui::Spacing();

    // ---- Buttons row ----
    bool worker_active = q.worker && !q.worker_done.load();

    // Add config picker
    if (ImGui::Button("Add config...", ImVec2(140, 30))) {
        q.picker.open("Add to queue", false, "", ".json");
    }
    std::string picked;
    if (q.picker.draw(picked)) {
        QueueItem item;
        item.path = picked;
        q.items.push_back(item);
    }

    ImGui::SameLine();

    // Run / Cancel
    if (worker_active) {
        ImGui::PushStyleColor(ImGuiCol_Button,        IM_COL32(80, 80, 80, 255));
        ImGui::PushStyleColor(ImGuiCol_ButtonHovered, IM_COL32(80, 80, 80, 255));
        ImGui::PushStyleColor(ImGuiCol_ButtonActive,  IM_COL32(80, 80, 80, 255));
        ImGui::BeginDisabled();
        ImGui::Button("Running...", ImVec2(160, 30));
        ImGui::EndDisabled();
        ImGui::PopStyleColor(3);
        ImGui::SameLine();
        if (ImGui::Button("Cancel queue", ImVec2(140, 30))) {
            q.cancel_requested.store(true);
        }
    } else {
        // Only enable Run if there's at least one queued item
        bool any_queued = false;
        for (const auto& it : q.items) {
            if (it.status == QueueItem::Queued) { any_queued = true; break; }
        }
        if (!any_queued) ImGui::BeginDisabled();
        ImGui::PushStyleColor(ImGuiCol_Button,        IM_COL32(40, 110, 60, 255));
        ImGui::PushStyleColor(ImGuiCol_ButtonHovered, IM_COL32(60, 150, 80, 255));
        ImGui::PushStyleColor(ImGuiCol_ButtonActive,  IM_COL32(30,  90, 50, 255));
        if (ImGui::Button("Run queue", ImVec2(160, 30))) {
            run_clicked = true;
        }
        ImGui::PopStyleColor(3);
        if (!any_queued) ImGui::EndDisabled();

        ImGui::SameLine();
        if (ImGui::Button("Clear all", ImVec2(140, 30))) {
            q.items.clear();
            q.summary_text.clear();
        }
        ImGui::SameLine();
        if (ImGui::Button("Reset statuses", ImVec2(140, 30))) {
            // Set every item back to Queued so Run will pick them up again.
            for (auto& it : q.items) {
                it.status      = QueueItem::Queued;
                it.elapsed_sec = 0.0;
                it.rc          = 0;
                it.final_alt_km = -1.0;
            }
            q.summary_text.clear();
            q.cancel_requested.store(false);
        }
    }

    if (!q.summary_text.empty()) {
        ImGui::Spacing();
        ImGui::TextDisabled("%s", q.summary_text.c_str());
    }

    ImGui::Spacing();
    ImGui::Separator();
    ImGui::Spacing();

    // ---- Items list ----
    if (q.items.empty()) {
        ImGui::TextDisabled("Queue is empty.  Click \"Add config...\" to "
                            "add a mission config.");
        return run_clicked;
    }

    // Layout: status indicator | index | path | elapsed | result | Remove btn
    if (ImGui::BeginTable("queue_items", 6,
                          ImGuiTableFlags_RowBg |
                          ImGuiTableFlags_BordersInnerH)) {
        ImGui::TableSetupColumn("##idx",     ImGuiTableColumnFlags_WidthFixed, 30.0f);
        ImGui::TableSetupColumn("Status",    ImGuiTableColumnFlags_WidthFixed, 90.0f);
        ImGui::TableSetupColumn("Config",    ImGuiTableColumnFlags_WidthStretch);
        ImGui::TableSetupColumn("Elapsed",   ImGuiTableColumnFlags_WidthFixed, 100.0f);
        ImGui::TableSetupColumn("Final alt", ImGuiTableColumnFlags_WidthFixed, 120.0f);
        ImGui::TableSetupColumn("##remove",  ImGuiTableColumnFlags_WidthFixed, 100.0f);

        // Use an index loop because we may want to remove an item
        // mid-iteration.
        int remove_idx = -1;
        for (size_t i = 0; i < q.items.size(); ++i) {
            const auto& item = q.items[i];

            ImGui::TableNextRow();
            ImGui::TableNextColumn();
            ImGui::Text("%zu", i + 1);

            ImGui::TableNextColumn();
            ImGui::PushStyleColor(ImGuiCol_Text, status_color(item.status));
            ImGui::TextUnformatted(status_text(item.status));
            ImGui::PopStyleColor();

            ImGui::TableNextColumn();
            ImGui::TextUnformatted(item.path.c_str());

            ImGui::TableNextColumn();
            if (item.status == QueueItem::Done ||
                item.status == QueueItem::Failed ||
                item.status == QueueItem::Cancelled) {
                ImGui::Text("%.2f s", item.elapsed_sec);
            } else {
                ImGui::TextDisabled("--");
            }

            ImGui::TableNextColumn();
            if (item.status == QueueItem::Done && item.final_alt_km >= 0.0) {
                ImGui::Text("%.2f km", item.final_alt_km);
            } else if (item.status == QueueItem::Failed) {
                ImGui::Text("rc=%d", item.rc);
            } else {
                ImGui::TextDisabled("--");
            }

            ImGui::TableNextColumn();
            // Only show Remove for items not currently running
            if (item.status != QueueItem::Running) {
                ImGui::PushID(static_cast<int>(i));
                if (ImGui::SmallButton("Remove")) {
                    remove_idx = static_cast<int>(i);
                }
                ImGui::PopID();
            }
        }

        ImGui::EndTable();

        if (remove_idx >= 0 && remove_idx < static_cast<int>(q.items.size())) {
            q.items.erase(q.items.begin() + remove_idx);
        }
    }

    return run_clicked;
}

}  // namespace gui
