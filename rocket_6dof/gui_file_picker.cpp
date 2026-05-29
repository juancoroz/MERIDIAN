//  gui_file_picker.cpp -- Minimal ImGui-based file open/save dialog
//
//  Pure ImGui; no system file dialog (Zenity/KDialog) required.  POSIX
//  dirent + stat for directory listing.  Looks consistent with the rest
//  of the GUI since it's drawn by the same renderer.
//
//  Usage from caller code:
//
//      static FilePicker picker;
//      if (ImGui::Button("Load")) picker.open("Load config", false);
//      std::string chosen;
//      if (picker.draw(chosen)) {
//          // chosen is the absolute path the user picked
//          load_config_file(chosen, cfg);
//      }
//
//  Behavior:
//    - The dialog starts in the current working directory the first
//      time it's opened.  Subsequent opens remember the last directory
//      the user navigated to.
//    - Files are filtered by extension (default ".json").
//    - Directories show first, sorted; files second, sorted.
//    - "Open" / "Save" button is only enabled when a file is selected
//      or (in Save mode) a non-empty filename has been typed.
//    - For Save mode, if the typed name doesn't end in ".json", the
//      caller is responsible for appending it -- this dialog returns
//      the raw path the user chose.

#include "imgui.h"
#include "gui_file_picker.h"

#include <algorithm>
#include <cstring>
#include <cstdio>
#include <dirent.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <unistd.h>
#include <climits>

namespace gui {

// Helpers

// Strip trailing slashes from a path (except for "/").  Helps when
// concatenating "dir" + "/" + "name".
static std::string strip_trailing_slash(const std::string& p) {
    if (p.size() <= 1) return p;
    std::string s = p;
    while (s.size() > 1 && s.back() == '/') s.pop_back();
    return s;
}

// Return the parent directory of a path.  Examples:
//   "/foo/bar/baz" -> "/foo/bar"
//   "/foo"         -> "/"
//   "/"            -> "/"
static std::string parent_dir(const std::string& p) {
    std::string s = strip_trailing_slash(p);
    size_t pos = s.find_last_of('/');
    if (pos == std::string::npos) return ".";
    if (pos == 0) return "/";
    return s.substr(0, pos);
}

// Convert a path to canonical absolute form, resolving "." and "..".
// Returns the input unchanged if realpath() fails.  realpath() requires
// the path to exist for the full chain, so this is only safe for
// directories we've just successfully listed.
static std::string canonicalize(const std::string& p) {
    char buf[PATH_MAX];
    if (realpath(p.c_str(), buf)) {
        return std::string(buf);
    }
    return p;
}

// Get current working directory as a std::string.
static std::string get_cwd() {
    char buf[PATH_MAX];
    if (getcwd(buf, sizeof(buf))) return std::string(buf);
    return ".";
}

// Does the path end with the given suffix?  Case-sensitive.
static bool ends_with(const std::string& s, const std::string& suffix) {
    if (s.size() < suffix.size()) return false;
    return std::equal(suffix.rbegin(), suffix.rend(), s.rbegin());
}

// Split a comma-separated list of extensions ("a,b,c") into individual
// strings.  Used by the filter to accept multiple acceptable
// extensions like ".asc,.txt".
static std::vector<std::string> split_csv(const std::string& s) {
    std::vector<std::string> out;
    if (s.empty()) return out;
    size_t i = 0;
    while (i < s.size()) {
        size_t j = s.find(',', i);
        if (j == std::string::npos) j = s.size();
        // Strip surrounding whitespace
        size_t a = i, b = j;
        while (a < b && (s[a] == ' ' || s[a] == '\t')) ++a;
        while (b > a && (s[b-1] == ' ' || s[b-1] == '\t')) --b;
        if (a < b) out.push_back(s.substr(a, b - a));
        i = j + 1;
    }
    return out;
}

// FilePicker implementation

FilePicker::FilePicker() {
    // current_dir_ stays empty until first open(); we set it then.
}

void FilePicker::open(const std::string& title, bool save_mode,
                      const std::string& default_name,
                      const std::string& filter_extension) {
    title_ = title;
    save_mode_ = save_mode;
    filter_ext_ = filter_extension;
    selected_name_.clear();
    type_buffer_[0] = '\0';
    if (!default_name.empty()) {
        std::snprintf(type_buffer_, sizeof(type_buffer_), "%s",
                      default_name.c_str());
        selected_name_ = default_name;
    }
    if (current_dir_.empty()) {
        current_dir_ = get_cwd();
    }
    refresh_entries();
    should_open_popup_ = true;
}

void FilePicker::refresh_entries() {
    entries_.clear();

    DIR* d = opendir(current_dir_.c_str());
    if (!d) return;

    struct dirent* de;
    while ((de = readdir(d)) != nullptr) {
        std::string name = de->d_name;
        // Skip "." (always; "." is current dir, no need to navigate to
        // it).  Keep ".." (so user can go up).
        if (name == ".") continue;
        // Skip hidden files except "..".  Engineers generally don't
        // store configs in dotfiles.
        if (name.size() > 1 && name[0] == '.' && name != "..") continue;

        std::string full = current_dir_ + "/" + name;
        struct stat st;
        if (stat(full.c_str(), &st) != 0) continue;

        Entry e;
        e.name = name;
        e.is_dir = S_ISDIR(st.st_mode);
        // Filter files by extension(s); always keep directories
        if (!e.is_dir) {
            if (!filter_ext_.empty()) {
                auto exts = split_csv(filter_ext_);
                bool match = false;
                for (const auto& ext : exts) {
                    if (ends_with(name, ext)) { match = true; break; }
                }
                if (!match) continue;
            }
        }
        entries_.push_back(std::move(e));
    }
    closedir(d);

    // Directories first, then files, both alphabetical
    std::sort(entries_.begin(), entries_.end(),
        [](const Entry& a, const Entry& b) {
            if (a.is_dir != b.is_dir) return a.is_dir > b.is_dir;
            return a.name < b.name;
        });
}

bool FilePicker::draw(std::string& out_path) {
    bool confirmed = false;

    if (should_open_popup_) {
        ImGui::OpenPopup(title_.c_str());
        should_open_popup_ = false;
    }

    // Center popup on the main viewport.
    ImVec2 center = ImGui::GetMainViewport()->GetCenter();
    ImGui::SetNextWindowPos(center, ImGuiCond_Appearing, ImVec2(0.5f, 0.5f));
    ImGui::SetNextWindowSize(ImVec2(620, 480), ImGuiCond_Appearing);

    if (ImGui::BeginPopupModal(title_.c_str(), nullptr,
                               ImGuiWindowFlags_NoSavedSettings)) {

        // ---- Current directory label + Up button ----
        ImGui::Text("Path: %s", current_dir_.c_str());
        ImGui::SameLine();
        if (ImGui::SmallButton("Up")) {
            current_dir_ = canonicalize(parent_dir(current_dir_));
            refresh_entries();
        }

        ImGui::Separator();

        // ---- Entry list ----
        // Reserve space for the typed-name row + buttons below.
        float footer_h = ImGui::GetFrameHeightWithSpacing() * 2 + 8.0f;
        if (ImGui::BeginChild("file_list",
                              ImVec2(0, -footer_h),
                              true)) {
            for (size_t i = 0; i < entries_.size(); ++i) {
                const Entry& e = entries_[i];
                // Prefix directories with "/" so they're visually
                // distinct from files.  ImGui doesn't ship with
                // folder icons by default; this is the lightest visual
                // cue that costs no extra assets.
                std::string label = e.is_dir ? ("[ ] " + e.name + "/")
                                             : ("    " + e.name);

                bool is_selected = (!e.is_dir &&
                                    selected_name_ == e.name);
                if (ImGui::Selectable(label.c_str(), is_selected,
                                      ImGuiSelectableFlags_AllowDoubleClick)) {
                    if (e.is_dir) {
                        // Single click navigates into directories.
                        // A "click to select dir, double-click to enter"
                        // distinction adds complexity for no benefit
                        // in this UX.
                        std::string new_dir;
                        if (e.name == "..") {
                            new_dir = parent_dir(current_dir_);
                        } else {
                            new_dir = current_dir_ + "/" + e.name;
                        }
                        current_dir_ = canonicalize(new_dir);
                        refresh_entries();
                        // Clear selection on directory change
                        selected_name_.clear();
                        type_buffer_[0] = '\0';
                    } else {
                        // File selection: store the name, populate the
                        // text field so the user can edit before save.
                        selected_name_ = e.name;
                        std::snprintf(type_buffer_, sizeof(type_buffer_),
                                      "%s", e.name.c_str());
                        // Double-click on a file confirms the dialog
                        // (Open mode only -- in Save mode it would be
                        // surprising to overwrite without an explicit
                        // confirm).
                        if (!save_mode_ && ImGui::IsMouseDoubleClicked(0)) {
                            confirmed = true;
                        }
                    }
                }
            }
        }
        ImGui::EndChild();

        // ---- Typed name field + buttons ----
        // In Save mode, the field is freely editable.  In Open mode,
        // it's a display of the current selection (still editable in
        // case the user wants to type an exact filename).
        ImGui::SetNextItemWidth(-1.0f);
        if (ImGui::InputText("##typed_name", type_buffer_,
                             sizeof(type_buffer_))) {
            selected_name_ = type_buffer_;
        }

        // OK button label depends on mode.
        const char* ok_label = save_mode_ ? "Save" : "Open";
        bool can_confirm = (selected_name_.length() > 0);

        if (!can_confirm) {
            ImGui::BeginDisabled();
        }
        if (ImGui::Button(ok_label, ImVec2(120, 0))) {
            confirmed = true;
        }
        if (!can_confirm) {
            ImGui::EndDisabled();
        }
        ImGui::SameLine();
        if (ImGui::Button("Cancel", ImVec2(120, 0))) {
            ImGui::CloseCurrentPopup();
        }

        if (confirmed) {
            out_path = current_dir_ + "/" + selected_name_;
            ImGui::CloseCurrentPopup();
        }

        ImGui::EndPopup();
    }

    return confirmed;
}

}  // namespace gui
