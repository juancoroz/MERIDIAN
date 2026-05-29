//  gui_file_picker.h -- header for gui_file_picker.cpp
#ifndef ROCKET6DOF_GUI_FILE_PICKER_H
#define ROCKET6DOF_GUI_FILE_PICKER_H

#include <string>
#include <vector>

namespace gui {

class FilePicker {
public:
    FilePicker();

    // Open the dialog.  Call this from a button handler.  The popup
    // will be drawn on subsequent draw() calls until the user confirms
    // or cancels.
    //
    // Parameters:
    //   title              -- popup window title (used as ImGui ID; must
    //                         be stable across frames for the same logical
    //                         dialog)
    //   save_mode          -- false: pick existing file; true: pick a
    //                         destination (typed names allowed)
    //   default_name       -- pre-populate the typed-name field; pass ""
    //                         to leave it blank
    //   filter_extension   -- only show files matching any of these
    //                         extensions.  Accepts a single extension
    //                         (".json") or a comma-separated list
    //                         (".asc,.txt").  Pass "" to show all files.
    void open(const std::string& title,
              bool save_mode,
              const std::string& default_name = "",
              const std::string& filter_extension = ".json");

    // Draw the popup if it's open.  Returns true on the frame the user
    // confirms; out_path is then set to the full chosen path.  Returns
    // false otherwise.
    //
    // Call this every frame from the main loop (cheap when the popup
    // is closed).
    bool draw(std::string& out_path);

private:
    struct Entry {
        std::string name;
        bool        is_dir;
    };

    void refresh_entries();

    std::string         title_;
    bool                save_mode_         = false;
    std::string         current_dir_;
    std::string         filter_ext_        = ".json";
    std::vector<Entry>  entries_;
    std::string         selected_name_;
    char                type_buffer_[512]  = {0};

    // OpenPopup must be called between BeginFrame and the popup's own
    // BeginPopupModal.  We set this from open() and consume it on the
    // next draw() call.
    bool                should_open_popup_ = false;
};

}  // namespace gui

#endif
