//  gui_plot_panel.h -- types shared between gui_config.cpp's plot
//                      panel renderer and mission_gui.cpp's main loop
//
//  PlotPanelActions is the bundle of "what did the user click this
//  frame" booleans returned by render_plot_panel(). Keeping it in a
//  shared header avoids two separate struct definitions in two TUs
//  (which would be an ODR violation).
#ifndef ROCKET6DOF_GUI_PLOT_PANEL_H
#define ROCKET6DOF_GUI_PLOT_PANEL_H

namespace gui {

struct PlotPanelActions {
    bool export_clicked         = false;  // user clicked "Export as PNG"
    bool load_compare_clicked   = false;  // user clicked "Load comparison..."
    bool clear_compare_clicked  = false;  // user clicked "Clear comparison"
};

}  // namespace gui

#endif
