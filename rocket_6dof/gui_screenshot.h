//  gui_screenshot.h -- save the GL framebuffer (or a region of it) to PNG
//
//  Used by the Plots tab to implement "Export plots as PNG".  Reads
//  the current framebuffer with glReadPixels and writes via
//  stb_image_write.
//
//  Coordinate convention: (x, y) is the BOTTOM-LEFT of the region in
//  OpenGL pixel space (origin at lower-left).  ImPlot::GetPlotPos()
//  returns top-left in ImGui space (origin at upper-left), so callers
//  need to convert.
//
//  Usage from Plots tab:
//      // Get plot rect from ImPlot::GetPlotPos / GetPlotSize
//      // OR just save the whole window: x=0, y=0, w=fb_w, h=fb_h
//      gui::save_framebuffer_png("altitude.png", x, y, w, h);
//
//  Save calls must happen AFTER ImGui::Render() but BEFORE
//  glfwSwapBuffers(); the back buffer is the one read at that point.
#ifndef ROCKET6DOF_GUI_SCREENSHOT_H
#define ROCKET6DOF_GUI_SCREENSHOT_H

#include <string>

namespace gui {

// Save a rectangular region of the current framebuffer to a PNG file.
// (x, y) is the bottom-left of the region in framebuffer pixels.
// Returns true on success.
bool save_framebuffer_png(const std::string& path,
                          int x, int y, int width, int height);

}  // namespace gui

#endif
