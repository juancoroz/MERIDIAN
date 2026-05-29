//  gui_screenshot.cpp -- glReadPixels -> PNG via stb_image_write
//
//  glReadPixels returns rows bottom-to-top (origin at lower-left).
//  stb_image_write expects rows top-to-bottom (image-file convention).
//  We flip the row order after reading.
#include "gui_screenshot.h"

#define STB_IMAGE_WRITE_IMPLEMENTATION
#include "stb_image_write.h"

#include <GL/gl.h>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <vector>

namespace gui {

bool save_framebuffer_png(const std::string& path,
                          int x, int y, int width, int height) {
    if (width <= 0 || height <= 0) return false;

    // Read the framebuffer.  GL_RGBA -> 4 bytes per pixel; PNG handles
    // transparency just fine even if we don't really need it.
    const int channels = 4;
    std::vector<unsigned char> pixels(
        static_cast<size_t>(width) * static_cast<size_t>(height) * channels);

    // glReadPixels reads from the currently-bound read framebuffer.
    // After ImGui_ImplOpenGL3_RenderDrawData but before SwapBuffers,
    // the back buffer holds the just-rendered frame.
    glReadPixels(x, y, width, height, GL_RGBA, GL_UNSIGNED_BYTE, pixels.data());

    // Flip rows: glReadPixels gives bottom-to-top; PNG wants top-to-bottom.
    // Do it in-place by swapping row pairs (saves a second buffer).
    const size_t row_bytes = static_cast<size_t>(width) * channels;
    std::vector<unsigned char> tmp_row(row_bytes);
    for (int row = 0; row < height / 2; ++row) {
        unsigned char* top = pixels.data() + row * row_bytes;
        unsigned char* bot = pixels.data() + (height - 1 - row) * row_bytes;
        std::memcpy(tmp_row.data(), top, row_bytes);
        std::memcpy(top, bot, row_bytes);
        std::memcpy(bot, tmp_row.data(), row_bytes);
    }

    // stb_image_write_png returns nonzero on success.
    int rc = stbi_write_png(path.c_str(), width, height, channels,
                            pixels.data(), static_cast<int>(row_bytes));
    if (rc == 0) {
        std::fprintf(stderr, "save_framebuffer_png: stbi_write_png failed for %s\n",
                     path.c_str());
        return false;
    }
    return true;
}

}  // namespace gui
