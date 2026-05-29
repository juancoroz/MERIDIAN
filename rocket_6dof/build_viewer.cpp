//  build_viewer.cpp -- Embed mission_log.csv into mission_viewer.html
//
//  Reads `mission_viewer_template.html`, finds the marker
//      <!-- CSV_DATA_INJECTION_POINT -->
//  and replaces it with a <script> block defining window.CSV_DATA.
//  The CSV is wrapped in a JS template literal (backticks); backslash,
//  backtick, and dollar are escaped.
//
//  Usage:
//    build_viewer
//
//  Inputs (defaults, current dir):
//    mission_viewer_template.html
//    mission_log.csv
//  Output:
//    mission_viewer.html

#include <cstdio>
#include <fstream>
#include <sstream>
#include <string>
#include <sys/stat.h>

namespace {

bool file_exists(const std::string& path) {
    struct stat st;
    return ::stat(path.c_str(), &st) == 0;
}

bool read_file(const std::string& path, std::string& out) {
    std::ifstream f(path, std::ios::binary);
    if (!f) return false;
    std::ostringstream ss;
    ss << f.rdbuf();
    out = ss.str();
    return true;
}

std::string escape_for_js(const std::string& s) {
    std::string out;
    out.reserve(s.size() + s.size() / 16);
    for (char c : s) {
        switch (c) {
            case '\\': out += "\\\\"; break;
            case '`':  out += "\\`";  break;
            case '$':  out += "\\$";  break;
            default:   out += c;      break;
        }
    }
    return out;
}

}  // anon

int main() {
    const std::string template_path = "mission_viewer_template.html";
    const std::string csv_path      = "mission_log.csv";
    const std::string output_path   = "mission_viewer.html";

    if (!file_exists(csv_path)) {
        std::fprintf(stderr, "error: %s not found; run ./mission first\n", csv_path.c_str());
        return 1;
    }
    if (!file_exists(template_path)) {
        std::fprintf(stderr, "error: %s not found\n", template_path.c_str());
        return 1;
    }

    std::string templ, csv_raw;
    if (!read_file(template_path, templ))  { std::fprintf(stderr, "error reading %s\n", template_path.c_str()); return 1; }
    if (!read_file(csv_path,      csv_raw)) { std::fprintf(stderr, "error reading %s\n", csv_path.c_str());     return 1; }

    std::string csv_js = escape_for_js(csv_raw);
    std::ostringstream injection;
    injection << "<script>\nwindow.CSV_DATA = `" << csv_js << "`;\n</script>\n";

    const std::string marker = "<!-- CSV_DATA_INJECTION_POINT -->";
    size_t pos = templ.find(marker);
    if (pos == std::string::npos) {
        std::fprintf(stderr, "error: marker '%s' not found in template\n", marker.c_str());
        return 1;
    }
    templ.replace(pos, marker.size(), injection.str());

    std::ofstream fout(output_path, std::ios::binary);
    if (!fout) {
        std::fprintf(stderr, "error: cannot open %s for writing\n", output_path.c_str());
        return 1;
    }
    fout.write(templ.data(), static_cast<std::streamsize>(templ.size()));
    std::printf("Wrote %s (%zu bytes, %zu bytes CSV embedded)\n",
                output_path.c_str(), templ.size(), csv_js.size());
    return 0;
}
