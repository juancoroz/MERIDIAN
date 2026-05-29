//  build_mc_viewer.cpp -- Embed MC/Sobol CSVs into mc_viewer.html
//
//  Reads `mc_viewer_template.html`, finds the marker
//      <!-- CSV_DATA_INJECTION_POINT -->
//  and replaces it with a <script> block that defines
//      window.MC_OUTPUTS_CSV
//      window.MC_INPUTS_CSV
//      window.MC_SOBOL_CSV   (only if sobol_indices.csv exists)
//  Strings are wrapped in JS template literals (backticks), so backslash,
//  backtick, and dollar must be escaped.
//
//  Usage:
//    build_mc_viewer                                              # defaults
//    build_mc_viewer <outputs.csv> <inputs.csv>
//    build_mc_viewer <outputs.csv> <inputs.csv> <output.html>
//    build_mc_viewer <outputs.csv> <inputs.csv> <output.html> <sobol.csv>
//
//  Defaults match the previous Python script:
//    outputs csv : monte_carlo_outputs.csv
//    inputs  csv : monte_carlo_inputs.csv
//    output html : mc_viewer.html
//    sobol   csv : sobol_indices.csv   (embedded if present; otherwise
//                                       the viewer hides the Sobol section)

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

// JS-escape a string so it can sit inside a template literal between
// backticks: \, `, and ${ all need escaping.  We escape just $ for
// simplicity (slightly over-escapes but works).
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

int main(int argc, char** argv) {
    const std::string template_path = "mc_viewer_template.html";
    std::string out_csv   = (argc >= 2) ? argv[1] : "monte_carlo_outputs.csv";
    std::string in_csv    = (argc >= 3) ? argv[2] : "monte_carlo_inputs.csv";
    std::string output    = (argc >= 4) ? argv[3] : "mc_viewer.html";
    std::string sobol_csv = (argc >= 5) ? argv[4] : "sobol_indices.csv";

    // Pairs CSV is derived from sobol_csv by inserting "_pairs" before
    // the extension.  Matches the sobol_runner naming convention.
    std::string pairs_csv, triplets_csv;
    {
        size_t dot = sobol_csv.find_last_of('.');
        if (dot != std::string::npos) {
            pairs_csv    = sobol_csv.substr(0, dot) + "_pairs"    + sobol_csv.substr(dot);
            triplets_csv = sobol_csv.substr(0, dot) + "_triplets" + sobol_csv.substr(dot);
        } else {
            pairs_csv    = sobol_csv + "_pairs";
            triplets_csv = sobol_csv + "_triplets";
        }
    }

    if (!file_exists(template_path)) {
        std::fprintf(stderr, "error: %s not found\n", template_path.c_str());
        return 1;
    }
    if (!file_exists(out_csv)) {
        std::fprintf(stderr, "error: %s not found\n", out_csv.c_str());
        return 1;
    }
    if (!file_exists(in_csv)) {
        std::fprintf(stderr, "error: %s not found\n", in_csv.c_str());
        return 1;
    }

    std::string templ, out_text_raw, in_text_raw, sobol_text_raw, pairs_text_raw, triplets_text_raw;
    if (!read_file(template_path, templ))   { std::fprintf(stderr, "error reading %s\n", template_path.c_str()); return 1; }
    if (!read_file(out_csv,        out_text_raw)) { std::fprintf(stderr, "error reading %s\n", out_csv.c_str()); return 1; }
    if (!read_file(in_csv,         in_text_raw))  { std::fprintf(stderr, "error reading %s\n", in_csv.c_str());  return 1; }

    bool have_sobol = file_exists(sobol_csv);
    if (have_sobol) {
        if (!read_file(sobol_csv, sobol_text_raw)) {
            std::fprintf(stderr, "warn: failed to read %s; skipping Sobol section\n", sobol_csv.c_str());
            have_sobol = false;
        }
    }

    bool have_pairs = file_exists(pairs_csv);
    if (have_pairs) {
        if (!read_file(pairs_csv, pairs_text_raw)) {
            std::fprintf(stderr, "warn: failed to read %s; skipping pairs section\n", pairs_csv.c_str());
            have_pairs = false;
        }
    }

    bool have_triplets = file_exists(triplets_csv);
    if (have_triplets) {
        if (!read_file(triplets_csv, triplets_text_raw)) {
            std::fprintf(stderr, "warn: failed to read %s; skipping triplets section\n", triplets_csv.c_str());
            have_triplets = false;
        }
    }

    // Optional targets sidecar (emitted by the runner when
    // target_ci_width is set in the config).  Derived from the indices
    // CSV path: sobol_indices.csv -> sobol_indices_targets.json.
    std::string targets_path;
    {
        size_t dot = sobol_csv.find_last_of('.');
        if (dot != std::string::npos)
            targets_path = sobol_csv.substr(0, dot) + "_targets.json";
        else
            targets_path = sobol_csv + "_targets.json";
    }
    std::string targets_text_raw;
    bool have_targets = file_exists(targets_path);
    if (have_targets) {
        if (!read_file(targets_path, targets_text_raw)) {
            std::fprintf(stderr, "warn: failed to read %s; skipping targets section\n",
                         targets_path.c_str());
            have_targets = false;
        }
    }

    std::string out_text = escape_for_js(out_text_raw);
    std::string in_text  = escape_for_js(in_text_raw);
    std::string sobol_text, pairs_text, triplets_text, targets_text;
    if (have_sobol)    sobol_text    = escape_for_js(sobol_text_raw);
    if (have_pairs)    pairs_text    = escape_for_js(pairs_text_raw);
    if (have_triplets) triplets_text = escape_for_js(triplets_text_raw);
    if (have_targets)  targets_text  = escape_for_js(targets_text_raw);

    // Build the <script> injection block
    std::ostringstream injection;
    injection << "<script>\n";
    injection << "window.MC_OUTPUTS_CSV = `" << out_text << "`;\n";
    injection << "window.MC_INPUTS_CSV  = `" << in_text  << "`;\n";
    if (have_sobol)    injection << "window.MC_SOBOL_CSV         = `" << sobol_text    << "`;\n";
    if (have_pairs)    injection << "window.MC_SOBOL_PAIRS_CSV    = `" << pairs_text    << "`;\n";
    if (have_triplets) injection << "window.MC_SOBOL_TRIPLETS_CSV = `" << triplets_text << "`;\n";
    if (have_targets)  injection << "window.MC_SOBOL_TARGETS_JSON = `" << targets_text  << "`;\n";
    injection << "</script>\n";

    // Substitute the marker
    const std::string marker = "<!-- CSV_DATA_INJECTION_POINT -->";
    size_t pos = templ.find(marker);
    if (pos == std::string::npos) {
        std::fprintf(stderr, "error: marker '%s' not found in template\n", marker.c_str());
        return 1;
    }
    templ.replace(pos, marker.size(), injection.str());

    // Write output
    std::ofstream fout(output, std::ios::binary);
    if (!fout) {
        std::fprintf(stderr, "error: cannot open %s for writing\n", output.c_str());
        return 1;
    }
    fout.write(templ.data(), static_cast<std::streamsize>(templ.size()));
    fout.close();

    std::printf("Wrote %s (%zu bytes)\n", output.c_str(), templ.size());
    std::printf("  outputs from %s (%zu bytes)\n", out_csv.c_str(), out_text.size());
    std::printf("  inputs  from %s (%zu bytes)\n", in_csv.c_str(),  in_text.size());
    if (have_sobol) {
        std::printf("  sobol   from %s (%zu bytes)\n", sobol_csv.c_str(), sobol_text.size());
    } else {
        std::printf("  sobol   not found at %s (section will be hidden)\n", sobol_csv.c_str());
    }
    if (have_pairs) {
        std::printf("  pairs   from %s (%zu bytes)\n", pairs_csv.c_str(), pairs_text.size());
    } else {
        std::printf("  pairs   not found at %s (section will be hidden)\n", pairs_csv.c_str());
    }
    if (have_triplets) {
        std::printf("  triplets from %s (%zu bytes)\n", triplets_csv.c_str(), triplets_text.size());
    }
    if (have_targets) {
        std::printf("  targets from %s (%zu bytes)\n", targets_path.c_str(), targets_text.size());
    }
    return 0;
}
