//  aero_deck.cpp  --  Aerodynamic table loader implementation
//
//  Two parsers:
//    - parse_osk_*  : OSK Table1/Table2 format (the legacy aero.txt)
//    - parse_datcom_*: Missile DATCOM .asc format (Zipfel-compatible)
//
//  Both use a shared lookahead-tokenizer that handles `//` line
//  comments and skips free-form text.  The tokenizer normalizes line
//  comments by simply truncating at "//".

#include "aero_deck.h"

#include <algorithm>
#include <cctype>
#include <cstdio>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>

namespace rocket6dof {

namespace {

// ---- Format detection ----
bool ends_with_ci(const std::string& s, const std::string& suffix) {
    if (s.size() < suffix.size()) return false;
    for (size_t i = 0; i < suffix.size(); ++i) {
        char a = s[s.size() - suffix.size() + i];
        char b = suffix[i];
        if (std::tolower(static_cast<unsigned char>(a)) !=
            std::tolower(static_cast<unsigned char>(b))) return false;
    }
    return true;
}

// Strip everything from "//" to end-of-line.
std::string strip_line_comment(const std::string& line) {
    auto pos = line.find("//");
    if (pos == std::string::npos) return line;
    return line.substr(0, pos);
}

// Read the whole file with line comments stripped.  Returns true on
// success.  We don't try to be efficient -- aero decks are small.
bool slurp_file(const std::string& path, std::string& out) {
    std::ifstream f(path);
    if (!f) {
        std::cerr << "AeroDeck: cannot open '" << path << "'\n";
        return false;
    }
    std::ostringstream o;
    std::string line;
    while (std::getline(f, line)) {
        o << strip_line_comment(line) << '\n';
    }
    out = o.str();
    return true;
}

bool is_number_token(const std::string& tok) {
    if (tok.empty()) return false;
    size_t i = 0;
    if (tok[0] == '+' || tok[0] == '-') i = 1;
    bool seen_digit = false;
    for (; i < tok.size(); ++i) {
        char c = tok[i];
        if (std::isdigit(static_cast<unsigned char>(c))) seen_digit = true;
        else if (c == '.' || c == 'e' || c == 'E' || c == '+' || c == '-') {
            continue;
        } else {
            return false;
        }
    }
    return seen_digit;
}

// ---- OSK format parser ----
// After locating the tag line, read:
//   <n> ints (for 1D: 1 size; for 2D: 2 sizes)
//   <label1>
//   axis1 values
//   [for 2D: <label2>, axis2 values, <label3>]
//   data values
//
// Labels are arbitrary non-numeric tokens.  We consume them by reading
// a single token and discarding if not a number.

bool find_tag(std::istringstream& iss, const std::string& tag) {
    std::string line;
    iss.clear();   // reset
    iss.seekg(0);
    while (std::getline(iss, line)) {
        std::istringstream ls(line);
        std::string first;
        if (ls >> first && first == tag) return true;
    }
    return false;
}

// After find_tag, skip non-numeric "label" lines until we hit data.
// We do this by peeking at the next token.
bool skip_labels(std::istringstream& iss) {
    std::streampos pos = iss.tellg();
    std::string tok;
    while (iss >> tok) {
        if (is_number_token(tok)) {
            // Rewind one token
            iss.clear();
            iss.seekg(pos);
            return true;
        }
        pos = iss.tellg();
    }
    return false;
}

bool parse_osk_1d(const std::string& body, const std::string& tag,
                  std::vector<double>& x_out,
                  std::vector<double>& y_out)
{
    std::istringstream iss(body);
    if (!find_tag(iss, tag)) return false;
    int n = 0;
    if (!(iss >> n) || n <= 0) return false;
    if (!skip_labels(iss)) return false;
    x_out.resize(n);
    y_out.resize(n);
    for (int i = 0; i < n; ++i) {
        if (!(iss >> x_out[i] >> y_out[i])) return false;
    }
    return true;
}

bool parse_osk_2d(const std::string& body, const std::string& tag,
                  std::vector<double>& x1_out,
                  std::vector<double>& x2_out,
                  std::vector<double>& z_out)
{
    std::istringstream iss(body);
    if (!find_tag(iss, tag)) return false;
    int n1 = 0, n2 = 0;
    if (!(iss >> n1 >> n2) || n1 <= 0 || n2 <= 0) return false;
    if (!skip_labels(iss)) return false;
    x1_out.resize(n1);
    for (int i = 0; i < n1; ++i) {
        if (!(iss >> x1_out[i])) return false;
    }
    if (!skip_labels(iss)) return false;
    x2_out.resize(n2);
    for (int j = 0; j < n2; ++j) {
        if (!(iss >> x2_out[j])) return false;
    }
    if (!skip_labels(iss)) return false;
    z_out.resize(static_cast<size_t>(n1) * n2);
    for (int i = 0; i < n1; ++i) {
        for (int j = 0; j < n2; ++j) {
            if (!(iss >> z_out[i*n2 + j])) return false;
        }
    }
    return true;
}

// ---- DATCOM format parser ----
//
// Header is one of:
//   1DIM tag
//   2DIM tag
// followed by:
//   NX1 n
//     or
//   NX1 n1 NX2 n2
//
// For 1D data: pairs of (x, y) until next header or EOF.
//
// For 2D data: lines of (x1_value, [x2_value,]  z_0, z_1, ..., z_{n2-1}).
// The x2 value column is filled only for the first n2 rows; subsequent
// rows have only x1 and the n2 data columns.  We disambiguate by
// counting tokens on each row.

// Find the line that begins "1DIM tag" or "2DIM tag".  Returns the
// dimension (1 or 2) and stream positioned just after the header.
int find_datcom_tag(std::istringstream& iss, const std::string& tag) {
    iss.clear();
    iss.seekg(0);
    std::string line;
    while (std::getline(iss, line)) {
        std::istringstream ls(line);
        std::string marker, name;
        if (ls >> marker >> name) {
            if (name == tag) {
                if (marker == "1DIM") return 1;
                if (marker == "2DIM") return 2;
            }
        }
    }
    return 0;
}

// Read tokens, skipping any non-number tokens (NX1, NX2, alpha labels
// like "mach", etc).  Stops at end of stream.  Doesn't help us locate
// section boundaries -- we use a different approach: read until we've
// consumed exactly the expected number of numbers.
bool next_number(std::istringstream& iss, double& out) {
    std::string tok;
    while (iss >> tok) {
        if (is_number_token(tok)) {
            out = std::strtod(tok.c_str(), nullptr);
            return true;
        }
        // Not a number; check if it's a stop marker (next section)
        if (tok == "1DIM" || tok == "2DIM") {
            // We've run past our section; rewind by length of token
            // so the caller can handle EOF
            iss.clear();
            // Push back: easiest just to leave the stream at EOF for
            // this section and return false
            iss.setstate(std::ios::eofbit);
            return false;
        }
        // Else: just skip the token (could be NX1, NX2, axis name, etc)
    }
    return false;
}

// Read n_x and n_y from the NX1/NX2 markers.  We assume they appear
// in order: NX1 followed by its number, optionally NX2 followed by
// its number.
bool read_datcom_sizes(std::istringstream& iss, int& n1, int& n2) {
    n1 = 0; n2 = 0;
    std::string tok;
    std::streampos save = iss.tellg();
    while (iss >> tok) {
        if (tok == "NX1") {
            if (!(iss >> n1) || n1 <= 0) return false;
        } else if (tok == "NX2") {
            if (!(iss >> n2) || n2 <= 0) return false;
        } else if (tok == "1DIM" || tok == "2DIM") {
            // Started next section: rewind to just before this token
            iss.clear();
            iss.seekg(save);
            break;
        } else if (is_number_token(tok)) {
            // First data value -- rewind and stop
            iss.clear();
            iss.seekg(save);
            break;
        }
        save = iss.tellg();
    }
    return n1 > 0;
}

bool parse_datcom_1d(const std::string& body, const std::string& tag,
                     std::vector<double>& x_out,
                     std::vector<double>& y_out)
{
    std::istringstream iss(body);
    int dim = find_datcom_tag(iss, tag);
    if (dim != 1) return false;

    int n1 = 0, n2 = 0;
    if (!read_datcom_sizes(iss, n1, n2) || n1 <= 0) return false;

    x_out.resize(n1);
    y_out.resize(n1);
    for (int i = 0; i < n1; ++i) {
        if (!next_number(iss, x_out[i])) return false;
        if (!next_number(iss, y_out[i])) return false;
    }
    return true;
}

bool parse_datcom_2d(const std::string& body, const std::string& tag,
                     std::vector<double>& x1_out,
                     std::vector<double>& x2_out,
                     std::vector<double>& z_out)
{
    // Read the file's lines starting after the tag.  For each row of
    // data we read the row's tokens, examine the token count, and
    // decide whether row[1] is X2 or part of the data.

    // First locate the tag line and consume the size header.
    int dim;
    int n1 = 0, n2 = 0;

    // Find tag line and remember position.
    std::istringstream iss_scan(body);
    {
        std::string line;
        bool found = false;
        size_t after_tag_offset = 0;
        while (std::getline(iss_scan, line)) {
            std::istringstream ls(line);
            std::string marker, name;
            if (ls >> marker >> name && name == tag) {
                if (marker == "2DIM") {
                    found = true;
                    after_tag_offset = iss_scan.tellg();
                    break;
                }
                return false;
            }
        }
        if (!found) return false;
        // Reset to position right after the tag line
        iss_scan.clear();
        iss_scan.seekg(after_tag_offset);
        dim = 2;
        (void)dim;
    }

    // Read NX1, NX2 from the next non-blank line(s)
    {
        std::string line;
        while (std::getline(iss_scan, line)) {
            std::istringstream ls(line);
            std::string tok;
            while (ls >> tok) {
                if (tok == "NX1") {
                    if (!(ls >> n1) || n1 <= 0) return false;
                } else if (tok == "NX2") {
                    if (!(ls >> n2) || n2 <= 0) return false;
                }
            }
            if (n1 > 0 && n2 > 0) break;
        }
        if (n1 <= 0 || n2 <= 0) return false;
    }

    x1_out.assign(n1, 0.0);
    x2_out.assign(n2, 0.0);
    z_out.assign(static_cast<size_t>(n1) * n2, 0.0);

    // Now read n1 data rows.  For row i:
    //   if i < n2:    n1_val  x2_val  z[i][0..n2-1]   (2 + n2 tokens)
    //   else:         n1_val          z[i][0..n2-1]   (1 + n2 tokens)
    // We stop reading on the next "1DIM"/"2DIM" marker or EOF.
    int rows_read = 0;
    while (rows_read < n1) {
        std::string line;
        if (!std::getline(iss_scan, line)) break;
        // Strip leading whitespace; skip empty lines
        std::istringstream ls(line);
        std::vector<double> nums;
        std::string tok;
        bool saw_marker = false;
        while (ls >> tok) {
            if (tok == "1DIM" || tok == "2DIM") { saw_marker = true; break; }
            if (is_number_token(tok)) nums.push_back(std::strtod(tok.c_str(), nullptr));
        }
        if (saw_marker) break;
        if (nums.empty()) continue;     // blank / comment-only line

        // Decide whether this row has X2 in column 1
        int expected_with_x2    = 2 + n2;
        int expected_without_x2 = 1 + n2;

        if (rows_read < n2 && static_cast<int>(nums.size()) == expected_with_x2) {
            x1_out[rows_read] = nums[0];
            x2_out[rows_read] = nums[1];
            for (int j = 0; j < n2; ++j) {
                z_out[rows_read * n2 + j] = nums[2 + j];
            }
        }
        else if (static_cast<int>(nums.size()) == expected_without_x2) {
            x1_out[rows_read] = nums[0];
            for (int j = 0; j < n2; ++j) {
                z_out[rows_read * n2 + j] = nums[1 + j];
            }
        }
        else {
            std::cerr << "AeroDeck DATCOM 2D: unexpected row width "
                      << nums.size() << " for tag '" << tag
                      << "' row " << rows_read
                      << " (expected " << expected_with_x2
                      << " or " << expected_without_x2 << ")\n";
            return false;
        }
        ++rows_read;
    }

    return rows_read == n1;
}

// ---- Interpolation helpers (shared with osk::Table internals) ----
void locate(const std::vector<double>& x, double xq,
            size_t& i, double& t, bool& clamped)
{
    clamped = false;
    if (x.size() < 2) { i = 0; t = 0.0; clamped = true; return; }
    if (xq <= x.front()) { i = 0;          t = 0.0; clamped = true; return; }
    if (xq >= x.back())  { i = x.size()-1; t = 0.0; clamped = true; return; }
    auto it = std::upper_bound(x.begin(), x.end(), xq);
    i = static_cast<size_t>(std::distance(x.begin(), it) - 1);
    t = (xq - x[i]) / (x[i+1] - x[i]);
}

} // anon

AeroFormat detect_aero_format(const std::string& path) {
    return ends_with_ci(path, ".asc") ? AERO_FORMAT_DATCOM : AERO_FORMAT_OSK;
}

// ---- AeroTable1 ----
AeroTable1::AeroTable1(const std::string& path) : path_(path), loaded_(false) {}

bool AeroTable1::read(const std::string& tag) {
    tag_ = tag;
    std::string body;
    if (!slurp_file(path_, body)) return false;

    bool ok = false;
    if (detect_aero_format(path_) == AERO_FORMAT_DATCOM) {
        ok = parse_datcom_1d(body, tag, x_, y_);
    } else {
        ok = parse_osk_1d(body, tag, x_, y_);
    }
    if (!ok) {
        std::cerr << "AeroTable1: failed to load tag '" << tag
                  << "' from '" << path_ << "'\n";
    }
    loaded_ = ok;
    return ok;
}

double AeroTable1::interp(double x) const {
    if (x_.empty()) return 0.0;
    size_t i; double t; bool clamped;
    locate(x_, x, i, t, clamped);
    if (clamped) return y_[i];
    return y_[i] + t * (y_[i+1] - y_[i]);
}

// ---- AeroTable2 ----
AeroTable2::AeroTable2(const std::string& path) : path_(path), loaded_(false) {}

bool AeroTable2::read(const std::string& tag) {
    tag_ = tag;
    std::string body;
    if (!slurp_file(path_, body)) return false;

    bool ok = false;
    if (detect_aero_format(path_) == AERO_FORMAT_DATCOM) {
        ok = parse_datcom_2d(body, tag, x1_, x2_, z_);
    } else {
        ok = parse_osk_2d(body, tag, x1_, x2_, z_);
    }
    if (!ok) {
        std::cerr << "AeroTable2: failed to load tag '" << tag
                  << "' from '" << path_ << "'\n";
    }
    loaded_ = ok;
    return ok;
}

double AeroTable2::interp(double x1, double x2) const {
    if (x1_.empty() || x2_.empty()) return 0.0;
    size_t i, j; double ti, tj; bool ci, cj;
    locate(x1_, x1, i, ti, ci);
    locate(x2_, x2, j, tj, cj);
    const int N2 = n2();

    auto Z = [&](size_t ii, size_t jj) {
        return z_[ii * static_cast<size_t>(N2) + jj];
    };
    size_t i1 = ci ? i : i + 1;
    size_t j1 = cj ? j : j + 1;

    double z00 = Z(i,  j );
    double z01 = Z(i,  j1);
    double z10 = Z(i1, j );
    double z11 = Z(i1, j1);

    double z0 = z00 + tj * (z01 - z00);
    double z1 = z10 + tj * (z11 - z10);
    return z0 + ti * (z1 - z0);
}

} // namespace rocket6dof
