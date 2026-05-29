//  aero_deck.h  --  Aerodynamic table loader, OSK + DATCOM formats
//
//  Holds 1-D and 2-D interpolation tables loaded from one of two file
//  formats:
//
//    1. The original OSK Table1/Table2 format used by aero.txt:
//
//         tag_name
//         <n>
//         label
//          x1  y1
//          x2  y2
//          ...
//
//       and for 2-D:
//
//         tag_name
//         <n1> <n2>
//         label_x1
//           x1[0] x1[1] ... x1[n1-1]
//         label_x2
//           x2[0] x2[1] ... x2[n2-1]
//         label_z
//           z[0][0] z[0][1] ... z[n1-1][n2-1]
//
//    2. Missile DATCOM .asc format used by Zipfel's aero_deck_SLV.asc:
//
//         1DIM tag_name
//         NX1 <n>  // optional comment
//          x1  y1
//          x2  y2
//          ...
//
//       and for 2-D:
//
//         2DIM tag_name
//         NX1 <n1> NX2 <n2>  // optional comment
//          x1[0]  x2[0]   z[0][0]  z[0][1]  ...  z[0][n2-1]
//          x1[1]  x2[1]   z[1][0]  z[1][1]  ...  z[1][n2-1]
//          ...
//          x1[n2-1] x2[n2-1] z[n2-1][0] ...
//          x1[n2]            z[n2][0]   ...  (x2 column blank past row n2-1)
//          ...
//          x1[n1-1]          z[n1-1][0] ...
//
//       The first NX2 rows include the X2 axis values in column 2; the
//       remaining rows have only X1 plus data.
//
//  The file extension determines the format: `.asc` -> DATCOM, anything
//  else -> OSK Table.  Both formats tolerate `//` line comments and
//  free-form text outside data blocks.
//
//  Once loaded, an AeroTable1 or AeroTable2 supports the same interp()
//  operator() interface as osk::Table1 / osk::Table2, with the same
//  linear interpolation + clamp-at-boundary semantics.

#ifndef ROCKET6DOF_AERO_DECK_H
#define ROCKET6DOF_AERO_DECK_H

#include <string>
#include <vector>

namespace rocket6dof {

enum AeroFormat {
    AERO_FORMAT_OSK,     // tag\n n\n label\n x y x y ...
    AERO_FORMAT_DATCOM   // 1DIM tag, NX1 n, ...
};

// Detect the format from the file extension.  Default OSK.
AeroFormat detect_aero_format(const std::string& path);

// 1-D table: y(x).  Loaded by tag name from a file in either format.
class AeroTable1 {
public:
    explicit AeroTable1(const std::string& path);

    bool read(const std::string& tag);

    double interp(double x) const;
    double operator()(double x) const { return interp(x); }

    int size() const { return static_cast<int>(x_.size()); }
    bool loaded() const { return loaded_; }

private:
    std::string path_;
    std::string tag_;
    std::vector<double> x_, y_;
    bool loaded_;
};

// 2-D table: z(x1, x2).  Bilinear interpolation in the box, nearest-
// neighbour clamp outside.  Row-major: z_[i*n2 + j] is the value at
// (x1_[i], x2_[j]).
class AeroTable2 {
public:
    explicit AeroTable2(const std::string& path);

    bool read(const std::string& tag);

    double interp(double x1, double x2) const;
    double operator()(double x1, double x2) const { return interp(x1, x2); }

    int n1() const { return static_cast<int>(x1_.size()); }
    int n2() const { return static_cast<int>(x2_.size()); }
    bool loaded() const { return loaded_; }

private:
    std::string path_;
    std::string tag_;
    std::vector<double> x1_, x2_;
    std::vector<double> z_;
    bool loaded_;
};

} // namespace rocket6dof

#endif
