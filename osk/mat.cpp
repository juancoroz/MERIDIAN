
#include "mat.h"
#include "quat.h"
#include <cmath>
#include <iostream>
#include <iomanip>
#include <cstdlib>

namespace osk {

Mat::Mat() {
    rows_[0] = Vec(); rows_[1] = Vec(); rows_[2] = Vec();
}
Mat::Mat(double a00, double a01, double a02,
         double a10, double a11, double a12,
         double a20, double a21, double a22) {
    rows_[0] = Vec(a00, a01, a02);
    rows_[1] = Vec(a10, a11, a12);
    rows_[2] = Vec(a20, a21, a22);
}
Mat::Mat(const Vec& r0, const Vec& r1, const Vec& r2) {
    rows_[0] = r0; rows_[1] = r1; rows_[2] = r2;
}

Vec&       Mat::operator[](int i)       {
    if (i < 0 || i > 2) i = 2;
    return rows_[i];
}
const Vec& Mat::operator[](int i) const {
    if (i < 0 || i > 2) i = 2;
    return rows_[i];
}

Mat& Mat::operator()(double a00, double a01, double a02,
                     double a10, double a11, double a12,
                     double a20, double a21, double a22) {
    rows_[0](a00, a01, a02);
    rows_[1](a10, a11, a12);
    rows_[2](a20, a21, a22);
    return *this;
}
Mat& Mat::operator()(const Vec& r0, const Vec& r1, const Vec& r2) {
    rows_[0] = r0; rows_[1] = r1; rows_[2] = r2;
    return *this;
}

Mat Mat::operator+(const Mat& o) const {
    return Mat(rows_[0] + o.rows_[0],
               rows_[1] + o.rows_[1],
               rows_[2] + o.rows_[2]);
}
Mat Mat::operator-(const Mat& o) const {
    return Mat(rows_[0] - o.rows_[0],
               rows_[1] - o.rows_[1],
               rows_[2] - o.rows_[2]);
}

Mat Mat::operator*(const Mat& o) const {
    Mat r;
    for (int i = 0; i < 3; ++i) {
        for (int j = 0; j < 3; ++j) {
            double s = 0.0;
            for (int k = 0; k < 3; ++k) s += rows_[i][k] * o.rows_[k][j];
            r.rows_[i][j] = s;
        }
    }
    return r;
}

Vec Mat::operator*(const Vec& v) const {
    return Vec(rows_[0].dot(v),
               rows_[1].dot(v),
               rows_[2].dot(v));
}

Mat Mat::operator*(double s) const {
    return Mat(rows_[0] * s, rows_[1] * s, rows_[2] * s);
}
Mat Mat::operator/(double s) const {
    return Mat(rows_[0] / s, rows_[1] / s, rows_[2] / s);
}

Mat& Mat::operator*=(double s) {
    rows_[0] *= s; rows_[1] *= s; rows_[2] *= s;
    return *this;
}
Mat& Mat::operator/=(double s) {
    rows_[0] /= s; rows_[1] /= s; rows_[2] /= s;
    return *this;
}

Mat Mat::scale(double s) const {
    return Mat(rows_[0] * s, rows_[1] * s, rows_[2] * s);
}

Mat Mat::transpose() const {
    return Mat(rows_[0][0], rows_[1][0], rows_[2][0],
               rows_[0][1], rows_[1][1], rows_[2][1],
               rows_[0][2], rows_[1][2], rows_[2][2]);
}

double Mat::det() const {
    const Vec& r0 = rows_[0]; const Vec& r1 = rows_[1]; const Vec& r2 = rows_[2];
    return r0[0] * (r1[1] * r2[2] - r1[2] * r2[1])
         - r0[1] * (r1[0] * r2[2] - r1[2] * r2[0])
         + r0[2] * (r1[0] * r2[1] - r1[1] * r2[0]);
}

Mat Mat::inv() const {
    double D = det();
    if (std::fabs(D) < 1.0e-300) {
        std::cerr << "Mat::inv: singular matrix\n";
        return *this;
    }
    const Vec& r0 = rows_[0]; const Vec& r1 = rows_[1]; const Vec& r2 = rows_[2];

    Mat adj(
        (r1[1]*r2[2] - r1[2]*r2[1]), -(r0[1]*r2[2] - r0[2]*r2[1]),  (r0[1]*r1[2] - r0[2]*r1[1]),
       -(r1[0]*r2[2] - r1[2]*r2[0]),  (r0[0]*r2[2] - r0[2]*r2[0]), -(r0[0]*r1[2] - r0[2]*r1[0]),
        (r1[0]*r2[1] - r1[1]*r2[0]), -(r0[0]*r2[1] - r0[1]*r2[0]),  (r0[0]*r1[1] - r0[1]*r1[0])
    );
    return adj / D;
}

Mat Mat::apply(double (*f)(double)) const {
    return Mat(rows_[0].apply(f),
               rows_[1].apply(f),
               rows_[2].apply(f));
}

Vec Mat::getEuler() const {
    const double s2 = -rows_[0][2];
    double theta, phi, psi;
    if (std::fabs(s2) >= 1.0 - 1.0e-12) {

        theta = (s2 > 0.0 ? 1.0 : -1.0) * PI / 2.0;
        psi   = 0.0;
        phi   = std::atan2(rows_[1][0], rows_[1][1]);
    } else {
        theta = std::asin(s2);
        phi   = std::atan2(rows_[1][2], rows_[2][2]);
        psi   = std::atan2(rows_[0][1], rows_[0][0]);
    }
    return Vec(phi, theta, psi);
}

Quat Mat::getQuat() const {
    const Mat& M = *this;
    double tr = M[0][0] + M[1][1] + M[2][2];

    double s, x, y, z;
    if (tr > 0.0) {

        double S = std::sqrt(tr + 1.0) * 2.0;
        s = 0.25 * S;
        x = (M[1][2] - M[2][1]) / S;
        y = (M[2][0] - M[0][2]) / S;
        z = (M[0][1] - M[1][0]) / S;
    } else if (M[0][0] > M[1][1] && M[0][0] > M[2][2]) {
        double S = std::sqrt(1.0 + M[0][0] - M[1][1] - M[2][2]) * 2.0;
        s = (M[1][2] - M[2][1]) / S;
        x = 0.25 * S;
        y = (M[1][0] + M[0][1]) / S;
        z = (M[2][0] + M[0][2]) / S;
    } else if (M[1][1] > M[2][2]) {
        double S = std::sqrt(1.0 + M[1][1] - M[0][0] - M[2][2]) * 2.0;
        s = (M[2][0] - M[0][2]) / S;
        x = (M[1][0] + M[0][1]) / S;
        y = 0.25 * S;
        z = (M[2][1] + M[1][2]) / S;
    } else {
        double S = std::sqrt(1.0 + M[2][2] - M[0][0] - M[1][1]) * 2.0;
        s = (M[0][1] - M[1][0]) / S;
        x = (M[2][0] + M[0][2]) / S;
        y = (M[2][1] + M[1][2]) / S;
        z = 0.25 * S;
    }
    return Quat(s, x, y, z);
}

std::ostream& operator<<(std::ostream& os, const Mat& m) {
    os << "[" << m[0] << ",\n"
       << " " << m[1] << ",\n"
       << " " << m[2] << "]";
    return os;
}

Mat operator*(double s, const Mat& m) { return m * s; }

}
