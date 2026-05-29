
#include "vec.h"
#include "mat.h"
#include "quat.h"
#include <cmath>
#include <iostream>
#include <iomanip>

namespace osk {

Vec::Vec() : x(0.0), y(0.0), z(0.0), m(0.0) {}
Vec::Vec(double x_, double y_, double z_)
    : x(x_), y(y_), z(z_), m(0.0) {}

double& Vec::operator[](int i) {
    switch (i) {
        case 0: return x;
        case 1: return y;
        case 2: return z;
        default: return z;
    }
}
const double& Vec::operator[](int i) const {
    switch (i) {
        case 0: return x;
        case 1: return y;
        case 2: return z;
        default: return z;
    }
}

Vec& Vec::operator()(double x_, double y_, double z_) {
    x = x_; y = y_; z = z_;
    return *this;
}

Vec& Vec::operator=(double s) {
    x = y = z = s;
    return *this;
}

void Vec::extract(double& a, double& b, double& c) const {
    a = x; b = y; c = z;
}

Vec Vec::operator+(const Vec& o) const { return Vec(x + o.x, y + o.y, z + o.z); }
Vec Vec::operator-(const Vec& o) const { return Vec(x - o.x, y - o.y, z - o.z); }
Vec Vec::operator*(double s)     const { return Vec(x * s, y * s, z * s); }
Vec Vec::operator/(double s)     const { return Vec(x / s, y / s, z / s); }

Vec& Vec::operator*=(double s) { x *= s; y *= s; z *= s; return *this; }
Vec& Vec::operator/=(double s) { x /= s; y /= s; z /= s; return *this; }

Vec    Vec::scale(double s) const { return Vec(x * s, y * s, z * s); }
double Vec::mag() const            { return std::sqrt(x*x + y*y + z*z); }

Vec Vec::unit() const {
    double M = mag();
    if (M == 0.0) return Vec(0.0, 0.0, 0.0);
    return Vec(x / M, y / M, z / M);
}

Vec Vec::apply(double (*f)(double)) const {
    return Vec(f(x), f(y), f(z));
}

double Vec::dot(const Vec& o) const {
    return x * o.x + y * o.y + z * o.z;
}

Vec Vec::cross(const Vec& o) const {
    return Vec(y * o.z - z * o.y,
               z * o.x - x * o.z,
               x * o.y - y * o.x);
}

Mat Vec::getDCM() const {
    double phi = x, theta = y, psi = z;
    double cphi = std::cos(phi),   sphi = std::sin(phi);
    double cth  = std::cos(theta), sth  = std::sin(theta);
    double cpsi = std::cos(psi),   spsi = std::sin(psi);
    return Mat(
        cth*cpsi,                 cth*spsi,                -sth,
        sphi*sth*cpsi - cphi*spsi,sphi*sth*spsi + cphi*cpsi,sphi*cth,
        cphi*sth*cpsi + sphi*spsi,cphi*sth*spsi - sphi*cpsi,cphi*cth
    );
}

Quat Vec::getQuat() const {
    double cphi = std::cos(x / 2.0), sphi = std::sin(x / 2.0);
    double cth  = std::cos(y / 2.0), sth  = std::sin(y / 2.0);
    double cpsi = std::cos(z / 2.0), spsi = std::sin(z / 2.0);
    return Quat(
        cphi*cth*cpsi + sphi*sth*spsi,
        sphi*cth*cpsi - cphi*sth*spsi,
        cphi*sth*cpsi + sphi*cth*spsi,
        cphi*cth*spsi - sphi*sth*cpsi
    );
}

std::ostream& operator<<(std::ostream& os, const Vec& v) {
    os << "(" << std::setw(10) << v.x << ", "
              << std::setw(10) << v.y << ", "
              << std::setw(10) << v.z << ")";
    return os;
}

Vec operator*(double s, const Vec& v) { return v * s; }

}
