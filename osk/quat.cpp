
#include "quat.h"
#include "mat.h"
#include "vec.h"
#include <cmath>
#include <iostream>
#include <iomanip>

namespace osk {

Quat::Quat() : s(0.0), x(0.0), y(0.0), z(0.0) {}
Quat::Quat(double s_, double x_, double y_, double z_)
    : s(s_), x(x_), y(y_), z(z_) {}

double& Quat::operator[](int i) {
    switch (i) {
        case 0: return s;
        case 1: return x;
        case 2: return y;
        case 3: return z;
        default: return z;
    }
}
const double& Quat::operator[](int i) const {
    switch (i) {
        case 0: return s;
        case 1: return x;
        case 2: return y;
        case 3: return z;
        default: return z;
    }
}

Quat& Quat::operator()(double s_, double x_, double y_, double z_) {
    s = s_; x = x_; y = y_; z = z_;
    return *this;
}

double Quat::mag() const {
    return std::sqrt(s*s + x*x + y*y + z*z);
}

void Quat::normalize() {
    double M = mag();
    if (M == 0.0) return;
    s /= M; x /= M; y /= M; z /= M;
}

Mat Quat::getDCM() const {
    double xx = x*x, yy = y*y, zz = z*z;
    double xy = x*y, xz = x*z, yz = y*z;
    double sx = s*x, sy = s*y, sz = s*z;
    return Mat(
        1.0 - 2.0*(yy + zz),  2.0*(xy + sz),         2.0*(xz - sy),
        2.0*(xy - sz),        1.0 - 2.0*(xx + zz),   2.0*(yz + sx),
        2.0*(xz + sy),        2.0*(yz - sx),         1.0 - 2.0*(xx + yy)
    );
}

Vec Quat::getEuler() const {
    return getDCM().getEuler();
}

std::ostream& operator<<(std::ostream& os, const Quat& q) {
    os << "[" << std::setw(10) << q.s << "; "
              << std::setw(10) << q.x << ", "
              << std::setw(10) << q.y << ", "
              << std::setw(10) << q.z << "]";
    return os;
}

}
