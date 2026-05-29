
#ifndef OSK_QUAT_H
#define OSK_QUAT_H

#include "vec.h"
#include "mat.h"
#include <iosfwd>

namespace osk {

class Quat {
public:

    double s, x, y, z;

    Quat();
    Quat(double s_, double x_, double y_, double z_);

    double&       operator[](int i);
    const double& operator[](int i) const;

    Quat& operator()(double s_, double x_, double y_, double z_);

    void normalize();

    double mag() const;

    Mat getDCM()   const;

    Vec getEuler() const;
};

std::ostream& operator<<(std::ostream& os, const Quat& q);

}

#endif
