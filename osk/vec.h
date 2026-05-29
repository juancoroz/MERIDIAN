
#ifndef OSK_VEC_H
#define OSK_VEC_H

#include <iosfwd>

namespace osk {

constexpr double PI = 3.14159265358979323846;

class Mat;
class Quat;

class Vec {
public:

    double x, y, z;
    double m;

    Vec();
    Vec(double x_, double y_, double z_);

    double&       operator[](int i);
    const double& operator[](int i) const;

    Vec& operator()(double x_, double y_, double z_);

    Vec& operator=(double s);

    void extract(double& a, double& b, double& c) const;

    Vec operator+(const Vec& other) const;
    Vec operator-(const Vec& other) const;
    Vec operator*(double s) const;
    Vec operator/(double s) const;

    Vec& operator*=(double s);
    Vec& operator/=(double s);

    Vec    scale(double s) const;
    double mag() const;
    Vec    unit() const;
    Vec    apply(double (*f)(double)) const;

    double dot(const Vec& other) const;
    Vec    cross(const Vec& other) const;

    Mat  getDCM()  const;
    Quat getQuat() const;
};

std::ostream& operator<<(std::ostream& os, const Vec& v);

Vec operator*(double s, const Vec& v);

}

#endif
