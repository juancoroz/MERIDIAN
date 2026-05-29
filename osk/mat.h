
#ifndef OSK_MAT_H
#define OSK_MAT_H

#include "vec.h"
#include <iosfwd>

namespace osk {

class Quat;

class Mat {
public:

    Vec rows_[3];

    Mat();
    Mat(double a00, double a01, double a02,
        double a10, double a11, double a12,
        double a20, double a21, double a22);
    Mat(const Vec& r0, const Vec& r1, const Vec& r2);

    Vec&       operator[](int i);
    const Vec& operator[](int i) const;

    Mat& operator()(double a00, double a01, double a02,
                    double a10, double a11, double a12,
                    double a20, double a21, double a22);
    Mat& operator()(const Vec& r0, const Vec& r1, const Vec& r2);

    Mat operator+(const Mat& other) const;
    Mat operator-(const Mat& other) const;
    Mat operator*(const Mat& other) const;
    Vec operator*(const Vec& v) const;
    Mat operator*(double s) const;
    Mat operator/(double s) const;

    Mat& operator*=(double s);
    Mat& operator/=(double s);

    Mat    scale(double s) const;
    Mat    transpose() const;
    double det() const;
    Mat    inv() const;
    Mat    apply(double (*f)(double)) const;

    Vec  getEuler() const;
    Quat getQuat()  const;
};

std::ostream& operator<<(std::ostream& os, const Mat& m);

Mat operator*(double s, const Mat& m);

}

#endif
