//  kalman_matrix.h  --  Small dense-matrix utility for Kalman filters
//
//  The OSK kernel's osk::Mat is fixed 3x3.  Zipfel's GPS Kalman filter
//  uses 8x8 matrices.  This header provides a minimal dense matrix
//  class with the operations needed: construct, set/get, transpose,
//  matrix multiply, matrix-vector multiply, addition/subtraction,
//  identity, and inverse (Gauss-Jordan with partial pivoting).
//
//  All operations are pre-allocated (rows*cols stored in std::vector),
//  no aliasing-free guarantees: caller responsible for not aliasing
//  output with input where the routine reads in-place.

#ifndef ROCKET6DOF_KALMAN_MATRIX_H
#define ROCKET6DOF_KALMAN_MATRIX_H

#include <vector>
#include <cstddef>
#include <cmath>
#include <stdexcept>

namespace rocket6dof {

class MatN {
public:
    std::size_t rows, cols;
    std::vector<double> data;   // row-major

    MatN() : rows(0), cols(0) {}
    MatN(std::size_t r, std::size_t c)
        : rows(r), cols(c), data(r * c, 0.0) {}

    double  at(std::size_t i, std::size_t j) const { return data[i*cols + j]; }
    double& at(std::size_t i, std::size_t j)       { return data[i*cols + j]; }

    void set(std::size_t i, std::size_t j, double v) { data[i*cols + j] = v; }

    void zero() {
        std::fill(data.begin(), data.end(), 0.0);
    }

    static MatN identity(std::size_t n) {
        MatN m(n, n);
        for (std::size_t i = 0; i < n; i++) m.at(i, i) = 1.0;
        return m;
    }

    MatN transpose() const {
        MatN r(cols, rows);
        for (std::size_t i = 0; i < rows; i++)
            for (std::size_t j = 0; j < cols; j++)
                r.at(j, i) = at(i, j);
        return r;
    }

    MatN operator+(const MatN& o) const {
        MatN r(rows, cols);
        for (std::size_t i = 0; i < data.size(); i++) r.data[i] = data[i] + o.data[i];
        return r;
    }
    MatN operator-(const MatN& o) const {
        MatN r(rows, cols);
        for (std::size_t i = 0; i < data.size(); i++) r.data[i] = data[i] - o.data[i];
        return r;
    }
    MatN operator*(double s) const {
        MatN r(rows, cols);
        for (std::size_t i = 0; i < data.size(); i++) r.data[i] = data[i] * s;
        return r;
    }

    // Matrix-matrix product: (rows x inner) * (inner x cols) -> (rows x cols)
    MatN operator*(const MatN& o) const {
        MatN r(rows, o.cols);
        for (std::size_t i = 0; i < rows; i++) {
            for (std::size_t k = 0; k < cols; k++) {
                double a = at(i, k);
                if (a == 0.0) continue;
                for (std::size_t j = 0; j < o.cols; j++) {
                    r.at(i, j) += a * o.at(k, j);
                }
            }
        }
        return r;
    }

    // Inverse via Gauss-Jordan elimination with partial pivoting.
    // For our 8x8 Kalman PP matrix this is well-conditioned (positive
    // definite covariance), so partial pivoting is sufficient.
    MatN inverse() const {
        if (rows != cols) {
            throw std::runtime_error("MatN::inverse on non-square matrix");
        }
        std::size_t n = rows;
        MatN A = *this;
        MatN I = identity(n);
        for (std::size_t k = 0; k < n; k++) {
            // Partial pivot: find largest |A(i,k)| in column k, swap rows
            std::size_t piv = k;
            double maxv = std::fabs(A.at(k, k));
            for (std::size_t i = k + 1; i < n; i++) {
                double v = std::fabs(A.at(i, k));
                if (v > maxv) { maxv = v; piv = i; }
            }
            if (maxv < 1e-15) {
                throw std::runtime_error("MatN::inverse: singular matrix");
            }
            if (piv != k) {
                for (std::size_t j = 0; j < n; j++) {
                    std::swap(A.at(k, j), A.at(piv, j));
                    std::swap(I.at(k, j), I.at(piv, j));
                }
            }
            // Scale pivot row to make A(k,k) = 1
            double pivval = A.at(k, k);
            for (std::size_t j = 0; j < n; j++) {
                A.at(k, j) /= pivval;
                I.at(k, j) /= pivval;
            }
            // Eliminate other rows
            for (std::size_t i = 0; i < n; i++) {
                if (i == k) continue;
                double f = A.at(i, k);
                if (f == 0.0) continue;
                for (std::size_t j = 0; j < n; j++) {
                    A.at(i, j) -= f * A.at(k, j);
                    I.at(i, j) -= f * I.at(k, j);
                }
            }
        }
        return I;
    }
};

// Convenience: column-vector-as-matrix builder
inline MatN colvec(std::size_t n) { return MatN(n, 1); }

} // namespace rocket6dof

#endif
