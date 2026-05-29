
#include "table.h"
#include <fstream>
#include <sstream>
#include <iostream>
#include <iomanip>
#include <algorithm>
#include <cstdlib>

namespace osk {

namespace {

class TableParser {
public:
    TableParser(const std::string& path, const std::string& tag)
        : ok_(false)
    {
        std::ifstream f(path);
        if (!f) {
            std::cerr << "Table: cannot open '" << path << "'\n";
            return;
        }
        std::string line;

        bool found = false;
        while (std::getline(f, line)) {
            std::istringstream iss(line);
            std::string first;
            if (iss >> first && first == tag) { found = true; break; }
        }
        if (!found) {
            std::cerr << "Table: tag '" << tag << "' not found in '"
                      << path << "'\n";
            return;
        }

        std::string rest((std::istreambuf_iterator<char>(f)),
                          std::istreambuf_iterator<char>());
        iss_.str(rest);
        ok_ = true;
    }

    bool ok() const { return ok_; }

    bool next_token(std::string& out) {
        return static_cast<bool>(iss_ >> out);
    }
    bool next_int(int& out) {
        return static_cast<bool>(iss_ >> out);
    }
    bool next_double(double& out) {
        return static_cast<bool>(iss_ >> out);
    }
private:
    std::istringstream iss_;
    bool ok_;
};

double interp1(const std::vector<double>& x,
               const std::vector<double>& y, double xq)
{
    if (x.empty()) return 0.0;
    if (xq <= x.front()) return y.front();
    if (xq >= x.back())  return y.back();

    auto it = std::upper_bound(x.begin(), x.end(), xq);
    std::size_t i = std::distance(x.begin(), it) - 1;
    double t = (xq - x[i]) / (x[i+1] - x[i]);
    return y[i] + t * (y[i+1] - y[i]);
}

void locate(const std::vector<double>& x, double xq,
            std::size_t& i, double& t, bool& clamped)
{
    clamped = false;
    if (x.size() < 2) { i = 0; t = 0.0; clamped = true; return; }
    if (xq <= x.front()) { i = 0;          t = 0.0; clamped = true; return; }
    if (xq >= x.back())  { i = x.size()-1; t = 0.0; clamped = true; return; }
    auto it = std::upper_bound(x.begin(), x.end(), xq);
    i = std::distance(x.begin(), it) - 1;
    t = (xq - x[i]) / (x[i+1] - x[i]);
}

}

Table1::Table1(const std::string& path) : path_(path) {}

bool Table1::read(const std::string& tag, bool echo) {
    tag_ = tag;
    TableParser p(path_, tag);
    if (!p.ok()) return false;

    int n;
    if (!p.next_int(n) || n <= 0) {
        std::cerr << "Table1: bad size for '" << tag << "'\n"; return false;
    }

    std::string lbl1, lbl2;
    if (!p.next_token(lbl1) || !p.next_token(lbl2)) {
        std::cerr << "Table1: missing labels for '" << tag << "'\n"; return false;
    }
    x_.resize(n); y_.resize(n);
    for (int i = 0; i < n; ++i) {
        if (!p.next_double(x_[i]) || !p.next_double(y_[i])) {
            std::cerr << "Table1: short data for '" << tag << "'\n";
            return false;
        }
    }
    if (echo) std::cout << *this;
    return true;
}

double Table1::interp(double x) const {
    return interp1(x_, y_, x);
}

std::ostream& operator<<(std::ostream& os, const Table1& t) {
    os << "Table1 [" << t.tag_ << "] from " << t.path_
       << " (" << t.size() << " points)\n";
    for (int i = 0; i < t.size(); ++i) {
        os << "  " << std::setw(12) << t.x_[i]
           << "  " << std::setw(12) << t.y_[i] << "\n";
    }
    return os;
}

Table2::Table2(const std::string& path) : path_(path) {}

bool Table2::read(const std::string& tag, bool echo) {
    tag_ = tag;
    TableParser p(path_, tag);
    if (!p.ok()) return false;

    int n1, n2;
    if (!p.next_int(n1) || !p.next_int(n2) || n1 <= 0 || n2 <= 0) {
        std::cerr << "Table2: bad size for '" << tag << "'\n"; return false;
    }
    std::string lbl;

    if (!p.next_token(lbl)) return false;
    x1_.resize(n1);
    for (int i = 0; i < n1; ++i) {
        if (!p.next_double(x1_[i])) { std::cerr << "Table2: short x1\n"; return false; }
    }

    if (!p.next_token(lbl)) return false;
    x2_.resize(n2);
    for (int j = 0; j < n2; ++j) {
        if (!p.next_double(x2_[j])) { std::cerr << "Table2: short x2\n"; return false; }
    }

    if (!p.next_token(lbl)) return false;
    z_.resize(static_cast<std::size_t>(n1) * n2);
    for (int i = 0; i < n1; ++i) {
        for (int j = 0; j < n2; ++j) {
            if (!p.next_double(z_[i*n2 + j])) {
                std::cerr << "Table2: short data for '" << tag << "'\n";
                return false;
            }
        }
    }
    if (echo) std::cout << *this;
    return true;
}

double Table2::interp(double x1, double x2) const {
    if (x1_.empty() || x2_.empty()) return 0.0;
    std::size_t i, j; double ti, tj; bool ci, cj;
    locate(x1_, x1, i, ti, ci);
    locate(x2_, x2, j, tj, cj);
    const int N2 = n2();

    auto Z = [&](std::size_t ii, std::size_t jj) {
        return z_[ii * N2 + jj];
    };
    std::size_t i1 = ci ? i : i + 1;
    std::size_t j1 = cj ? j : j + 1;

    double z00 = Z(i,  j );
    double z01 = Z(i,  j1);
    double z10 = Z(i1, j );
    double z11 = Z(i1, j1);

    double z0 = z00 + tj * (z01 - z00);
    double z1 = z10 + tj * (z11 - z10);
    return z0 + ti * (z1 - z0);
}

std::ostream& operator<<(std::ostream& os, const Table2& t) {
    os << "Table2 [" << t.tag_ << "] from " << t.path_
       << " (" << t.n1() << "x" << t.n2() << ")\n";
    for (int i = 0; i < t.n1(); ++i) {
        for (int j = 0; j < t.n2(); ++j) {
            os << "  " << std::setw(10) << t.z_[i*t.n2() + j];
        }
        os << "\n";
    }
    return os;
}

Table3::Table3(const std::string& path) : path_(path) {}

bool Table3::read(const std::string& tag, bool echo) {
    tag_ = tag;
    TableParser p(path_, tag);
    if (!p.ok()) return false;

    int n1, n2, n3;
    if (!p.next_int(n1) || !p.next_int(n2) || !p.next_int(n3) ||
        n1 <= 0 || n2 <= 0 || n3 <= 0) {
        std::cerr << "Table3: bad size for '" << tag << "'\n"; return false;
    }
    std::string lbl;
    if (!p.next_token(lbl)) return false;
    x1_.resize(n1);
    for (int i = 0; i < n1; ++i)
        if (!p.next_double(x1_[i])) { std::cerr << "Table3: short x1\n"; return false; }
    if (!p.next_token(lbl)) return false;
    x2_.resize(n2);
    for (int j = 0; j < n2; ++j)
        if (!p.next_double(x2_[j])) { std::cerr << "Table3: short x2\n"; return false; }
    if (!p.next_token(lbl)) return false;
    x3_.resize(n3);
    for (int k = 0; k < n3; ++k)
        if (!p.next_double(x3_[k])) { std::cerr << "Table3: short x3\n"; return false; }
    if (!p.next_token(lbl)) return false;

    z_.resize(static_cast<std::size_t>(n1) * n2 * n3);
    for (int k = 0; k < n3; ++k) {
        for (int i = 0; i < n1; ++i) {
            for (int j = 0; j < n2; ++j) {
                if (!p.next_double(z_[k*n1*n2 + i*n2 + j])) {
                    std::cerr << "Table3: short data for '" << tag << "'\n";
                    return false;
                }
            }
        }
    }
    if (echo) std::cout << *this;
    return true;
}

double Table3::interp(double x1, double x2, double x3) const {
    if (x1_.empty() || x2_.empty() || x3_.empty()) return 0.0;
    std::size_t i, j, k; double ti, tj, tk; bool ci, cj, ck;
    locate(x1_, x1, i, ti, ci);
    locate(x2_, x2, j, tj, cj);
    locate(x3_, x3, k, tk, ck);
    const int N1 = n1();
    const int N2 = n2();

    auto Z = [&](std::size_t ii, std::size_t jj, std::size_t kk) {
        return z_[kk*N1*N2 + ii*N2 + jj];
    };
    std::size_t i1 = ci ? i : i + 1;
    std::size_t j1 = cj ? j : j + 1;
    std::size_t k1 = ck ? k : k + 1;

    auto bilin_at = [&](std::size_t kq) {
        double z00 = Z(i,  j,  kq);
        double z01 = Z(i,  j1, kq);
        double z10 = Z(i1, j,  kq);
        double z11 = Z(i1, j1, kq);
        double z0 = z00 + tj * (z01 - z00);
        double z1 = z10 + tj * (z11 - z10);
        return z0 + ti * (z1 - z0);
    };
    double a = bilin_at(k);
    double b = bilin_at(k1);
    return a + tk * (b - a);
}

std::ostream& operator<<(std::ostream& os, const Table3& t) {
    os << "Table3 [" << t.tag_ << "] from " << t.path_
       << " (" << t.n1() << "x" << t.n2() << "x" << t.n3() << ")\n";
    for (int k = 0; k < t.n3(); ++k) {
        os << " block k=" << k << " (x3=" << t.x3_[k] << ")\n";
        for (int i = 0; i < t.n1(); ++i) {
            for (int j = 0; j < t.n2(); ++j) {
                os << "  " << std::setw(10)
                   << t.z_[k*t.n1()*t.n2() + i*t.n2() + j];
            }
            os << "\n";
        }
    }
    return os;
}

}
