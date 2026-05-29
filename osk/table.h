
#ifndef OSK_TABLE_H
#define OSK_TABLE_H

#include <string>
#include <vector>
#include <iosfwd>

namespace osk {

class Table1 {
public:
    explicit Table1(const std::string& path);

    bool read(const std::string& tag, bool echo = false);

    double interp(double x) const;
    double operator()(double x) const { return interp(x); }

    int size() const { return static_cast<int>(x_.size()); }

    friend std::ostream& operator<<(std::ostream& os, const Table1& t);

private:
    std::string         path_;
    std::string         tag_;
    std::vector<double> x_, y_;
};

class Table2 {
public:
    explicit Table2(const std::string& path);

    bool read(const std::string& tag, bool echo = false);

    double interp(double x1, double x2) const;
    double operator()(double x1, double x2) const { return interp(x1, x2); }

    int n1() const { return static_cast<int>(x1_.size()); }
    int n2() const { return static_cast<int>(x2_.size()); }

    friend std::ostream& operator<<(std::ostream& os, const Table2& t);

private:
    std::string         path_;
    std::string         tag_;
    std::vector<double> x1_, x2_;
    std::vector<double> z_;
};

class Table3 {
public:
    explicit Table3(const std::string& path);

    bool read(const std::string& tag, bool echo = false);

    double interp(double x1, double x2, double x3) const;
    double operator()(double x1, double x2, double x3) const {
        return interp(x1, x2, x3);
    }

    int n1() const { return static_cast<int>(x1_.size()); }
    int n2() const { return static_cast<int>(x2_.size()); }
    int n3() const { return static_cast<int>(x3_.size()); }

    friend std::ostream& operator<<(std::ostream& os, const Table3& t);

private:
    std::string         path_;
    std::string         tag_;
    std::vector<double> x1_, x2_, x3_;
    std::vector<double> z_;
};

}

#endif
