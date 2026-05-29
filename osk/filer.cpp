
#include "filer.h"
#include <fstream>
#include <sstream>
#include <iostream>
#include <cstdlib>

namespace osk {

std::vector<std::string> Filer::tokenise(const std::string& raw) {
    std::vector<std::string> tokens;
    std::string cur;
    for (char c : raw) {
        if (c == '=' || c == ' ' || c == '\t' || c == '\r' || c == '\n') {
            if (!cur.empty()) { tokens.push_back(cur); cur.clear(); }
        } else {
            cur.push_back(c);
        }
    }
    if (!cur.empty()) tokens.push_back(cur);
    return tokens;
}

Filer::Filer(const std::string& path)
    : path_(path), line0_(0), found_(false)
{
    std::ifstream f(path);
    if (!f) {
        std::cerr << "Filer: cannot open '" << path << "'\n";
        return;
    }
    std::string line;
    while (std::getline(f, line)) {
        lines_.push_back(tokenise(line));
    }
}

void Filer::setLine0(const std::string& tag) {
    if (tag.empty()) {
        line0_ = 0;
        return;
    }
    for (std::size_t i = 0; i < lines_.size(); ++i) {
        if (!lines_[i].empty() && lines_[i][0] == tag) {
            line0_ = i + 1;
            return;
        }
    }
    std::cerr << "Filer: tag '" << tag << "' not found in '"
              << path_ << "'; scanning from top\n";
    line0_ = 0;
}

long Filer::find_key(const std::string& name) {
    for (std::size_t i = line0_; i < lines_.size(); ++i) {
        if (lines_[i].size() >= 2 && lines_[i][0] == name) {
            return static_cast<long>(i);
        }
    }
    return -1;
}

double Filer::getDouble(const std::string& name) {
    long i = find_key(name);
    if (i < 0) {
        std::cerr << "Filer: '" << name << "' not found in '"
                  << path_ << "'\n";
        found_ = false;
        return 0.0;
    }
    found_ = true;
    return std::atof(lines_[i][1].c_str());
}

int Filer::getInt(const std::string& name) {
    long i = find_key(name);
    if (i < 0) {
        std::cerr << "Filer: '" << name << "' not found in '"
                  << path_ << "'\n";
        found_ = false;
        return 0;
    }
    found_ = true;
    return std::atoi(lines_[i][1].c_str());
}

std::string Filer::getString(const std::string& name) {
    long i = find_key(name);
    if (i < 0) {
        std::cerr << "Filer: '" << name << "' not found in '"
                  << path_ << "'\n";
        found_ = false;
        return std::string();
    }
    found_ = true;
    return lines_[i][1];
}

}
