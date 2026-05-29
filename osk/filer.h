
#ifndef OSK_FILER_H
#define OSK_FILER_H

#include <string>
#include <vector>

namespace osk {

class Filer {
public:

    explicit Filer(const std::string& path);

    void setLine0(const std::string& tag);

    double      getDouble(const std::string& name);
    int         getInt   (const std::string& name);
    std::string getString(const std::string& name);

    bool found() const { return found_; }

    const std::string& path() const { return path_; }

private:
    std::string                                path_;
    std::vector<std::vector<std::string>>      lines_;
    std::size_t                                line0_;
    bool                                       found_;

    long find_key(const std::string& name);

    static std::vector<std::string> tokenise(const std::string& raw);
};

}

#endif
