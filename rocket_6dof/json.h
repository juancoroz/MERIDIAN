//  json.h  --  Minimal JSON parser for mission config files
//
//  A small, header-light, no-dependency JSON parser sufficient to drive
//  the rocket_6dof simulation from an external config file.  Supports
//  the standard JSON value types -- objects, arrays, strings, numbers,
//  booleans, null -- plus single-line `//` comments (JSON5-style; the
//  cost of supporting these is two lines of code and they are very
//  useful in config files).
//
//  The parser builds a tree of `Value` nodes.  Each Value carries its
//  type tag and a discriminated payload.  Lookups via the typed
//  accessors (asObject, asNumber, asString, etc) throw nothing -- they
//  return defaults if the type doesn't match, which makes config-file
//  loading straightforward: "fetch this key as a double, or use 1.0 if
//  it isn't there or isn't a number."
//
//  Parse errors throw a `ParseError` carrying a human-readable message
//  plus the line/column where the error was first detected.
//
//  Usage:
//
//      rocket6dof::json::Value v = rocket6dof::json::parse_file("mission.json");
//      double mass = v["propulsion"]["vmass0"].asNumber(1000.0);
//      std::string aero_file = v["aerodynamics"]["aero_file"].asString("aero.txt");
//
//  This is not a fast or feature-complete JSON library; it's a config
//  loader with line numbers on errors.  ~250 lines of implementation.

#ifndef ROCKET6DOF_JSON_H
#define ROCKET6DOF_JSON_H

#include <map>
#include <memory>
#include <stdexcept>
#include <string>
#include <vector>

namespace rocket6dof {
namespace json {

class Value;
using Object = std::map<std::string, Value>;
using Array  = std::vector<Value>;

enum Type {
    NUL,         // JSON null
    BOOL,        // true / false
    NUMBER,      // any numeric literal (always stored as double)
    STRING,
    ARRAY,
    OBJECT,
    MISSING      // sentinel: key not present (returned from object[k] for missing k)
};

class ParseError : public std::runtime_error {
public:
    int line, col;
    ParseError(const std::string& msg, int ln, int cn)
        : std::runtime_error(msg), line(ln), col(cn) {}
};

class Value {
public:
    Value() : type_(MISSING), num_(0.0), bool_(false) {}
    explicit Value(Type t) : type_(t), num_(0.0), bool_(false) {}

    // Factory helpers
    static Value makeNull()              { Value v(NUL);    return v; }
    static Value makeBool(bool b)        { Value v(BOOL);   v.bool_ = b; return v; }
    static Value makeNumber(double d)    { Value v(NUMBER); v.num_  = d; return v; }
    static Value makeString(std::string s) { Value v(STRING); v.str_ = std::move(s); return v; }
    static Value makeArray()             { Value v(ARRAY);  return v; }
    static Value makeObject()            { Value v(OBJECT); return v; }

    Type type() const { return type_; }
    bool exists() const { return type_ != MISSING; }
    bool isObject() const { return type_ == OBJECT; }
    bool isArray()  const { return type_ == ARRAY; }
    bool isNumber() const { return type_ == NUMBER; }
    bool isString() const { return type_ == STRING; }
    bool isBool()   const { return type_ == BOOL; }
    bool isNull()   const { return type_ == NUL; }

    // Typed accessors with default values.  Defaults are used when:
    //   - the key was missing from the parent object
    //   - the value exists but isn't the requested type
    // This is what makes partial config files easy: ask for what you
    // want, take the default otherwise.
    double      asNumber(double def = 0.0) const {
        return type_ == NUMBER ? num_ : def;
    }
    bool        asBool(bool def = false) const {
        return type_ == BOOL ? bool_ : def;
    }
    std::string asString(const std::string& def = "") const {
        return type_ == STRING ? str_ : def;
    }
    int         asInt(int def = 0) const {
        return type_ == NUMBER ? static_cast<int>(num_) : def;
    }

    // Container accessors.  `operator[](key)` returns a Value of type
    // MISSING if the key isn't found, so chained access is safe.
    const Value& operator[](const std::string& key) const;
    const Value& operator[](size_t idx) const;
    size_t size() const;

    // Mutating: build up Values programmatically (used by parser and
    // by tests).
    Value& set(const std::string& key, Value v);
    void   push(Value v);

    // Walk a dotted-and-bracket path and overwrite the leaf with a
    // numeric value.  Path syntax:
    //
    //   "propulsion.vmass0"          -- nested object key access
    //   "ins.bias_accel[0]"          -- index into an array element
    //   "guidance.theta_com_end"     -- multiple levels
    //
    // Intermediate objects are created on demand if they don't exist.
    // Array indices must be in bounds (path doesn't grow arrays).
    // Returns true on success, false if a non-numeric segment can't be
    // navigated (e.g. trying to index into a scalar).
    bool set_path(const std::string& path, double value);

    // Overload: write a string leaf at the same path syntax.  Useful
    // for fields like `aerodynamics.aero_file` that hold filenames
    // rather than numbers.  Same return-value semantics as the
    // numeric overload.
    bool set_path(const std::string& path, const std::string& value);

    // Get the object as a map (read-only).  Returns empty map if not
    // OBJECT.
    const Object& obj() const;
    const Array&  arr() const;

    // Pretty print (for debugging)
    std::string dump(int indent = 0) const;

private:
    Type type_;
    double num_;
    bool bool_;
    std::string str_;
    Object obj_;
    Array arr_;
    static const Value& missing();   // Singleton MISSING value
};

// ---- Top-level parse functions ----
// parse_string parses a string containing JSON text.  parse_file
// reads the file from disk and parses it.  Both throw ParseError on
// any syntactic problem.
Value parse_string(const std::string& text);
Value parse_file(const std::string& path);

} // namespace json
} // namespace rocket6dof

#endif
