//  json.cpp  --  Minimal JSON parser implementation
//
//  Recursive-descent parser with explicit line/column tracking.  No
//  external dependencies, no exotic C++ features.
//
//  Grammar (more permissive than strict JSON):
//      value   := object | array | string | number | bool | null
//      object  := '{' [pair (',' pair)* [',']] '}'
//      pair    := string ':' value
//      array   := '[' [value (',' value)* [',']] ']'
//      string  := '"' (any char except " and \\, or escape) '"'
//      number  := [-]? digits [. digits] [[eE] [+-]? digits]
//      bool    := 'true' | 'false'
//      null    := 'null'
//
//  Extensions: trailing commas in arrays/objects are tolerated, and
//  single-line `//` comments are stripped before parsing.

#include "json.h"

#include <cctype>
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <sstream>
#include <stdexcept>

namespace rocket6dof {
namespace json {

const Value& Value::missing() {
    static const Value m;     // type_ = MISSING by default
    return m;
}

const Value& Value::operator[](const std::string& key) const {
    if (type_ != OBJECT) return missing();
    auto it = obj_.find(key);
    if (it == obj_.end()) return missing();
    return it->second;
}

const Value& Value::operator[](size_t idx) const {
    if (type_ != ARRAY || idx >= arr_.size()) return missing();
    return arr_[idx];
}

size_t Value::size() const {
    if (type_ == ARRAY)  return arr_.size();
    if (type_ == OBJECT) return obj_.size();
    return 0;
}

Value& Value::set(const std::string& key, Value v) {
    if (type_ != OBJECT) {
        type_ = OBJECT;
        obj_.clear();
    }
    obj_[key] = std::move(v);
    return obj_[key];
}

void Value::push(Value v) {
    if (type_ != ARRAY) {
        type_ = ARRAY;
        arr_.clear();
    }
    arr_.push_back(std::move(v));
}

//  set_path: navigate a dotted-and-bracket path and write a number
//
//  Tokenizes the path into segments.  Each segment is either an object
//  key ("foo") or an array index ("[3]" attached to the previous key).
//  We accept paths like:
//
//    "a"            -> root["a"]
//    "a.b"          -> root["a"]["b"]
//    "a[2]"         -> root["a"][2]
//    "a.b[2].c"     -> root["a"]["b"][2]["c"]
//
//  Intermediate objects are auto-created (so we can write a path that
//  goes through keys not present in the config).  Arrays are NOT
//  auto-created; if a path segment expects an array index, the array
//  must already exist with at least that many elements.
namespace {

// Parse the next segment from `path` starting at `pos`.  On return,
// `key` holds the object-key part (may be empty), `has_index` is true
// if an index follows, and `index` is the index.  `pos` advances past
// the segment.
//
// Examples:
//   "a"        -> key="a", has_index=false
//   "a[3]"     -> key="a", has_index=true, index=3
//   "[3]"      -> key="",  has_index=true, index=3 (used after another [n])
bool parse_segment(const std::string& path, size_t& pos,
                   std::string& key, bool& has_index, size_t& index)
{
    key.clear();
    has_index = false;
    index = 0;
    // Read key chars until '.' or '[' or end
    while (pos < path.size() && path[pos] != '.' && path[pos] != '[') {
        key += path[pos++];
    }
    if (pos < path.size() && path[pos] == '[') {
        pos++;   // consume [
        std::string num;
        while (pos < path.size() && path[pos] != ']') num += path[pos++];
        if (pos >= path.size()) return false;
        pos++;   // consume ]
        if (num.empty()) return false;
        index = static_cast<size_t>(std::atoi(num.c_str()));
        has_index = true;
    }
    if (pos < path.size() && path[pos] == '.') pos++;   // consume .
    return true;
}

} // anon

bool Value::set_path(const std::string& path, double value) {
    // Walk the path, creating intermediate objects as needed.
    Value* node = this;
    size_t pos = 0;
    while (pos < path.size()) {
        std::string key;
        bool has_index;
        size_t index;
        if (!parse_segment(path, pos, key, has_index, index)) return false;

        // First descend into the object at this key (if any)
        if (!key.empty()) {
            if (node->type_ != OBJECT) {
                node->type_ = OBJECT;
                node->obj_.clear();
            }
            // Find or create the entry
            auto it = node->obj_.find(key);
            if (it == node->obj_.end()) {
                // Create as object by default; if next segment doesn't
                // descend, this gets overwritten below.
                node->obj_[key] = Value::makeObject();
                it = node->obj_.find(key);
            }
            node = &it->second;
        }

        // Then descend into an array element (if [index] was present)
        if (has_index) {
            if (node->type_ != ARRAY) return false;
            if (index >= node->arr_.size()) return false;
            node = &node->arr_[index];
        }
    }
    // Write the value
    node->type_ = NUMBER;
    node->num_  = value;
    return true;
}

bool Value::set_path(const std::string& path, const std::string& value) {
    // Identical walk to the numeric overload; only the leaf write
    // differs.  Duplicated rather than templated to keep this file
    // free of template machinery and to allow the numeric path to
    // stay the hot-loop case (called from the GUI's edit_double for
    // every InputDouble interaction).
    Value* node = this;
    size_t pos = 0;
    while (pos < path.size()) {
        std::string key;
        bool has_index;
        size_t index;
        if (!parse_segment(path, pos, key, has_index, index)) return false;

        if (!key.empty()) {
            if (node->type_ != OBJECT) {
                node->type_ = OBJECT;
                node->obj_.clear();
            }
            auto it = node->obj_.find(key);
            if (it == node->obj_.end()) {
                node->obj_[key] = Value::makeObject();
                it = node->obj_.find(key);
            }
            node = &it->second;
        }

        if (has_index) {
            if (node->type_ != ARRAY) return false;
            if (index >= node->arr_.size()) return false;
            node = &node->arr_[index];
        }
    }
    node->type_ = STRING;
    node->str_  = value;
    return true;
}

const Object& Value::obj() const {
    static const Object empty;
    return type_ == OBJECT ? obj_ : empty;
}

const Array& Value::arr() const {
    static const Array empty;
    return type_ == ARRAY ? arr_ : empty;
}

std::string Value::dump(int indent) const {
    std::string pad(indent * 2, ' ');
    std::string pad2((indent + 1) * 2, ' ');
    std::ostringstream o;
    switch (type_) {
        case NUL:    o << "null"; break;
        case BOOL:   o << (bool_ ? "true" : "false"); break;
        case NUMBER: o << num_; break;
        case STRING: o << "\"" << str_ << "\""; break;
        case ARRAY: {
            o << "[";
            bool first = true;
            for (const auto& v : arr_) {
                if (!first) o << ",";
                o << "\n" << pad2 << v.dump(indent + 1);
                first = false;
            }
            if (!first) o << "\n" << pad;
            o << "]";
            break;
        }
        case OBJECT: {
            o << "{";
            bool first = true;
            for (const auto& kv : obj_) {
                if (!first) o << ",";
                o << "\n" << pad2 << "\"" << kv.first << "\": " << kv.second.dump(indent + 1);
                first = false;
            }
            if (!first) o << "\n" << pad;
            o << "}";
            break;
        }
        case MISSING: o << "<missing>"; break;
    }
    return o.str();
}

//  Parser

namespace {

class Parser {
public:
    Parser(const std::string& text) : src_(text), pos_(0), line_(1), col_(1) {}

    Value parse() {
        skipWS();
        Value v = parseValue();
        skipWS();
        if (pos_ < src_.size()) {
            err("unexpected trailing characters");
        }
        return v;
    }

private:
    const std::string& src_;
    size_t pos_;
    int line_, col_;

    [[noreturn]] void err(const std::string& msg) {
        throw ParseError(msg, line_, col_);
    }

    char peek() const { return pos_ < src_.size() ? src_[pos_] : '\0'; }
    char peek_at(size_t off) const {
        return pos_ + off < src_.size() ? src_[pos_ + off] : '\0';
    }
    void advance() {
        if (pos_ >= src_.size()) return;
        if (src_[pos_] == '\n') { line_++; col_ = 1; }
        else col_++;
        pos_++;
    }

    void skipWS() {
        while (pos_ < src_.size()) {
            char c = src_[pos_];
            if (c == ' ' || c == '\t' || c == '\r' || c == '\n') {
                advance();
            } else if (c == '/' && peek_at(1) == '/') {
                // Line comment: skip to end of line
                while (pos_ < src_.size() && src_[pos_] != '\n') advance();
            } else {
                break;
            }
        }
    }

    Value parseValue() {
        skipWS();
        char c = peek();
        if (c == '{') return parseObject();
        if (c == '[') return parseArray();
        if (c == '"') return parseString();
        if (c == '-' || (c >= '0' && c <= '9')) return parseNumber();
        if (c == 't' || c == 'f') return parseBool();
        if (c == 'n') return parseNull();
        err(std::string("unexpected character '") + c + "'");
    }

    Value parseObject() {
        Value v = Value::makeObject();
        advance();    // consume '{'
        skipWS();
        if (peek() == '}') { advance(); return v; }
        for (;;) {
            skipWS();
            if (peek() != '"') err("expected string key in object");
            Value key = parseString();
            skipWS();
            if (peek() != ':') err("expected ':' after key");
            advance();
            skipWS();
            Value val = parseValue();
            v.set(key.asString(), std::move(val));
            skipWS();
            if (peek() == ',') { advance(); skipWS();
                if (peek() == '}') { advance(); return v; }  // trailing comma OK
                continue;
            }
            if (peek() == '}') { advance(); return v; }
            err("expected ',' or '}' in object");
        }
    }

    Value parseArray() {
        Value v = Value::makeArray();
        advance();    // consume '['
        skipWS();
        if (peek() == ']') { advance(); return v; }
        for (;;) {
            skipWS();
            v.push(parseValue());
            skipWS();
            if (peek() == ',') { advance(); skipWS();
                if (peek() == ']') { advance(); return v; }  // trailing comma OK
                continue;
            }
            if (peek() == ']') { advance(); return v; }
            err("expected ',' or ']' in array");
        }
    }

    Value parseString() {
        if (peek() != '"') err("expected string");
        advance();
        std::string out;
        while (pos_ < src_.size()) {
            char c = src_[pos_];
            if (c == '"') { advance(); return Value::makeString(out); }
            if (c == '\\') {
                advance();
                if (pos_ >= src_.size()) err("unterminated escape in string");
                char e = src_[pos_];
                switch (e) {
                    case '"': out += '"'; break;
                    case '\\': out += '\\'; break;
                    case '/':  out += '/'; break;
                    case 'n':  out += '\n'; break;
                    case 't':  out += '\t'; break;
                    case 'r':  out += '\r'; break;
                    case 'b':  out += '\b'; break;
                    case 'f':  out += '\f'; break;
                    default:
                        err(std::string("unknown escape \\") + e);
                }
                advance();
            } else if (c == '\n') {
                err("unterminated string (newline)");
            } else {
                out += c;
                advance();
            }
        }
        err("unterminated string");
    }

    Value parseNumber() {
        size_t start = pos_;
        if (peek() == '-') advance();
        while (peek() >= '0' && peek() <= '9') advance();
        if (peek() == '.') {
            advance();
            while (peek() >= '0' && peek() <= '9') advance();
        }
        if (peek() == 'e' || peek() == 'E') {
            advance();
            if (peek() == '+' || peek() == '-') advance();
            while (peek() >= '0' && peek() <= '9') advance();
        }
        std::string num = src_.substr(start, pos_ - start);
        if (num.empty() || num == "-") err("malformed number");
        return Value::makeNumber(std::atof(num.c_str()));
    }

    Value parseBool() {
        if (src_.compare(pos_, 4, "true") == 0) {
            for (int i = 0; i < 4; i++) advance();
            return Value::makeBool(true);
        }
        if (src_.compare(pos_, 5, "false") == 0) {
            for (int i = 0; i < 5; i++) advance();
            return Value::makeBool(false);
        }
        err("expected 'true' or 'false'");
    }

    Value parseNull() {
        if (src_.compare(pos_, 4, "null") == 0) {
            for (int i = 0; i < 4; i++) advance();
            return Value::makeNull();
        }
        err("expected 'null'");
    }
};

} // anon

Value parse_string(const std::string& text) {
    Parser p(text);
    return p.parse();
}

Value parse_file(const std::string& path) {
    std::ifstream f(path);
    if (!f.is_open()) {
        throw ParseError("cannot open file: " + path, 0, 0);
    }
    std::ostringstream ss;
    ss << f.rdbuf();
    return parse_string(ss.str());
}

} // namespace json
} // namespace rocket6dof
