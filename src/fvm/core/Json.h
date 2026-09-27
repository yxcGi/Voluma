#pragma once
// 轻量 JSON（算例文件格式）：支持对象、数组、数字、字符串、true/false/null，
// 另外允许 // 与 /* */ 注释和尾随逗号，便于手写算例。
//
//   Json j = Json::parseFile("case.json");
//   double nu = j["physics"].get("nu", 1e-3);
//   for (auto& [name, bc] : j["boundary"].items()) ...
//
// 读取过的键会被记录，reportUnused() 可提示拼写错误的键。

#include "fvm/core/Types.h"

#include <cstdlib>
#include <fstream>
#include <map>
#include <memory>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

namespace cfd {

class Json {
public:
    enum class Kind { Null, Bool, Number, String, Array, Object };

    Json() = default;
    static Json parse(const std::string& text, const std::string& source = "<string>") {
        Parser p{text, 0, source};
        p.ws();
        Json j = p.value();
        p.ws();
        if (p.i != text.size()) p.fail("trailing characters");
        return j;
    }
    static Json parseFile(const std::string& file) {
        std::ifstream is(file);
        if (!is) throw std::runtime_error("cannot open " + file);
        std::stringstream ss;
        ss << is.rdbuf();
        return parse(ss.str(), file);
    }

    Kind kind() const { return kind_; }
    bool isNull() const { return kind_ == Kind::Null; }
    bool isObject() const { return kind_ == Kind::Object; }
    bool isArray() const { return kind_ == Kind::Array; }
    bool isString() const { return kind_ == Kind::String; }
    bool isNumber() const { return kind_ == Kind::Number; }
    bool isBool() const { return kind_ == Kind::Bool; }

    bool has(const std::string& k) const {
        if (kind_ != Kind::Object) return false;
        for (auto& [key, v] : obj_)
            if (key == k) return true;
        return false;
    }
    // 对象成员（不存在时返回 null）
    const Json& operator[](const std::string& k) const {
        static const Json null;
        if (kind_ != Kind::Object) return null;
        for (std::size_t n = 0; n < obj_.size(); ++n)
            if (obj_[n].first == k) {
                used_->at(n) = true;
                return obj_[n].second;
            }
        return null;
    }
    const Json& operator[](std::size_t i) const { return arr_.at(i); }
    std::size_t size() const { return kind_ == Kind::Array ? arr_.size() : obj_.size(); }
    const std::vector<Json>& array() const { return arr_; }
    // 对象成员，按文件中的顺序（同时标记为已读）
    const std::vector<std::pair<std::string, Json>>& items() const {
        if (used_) std::fill(used_->begin(), used_->end(), true);
        return obj_;
    }

    double number() const {
        if (kind_ != Kind::Number) throw std::runtime_error(where() + ": expected a number");
        return num_;
    }
    const std::string& string() const {
        if (kind_ != Kind::String) throw std::runtime_error(where() + ": expected a string");
        return str_;
    }
    bool boolean() const {
        if (kind_ != Kind::Bool) throw std::runtime_error(where() + ": expected true/false");
        return b_;
    }
    Vec3 vec3() const {
        if (kind_ == Kind::Number) return {num_, num_, num_};
        if (kind_ != Kind::Array || arr_.size() != 3) throw std::runtime_error(where() + ": expected [x, y, z]");
        return {arr_[0].number(), arr_[1].number(), arr_[2].number()};
    }
    template <class T> T as() const;

    // 带默认值读取
    double get(const std::string& k, double def) const { return has(k) ? (*this)[k].number() : def; }
    int get(const std::string& k, int def) const { return has(k) ? int((*this)[k].number()) : def; }
    bool get(const std::string& k, bool def) const { return has(k) ? (*this)[k].boolean() : def; }
    std::string get(const std::string& k, const char* def) const { return has(k) ? (*this)[k].string() : def; }
    std::string get(const std::string& k, const std::string& def) const { return has(k) ? (*this)[k].string() : def; }
    Vec3 get(const std::string& k, const Vec3& def) const { return has(k) ? (*this)[k].vec3() : def; }

    // 列出未读取的键（递归），返回条数
    int reportUnused(std::ostream& os, const std::string& prefix = "") const {
        int n = 0;
        if (kind_ == Kind::Object) {
            for (std::size_t k = 0; k < obj_.size(); ++k) {
                const std::string path = prefix.empty() ? obj_[k].first : prefix + "." + obj_[k].first;
                if (!used_->at(k)) {
                    os << "Warning: unused case entry '" << path << "'\n";
                    ++n;
                } else {
                    n += obj_[k].second.reportUnused(os, path);
                }
            }
        } else if (kind_ == Kind::Array) {
            for (std::size_t k = 0; k < arr_.size(); ++k)
                n += arr_[k].reportUnused(os, prefix + "[" + std::to_string(k) + "]");
        }
        return n;
    }
    std::string where() const { return src_.empty() ? "json" : src_; }

private:
    struct Parser {
        const std::string& s;
        std::size_t i;
        std::string source;
        [[noreturn]] void fail(const std::string& msg) const {
            int line = 1;
            for (std::size_t k = 0; k < i && k < s.size(); ++k) line += s[k] == '\n';
            throw std::runtime_error(source + ":" + std::to_string(line) + ": " + msg);
        }
        void ws() {
            for (;;) {
                while (i < s.size() && std::isspace(static_cast<unsigned char>(s[i]))) ++i;
                if (i + 1 < s.size() && s[i] == '/' && s[i + 1] == '/') {
                    while (i < s.size() && s[i] != '\n') ++i;
                } else if (i + 1 < s.size() && s[i] == '/' && s[i + 1] == '*') {
                    const auto e = s.find("*/", i + 2);
                    if (e == std::string::npos) fail("unterminated comment");
                    i = e + 2;
                } else {
                    return;
                }
            }
        }
        std::string where() const {
            int line = 1;
            for (std::size_t k = 0; k < i && k < s.size(); ++k) line += s[k] == '\n';
            return source + ":" + std::to_string(line);
        }
        Json value() {
            if (i >= s.size()) fail("unexpected end of input");
            Json j;
            j.src_ = where();
            const char c = s[i];
            if (c == '{') {
                j.kind_ = Kind::Object;
                ++i;
                ws();
                while (i < s.size() && s[i] != '}') {
                    if (s[i] != '"') fail("expected a key string");
                    std::string k = str();
                    ws();
                    if (i >= s.size() || s[i] != ':') fail("expected ':'");
                    ++i;
                    ws();
                    j.obj_.emplace_back(std::move(k), value());
                    ws();
                    if (i < s.size() && s[i] == ',') {
                        ++i;
                        ws();
                    } else if (i < s.size() && s[i] != '}') {
                        fail("expected ',' or '}'");
                    }
                }
                if (i >= s.size()) fail("unterminated object");
                ++i;
                j.used_ = std::make_shared<std::vector<bool>>(j.obj_.size(), false);
            } else if (c == '[') {
                j.kind_ = Kind::Array;
                ++i;
                ws();
                while (i < s.size() && s[i] != ']') {
                    j.arr_.push_back(value());
                    ws();
                    if (i < s.size() && s[i] == ',') {
                        ++i;
                        ws();
                    } else if (i < s.size() && s[i] != ']') {
                        fail("expected ',' or ']'");
                    }
                }
                if (i >= s.size()) fail("unterminated array");
                ++i;
            } else if (c == '"') {
                j.kind_ = Kind::String;
                j.str_ = str();
            } else if (s.compare(i, 4, "true") == 0) {
                j.kind_ = Kind::Bool;
                j.b_ = true;
                i += 4;
            } else if (s.compare(i, 5, "false") == 0) {
                j.kind_ = Kind::Bool;
                i += 5;
            } else if (s.compare(i, 4, "null") == 0) {
                i += 4;
            } else {
                const char* b = s.c_str() + i;
                char* e = nullptr;
                j.num_ = std::strtod(b, &e);
                if (e == b) fail(std::string("unexpected character '") + c + "'");
                j.kind_ = Kind::Number;
                i += std::size_t(e - b);
            }
            return j;
        }
        std::string str() {
            ++i;
            std::string r;
            while (i < s.size() && s[i] != '"') {
                if (s[i] == '\\' && i + 1 < s.size()) {
                    ++i;
                    switch (s[i]) {
                    case 'n': r += '\n'; break;
                    case 't': r += '\t'; break;
                    default: r += s[i];
                    }
                } else {
                    r += s[i];
                }
                ++i;
            }
            if (i >= s.size()) fail("unterminated string");
            ++i;
            return r;
        }
    };

    Kind kind_ = Kind::Null;
    double num_ = 0;
    bool b_ = false;
    std::string str_, src_;
    std::vector<Json> arr_;
    std::vector<std::pair<std::string, Json>> obj_;
    std::shared_ptr<std::vector<bool>> used_;
};

template <> inline scalar Json::as<scalar>() const { return number(); }
template <> inline Vec3 Json::as<Vec3>() const { return vec3(); }

} // namespace cfd
