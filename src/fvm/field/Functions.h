#pragma once
// 边界条件 / 初值用的简单函数：
//   TimeTable<T>     分段线性时间表 v(t)。JSON 写法：常数（数或 [x,y,z]）、[[t0, v0], [t1, v1], ...]、{"file": "v.csv"}（列 t、分量…）
//   ProfileTable<T>  一维剖面表 v(s)，s 为沿某坐标轴或某方向（相对 origin）的坐标。
//                    JSON：{"file": "inlet.csv", "axis": "y"} 或 {"direction": [0,1,0], "origin": [0,0,0]}，
//                    也可内联 {"points": [[s0, v0...], ...]}；"scale" 乘到值上。文件列：s、分量…（空白/逗号分隔，# 注释）
//   coordHash        由坐标（及时间）散列得到的 [-1, 1) 伪随机数，与进程划分无关
// 表都做线性插值，超出范围取端点值。

#include "fvm/core/Json.h"

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <filesystem>
#include <fstream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

namespace cfd {

inline double coordHash(const Vec3& x, int comp) {
    std::uint64_t h = 1469598103934665603ull ^ std::uint64_t(comp);
    for (scalar v : {x.x, x.y, x.z}) {
        const std::int64_t q = std::llround(v * 1e9);
        h = (h ^ std::uint64_t(q)) * 1099511628211ull;
        h ^= h >> 29;
    }
    return double(h >> 11) / double(1ull << 53) * 2 - 1;
}
inline double coordHash(const Vec3& x, int comp, scalar t) {
    std::uint64_t h = 1469598103934665603ull ^ std::uint64_t(comp + 7);
    for (scalar v : {x.x, x.y, x.z, t}) {
        const std::int64_t q = std::llround(v * 1e9);
        h = (h ^ std::uint64_t(q)) * 1099511628211ull;
        h ^= h >> 29;
    }
    return double(h >> 11) / double(1ull << 53) * 2 - 1;
}

namespace detail {
// 读数值表文件：每行若干数（空白、逗号、分号分隔），# 后为注释，含非数字的行（表头）跳过
inline std::vector<std::vector<double>> readColumns(const std::string& file) {
    std::ifstream is(file);
    if (!is) throw std::runtime_error("cannot open table file " + file);
    std::vector<std::vector<double>> rows;
    std::string line;
    while (std::getline(is, line)) {
        const auto hash = line.find('#');
        if (hash != std::string::npos) line.erase(hash);
        for (char& c : line)
            if (c == ',' || c == ';' || c == '\t') c = ' ';
        std::istringstream ss(line);
        std::vector<double> row;
        std::string tok;
        bool bad = false;
        while (ss >> tok) {
            char* e = nullptr;
            const double v = std::strtod(tok.c_str(), &e);
            if (e == tok.c_str() || *e != '\0') {
                bad = true;
                break;
            }
            row.push_back(v);
        }
        if (!bad && !row.empty()) rows.push_back(row);
    }
    return rows;
}
inline std::string resolve(const std::string& dir, const std::string& p) {
    if (p.empty() || std::filesystem::path(p).is_absolute()) return p;
    return (std::filesystem::path(dir) / p).string();
}
template <class T> T fromRow(const std::vector<double>& row, std::size_t first);
template <> inline scalar fromRow<scalar>(const std::vector<double>& row, std::size_t first) {
    if (row.size() < first + 1) throw std::runtime_error("table row needs a value column");
    return row[first];
}
template <> inline Vec3 fromRow<Vec3>(const std::vector<double>& row, std::size_t first) {
    if (row.size() < first + 3) throw std::runtime_error("table row needs three value columns");
    return {row[first], row[first + 1], row[first + 2]};
}
template <class T> T lerp(const T& a, const T& b, scalar w) { return a + (b - a) * w; }
} // namespace detail

// ------------------------------------------------------------------ TimeTable
template <class T> class TimeTable {
public:
    TimeTable() = default;
    explicit TimeTable(const T& c) : constant_(true), c_(c) {}
    TimeTable(std::vector<scalar> t, std::vector<T> v) : t_(std::move(t)), v_(std::move(v)) { check(); }

    static TimeTable parse(const Json& j, const std::string& caseDir) {
        if (j.isNumber()) return TimeTable(j.as<T>());
        if (j.isObject()) {
            if (!j.has("file")) throw std::runtime_error(j.where() + ": time table needs \"file\"");
            const auto rows = detail::readColumns(detail::resolve(caseDir, j["file"].string()));
            std::vector<scalar> t;
            std::vector<T> v;
            for (auto& r : rows) {
                t.push_back(r.at(0));
                v.push_back(detail::fromRow<T>(r, 1));
            }
            return TimeTable(std::move(t), std::move(v));
        }
        if (j.isArray()) {
            if (j.size() > 0 && j[0].isArray()) {
                std::vector<scalar> t;
                std::vector<T> v;
                for (auto& e : j.array()) {
                    if (!e.isArray() || e.size() != 2) throw std::runtime_error(j.where() + ": time table entries are [t, value]");
                    t.push_back(e[0].number());
                    v.push_back(e[1].as<T>());
                }
                return TimeTable(std::move(t), std::move(v));
            }
            return TimeTable(j.as<T>());
        }
        throw std::runtime_error(j.where() + ": expected a constant, [[t, v], ...] or {\"file\": ...}");
    }

    T operator()(scalar t) const {
        if (constant_) return c_;
        if (t <= t_.front()) return v_.front();
        if (t >= t_.back()) return v_.back();
        const auto it = std::upper_bound(t_.begin(), t_.end(), t);
        const std::size_t i = std::size_t(it - t_.begin());
        const scalar w = (t - t_[i - 1]) / (t_[i] - t_[i - 1]);
        return detail::lerp(v_[i - 1], v_[i], w);
    }
    bool constant() const { return constant_; }

private:
    void check() {
        if (t_.empty() || t_.size() != v_.size()) throw std::runtime_error("time table: empty or mismatched columns");
        for (std::size_t i = 1; i < t_.size(); ++i)
            if (!(t_[i] > t_[i - 1])) throw std::runtime_error("time table: times must increase");
    }
    bool constant_ = false;
    T c_{};
    std::vector<scalar> t_;
    std::vector<T> v_;
};

// ------------------------------------------------------------------ ProfileTable
template <class T> class ProfileTable {
public:
    ProfileTable() = default;
    static ProfileTable parse(const Json& j, const std::string& caseDir) {
        ProfileTable p;
        if (!j.isObject()) throw std::runtime_error(j.where() + ": profile must be an object");
        if (j.has("axis")) {
            const std::string a = j["axis"].string();
            if (a == "x") p.dir_ = {1, 0, 0};
            else if (a == "y") p.dir_ = {0, 1, 0};
            else if (a == "z") p.dir_ = {0, 0, 1};
            else throw std::runtime_error(j.where() + ": profile.axis must be x, y or z");
        } else if (j.has("direction")) {
            p.dir_ = j["direction"].vec3();
            const scalar md = mag(p.dir_);
            if (md <= 0) throw std::runtime_error(j.where() + ": profile.direction must be non-zero");
            p.dir_ /= md;
        } else {
            throw std::runtime_error(j.where() + ": profile needs \"axis\" or \"direction\"");
        }
        p.origin_ = j.get("origin", Vec3{0, 0, 0});
        const scalar scale = j.get("scale", 1.0);
        std::vector<std::pair<scalar, T>> pts;
        if (j.has("file")) {
            for (auto& r : detail::readColumns(detail::resolve(caseDir, j["file"].string())))
                pts.emplace_back(r.at(0), detail::fromRow<T>(r, 1) * scale);
        } else if (j.has("points")) {
            for (auto& e : j["points"].array()) {
                std::vector<double> row;
                for (auto& x : e.array()) row.push_back(x.number());
                pts.emplace_back(row.at(0), detail::fromRow<T>(row, 1) * scale);
            }
        } else {
            throw std::runtime_error(j.where() + ": profile needs \"file\" or \"points\"");
        }
        if (pts.size() < 2) throw std::runtime_error(j.where() + ": profile needs at least two points");
        std::sort(pts.begin(), pts.end(), [](auto& a, auto& b) { return a.first < b.first; });
        for (auto& [s, v] : pts) {
            p.s_.push_back(s);
            p.v_.push_back(v);
        }
        return p;
    }
    T operator()(const Vec3& x) const {
        const scalar s = dot(x - origin_, dir_);
        if (s <= s_.front()) return v_.front();
        if (s >= s_.back()) return v_.back();
        const auto it = std::upper_bound(s_.begin(), s_.end(), s);
        const std::size_t i = std::size_t(it - s_.begin());
        const scalar w = (s - s_[i - 1]) / (s_[i] - s_[i - 1]);
        return detail::lerp(v_[i - 1], v_[i], w);
    }

private:
    Vec3 dir_{0, 1, 0}, origin_{0, 0, 0};
    std::vector<scalar> s_;
    std::vector<T> v_;
};

} // namespace cfd
