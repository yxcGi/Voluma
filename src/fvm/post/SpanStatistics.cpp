#include "fvm/post/SpanStatistics.h"

#include <algorithm>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <map>
#include <sstream>
#include <stdexcept>

namespace cfd {

namespace {

// 各进程的键（已排序去重）求全局并集：依次广播
std::vector<double> globalUnion(const std::vector<double>& local, int width) {
    std::vector<double> all;
    for (int r = 0; r < par::size(); ++r) {
        glabel n = r == par::rank() ? glabel(local.size()) : 0;
        par::broadcast(&n, 1, r);
        std::vector<double> buf(r == par::rank() ? local : std::vector<double>(std::size_t(n)));
        if (n > 0) par::broadcast(buf.data(), int(n), r);
        all.insert(all.end(), buf.begin(), buf.end());
    }
    // 按 width 个一组排序去重
    std::vector<std::vector<double>> rows;
    for (std::size_t i = 0; i + width <= all.size(); i += width) rows.emplace_back(all.begin() + i, all.begin() + i + width);
    std::sort(rows.begin(), rows.end());
    rows.erase(std::unique(rows.begin(), rows.end()), rows.end());
    std::vector<double> out;
    for (auto& r : rows) out.insert(out.end(), r.begin(), r.end());
    return out;
}

} // namespace

SpanStatistics::SpanStatistics(const IncompressibleFlow& flow, int streamwise, int normal, const std::string& wallPatch,
                               scalar Uref, scalar H)
    : flow_(flow), sx_(streamwise), nz_(normal), sy_(3 - streamwise - normal), wallName_(wallPatch), Uref_(Uref), H_(H) {
    if (sx_ == nz_ || sy_ < 0 || sy_ > 2) throw std::runtime_error("spanStatistics: streamwise and normal axes must differ");
    const Mesh& m = flow.mesh();
    wallPatch_ = m.findPatch(wallPatch);
    if (wallPatch_ < 0) throw std::runtime_error("spanStatistics: no patch " + wallPatch);
    // 坐标量化（相对域尺度 1e-8）
    scalar L = 0;
    for (label c = 0; c < m.nCells(); ++c) L = std::max({L, std::abs(m.C()[c][sx_]), std::abs(m.C()[c][nz_])});
    L = par::allMax(L);
    const double q = 1e-8 * std::max(L, SMALL);
    auto key = [q](scalar v) { return std::round(v / q); };
    // 单元列
    {
        std::vector<double> loc;
        for (label c = 0; c < m.nCells(); ++c) {
            loc.push_back(key(m.C()[c][sx_]));
            loc.push_back(key(m.C()[c][nz_]));
        }
        std::vector<std::pair<double, double>> pr;
        for (std::size_t i = 0; i < loc.size(); i += 2) pr.emplace_back(loc[i], loc[i + 1]);
        std::sort(pr.begin(), pr.end());
        pr.erase(std::unique(pr.begin(), pr.end()), pr.end());
        std::vector<double> flat;
        for (auto& p : pr) flat.insert(flat.end(), {p.first, p.second});
        const auto g = globalUnion(flat, 2);
        std::map<std::pair<double, double>, label> idx;
        for (std::size_t i = 0; i < g.size(); i += 2) {
            idx[{g[i], g[i + 1]}] = label(binX_.size());
            binX_.push_back(g[i] * q);
            binZ_.push_back(g[i + 1] * q);
        }
        cellBin_.resize(m.nCells());
        for (label c = 0; c < m.nCells(); ++c) cellBin_[c] = idx.at({key(m.C()[c][sx_]), key(m.C()[c][nz_])});
        acc_.assign(binX_.size() * NQ, 0.0);
    }
    // 壁面按流向坐标分组
    {
        const Patch& pt = m.patches()[wallPatch_];
        std::vector<double> loc;
        for (label f = pt.start; f < pt.end(); ++f) loc.push_back(key(m.Cf()[f][sx_]));
        std::sort(loc.begin(), loc.end());
        loc.erase(std::unique(loc.begin(), loc.end()), loc.end());
        const auto g = globalUnion(loc, 1);
        std::map<double, label> idx;
        for (double v : g) {
            idx[v] = label(wX_.size());
            wX_.push_back(v * q);
        }
        std::vector<double> zsum(wX_.size() * 2, 0.0);
        for (label f = pt.start; f < pt.end(); ++f) {
            const label b = idx.at(key(m.Cf()[f][sx_]));
            faceBin_.push_back(b);
            const Vec3 n = m.Sf()[f] / m.magSf()[f];
            Vec3 e{};
            e[sx_] = 1;
            Vec3 t = e - dot(e, n) * n;
            t[sy_] = 0;
            faceT_.push_back(t / std::max(mag(t), VSMALL));
            zsum[2 * b] += m.magSf()[f] * m.Cf()[f][nz_];
            zsum[2 * b + 1] += m.magSf()[f];
        }
        par::allSumInPlace(zsum.data(), int(zsum.size()));
        wZ_.resize(wX_.size());
        for (std::size_t b = 0; b < wX_.size(); ++b) wZ_[b] = zsum[2 * b] / std::max(zsum[2 * b + 1], VSMALL);
        wacc_.assign(wX_.size() * NW, 0.0);
    }
}

void SpanStatistics::sample(scalar dt) {
    const Mesh& m = flow_.mesh();
    const auto& U = flow_.U();
    const auto& p = flow_.p();
    const int ax[3] = {sx_, nz_, sy_};  // 流向 u、法向 v、展向 w（周期山文献惯例）
    for (label c = 0; c < m.nCells(); ++c) {
        double* a = acc_.data() + std::size_t(cellBin_[c]) * NQ;
        const scalar w = m.V()[c] * dt;
        const scalar ui[3] = {U[c][ax[0]], U[c][ax[1]], U[c][ax[2]]};
        a[0] += w;
        for (int k = 0; k < 3; ++k) a[1 + k] += w * ui[k];
        int q = 4;
        for (int i = 0; i < 3; ++i)
            for (int j = i; j < 3; ++j) a[q++] += w * ui[i] * ui[j];
        a[10] += w * p[c];
    }
    const auto nuEff = flow_.nuEffFaces();
    const Patch& pt = m.patches()[wallPatch_];
    for (label f = pt.start; f < pt.end(); ++f) {
        const label i = f - pt.start;
        double* a = wacc_.data() + std::size_t(faceBin_[i]) * NW;
        const scalar w = m.magSf()[f] * dt;
        const Vec3 n = m.Sf()[f] / m.magSf()[f];
        Vec3 du = U[m.owner()[f]] - U.bValue(f);
        du -= dot(du, n) * n;
        const Vec3 tau = nuEff[f] * m.nonOrthDeltaCoeffs()[f] * du;
        a[0] += w;
        a[1] += w * dot(tau, faceT_[i]);
        a[2] += w * p.bValue(f);
        a[3] += w * dot(U[m.owner()[f]], faceT_[i]);
    }
    time_ += dt;
}

std::pair<scalar, scalar> SpanStatistics::write(const std::string& dir, const std::vector<scalar>& stations) const {
    std::vector<double> g = acc_, wg = wacc_;
    par::allSumInPlace(g.data(), int(g.size()));
    par::allSumInPlace(wg.data(), int(wg.size()));
    std::pair<scalar, scalar> sepReat{-1, -1};
    if (time_ <= 0) return sepReat;
    // 壁面分布与分离 / 再附（所有进程都算，返回值一致）
    const scalar q2 = 2.0 / (Uref_ * Uref_);
    std::vector<scalar> Cf(wX_.size()), P(wX_.size()), Ut(wX_.size());
    for (std::size_t b = 0; b < wX_.size(); ++b) {
        const double* a = wg.data() + b * NW;
        const double inv = a[0] > 0 ? 1.0 / a[0] : 0.0;
        Cf[b] = a[1] * inv * q2;
        P[b] = a[2] * inv;
        Ut[b] = a[3] * inv;
    }
    for (std::size_t b = 1; b < wX_.size(); ++b) {
        if (wX_[b] / H_ < 0.1) continue;
        // 线性插值求零点
        const scalar x0 = wX_[b - 1] / H_, x1 = wX_[b] / H_;
        const scalar xz = Cf[b - 1] != Cf[b] ? x0 + (x1 - x0) * Cf[b - 1] / (Cf[b - 1] - Cf[b]) : x1;
        if (sepReat.first < 0 && Cf[b] < 0 && Cf[b - 1] >= 0) sepReat.first = xz;
        else if (sepReat.first >= 0 && sepReat.second < 0 && Cf[b] > 0 && Cf[b - 1] <= 0) sepReat.second = xz;
    }
    if (!par::master()) return sepReat;
    std::filesystem::create_directories(dir);
    {
        std::ofstream os(dir + "/wall_" + wallName_ + ".csv");
        os << "x_h;Cf;Cp;ut_Ub\n" << std::setprecision(10);
        const scalar pRef = P.empty() ? 0.0 : P[0];
        for (std::size_t b = 0; b < wX_.size(); ++b)
            os << wX_[b] / H_ << ';' << Cf[b] << ';' << (P[b] - pRef) * q2 << ';' << Ut[b] / Uref_ << '\n';
    }
    // 各站位：取流向坐标最接近的单元列
    std::vector<double> xs(binX_);
    std::sort(xs.begin(), xs.end());
    xs.erase(std::unique(xs.begin(), xs.end()), xs.end());
    for (scalar st : stations) {
        const scalar xt = st * H_;
        double best = xs.empty() ? 0.0 : xs[0];
        for (double x : xs)
            if (std::abs(x - xt) < std::abs(best - xt)) best = x;
        // 下壁高度：取最接近的壁面组
        std::size_t wb = 0;
        for (std::size_t b = 0; b < wX_.size(); ++b)
            if (std::abs(wX_[b] - best) < std::abs(wX_[wb] - best)) wb = b;
        std::vector<std::pair<double, std::size_t>> col;
        for (std::size_t b = 0; b < binX_.size(); ++b)
            if (binX_[b] == best) col.emplace_back(binZ_[b], b);
        std::sort(col.begin(), col.end());
        std::ostringstream name;
        name << dir << "/profile_xh_" << st << ".csv";
        std::ofstream os(name.str());
        os << "y_h;u_Ub;v_Ub;uu_Ub2;vv_Ub2;ww_Ub2;uv_Ub2;k_Ub2\n" << std::setprecision(10);
        const scalar iu = 1.0 / Uref_, iu2 = iu * iu;
        for (auto& [z, b] : col) {
            const double* a = g.data() + b * NQ;
            if (a[0] <= 0) continue;
            const double inv = 1.0 / a[0];
            const double mu[3] = {a[1] * inv, a[2] * inv, a[3] * inv};
            // a[4..9]：00 01 02 11 12 22（0 流向、1 法向、2 展向）
            const double uu = a[4] * inv - mu[0] * mu[0], uv = a[5] * inv - mu[0] * mu[1];
            const double vv = a[7] * inv - mu[1] * mu[1], ww = a[9] * inv - mu[2] * mu[2];
            os << (z - wZ_[wb]) / H_ << ';' << mu[0] * iu << ';' << mu[1] * iu << ';' << uu * iu2 << ';' << vv * iu2 << ';'
               << ww * iu2 << ';' << uv * iu2 << ';' << 0.5 * (uu + vv + ww) * iu2 << '\n';
        }
    }
    return sepReat;
}

} // namespace cfd
