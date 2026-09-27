#include "fvm/post/ChannelStatistics.h"

#include <algorithm>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <stdexcept>

namespace cfd {

ChannelStatistics::ChannelStatistics(const IncompressibleFlow& flow, int streamwise, int normal, scalar yWall0,
                                     scalar yWall1)
    : flow_(flow), sx_(streamwise), ny_(normal), sz_(3 - streamwise - normal) {
    if (sx_ == ny_ || sz_ < 0 || sz_ > 2) throw std::runtime_error("statistics: streamwise and normal axes must differ");
    const Mesh& m = flow.mesh();
    // 层：到最近壁面距离（按网格尺度量化后去重），全局一致
    const scalar H = std::abs(yWall1 - yWall0);
    const scalar tol = 1e-9 * H;
    std::vector<double> d;
    for (label c = 0; c < m.nCells(); ++c) {
        const scalar y = m.C()[c][ny_];
        d.push_back(std::min(std::abs(y - yWall0), std::abs(yWall1 - y)));
    }
    // 收集各进程的层坐标（先本地去重）
    std::vector<double> loc = d;
    std::sort(loc.begin(), loc.end());
    std::vector<double> uniq;
    for (double v : loc)
        if (uniq.empty() || v - uniq.back() > tol) uniq.push_back(v);
    // 全局并集：各进程依次广播
    std::vector<double> all;
    for (int r = 0; r < par::size(); ++r) {
        glabel n = r == par::rank() ? glabel(uniq.size()) : 0;
        par::broadcast(&n, 1, r);
        std::vector<double> buf(r == par::rank() ? uniq : std::vector<double>(std::size_t(n)));
        par::broadcast(buf.data(), int(n), r);
        all.insert(all.end(), buf.begin(), buf.end());
    }
    std::sort(all.begin(), all.end());
    for (double v : all)
        if (layerY_.empty() || v - layerY_.back() > tol) layerY_.push_back(v);
    cellLayer_.resize(m.nCells());
    upper_.resize(m.nCells());
    for (label c = 0; c < m.nCells(); ++c)
        upper_[c] = std::abs(yWall1 - m.C()[c][ny_]) < std::abs(m.C()[c][ny_] - yWall0);
    for (label c = 0; c < m.nCells(); ++c) {
        auto it = std::lower_bound(layerY_.begin(), layerY_.end(), d[c] - tol);
        cellLayer_[c] = label(it - layerY_.begin());
    }
    acc_.assign(layerY_.size() * NQ, 0.0);
}

void ChannelStatistics::sample(scalar dt) {
    const Mesh& m = flow_.mesh();
    const auto& U = flow_.U();
    const auto& p = flow_.p();
    const auto& nut = flow_.nut();
    const auto& G = flow_.gradU();
    const int ax[3] = {sx_, sz_, ny_};  // OpenLB 顺序：流向 u、展向 v、法向 w
    for (label c = 0; c < m.nCells(); ++c) {
        double* a = acc_.data() + std::size_t(cellLayer_[c]) * NQ;
        const scalar w = m.V()[c] * dt;
        const Vec3& u = U[c];
        // 上半通道镜像到下半：法向分量反号（否则两半的剪应力 uw 相互抵消）
        const scalar sg[3] = {1, 1, upper_[c] ? -1.0 : 1.0};
        const scalar ui[3] = {u[ax[0]], u[ax[1]], sg[2] * u[ax[2]]};
        a[0] += w;
        for (int k = 0; k < 3; ++k) a[1 + k] += w * ui[k];
        int q = 4;
        for (int i = 0; i < 3; ++i)
            for (int j = i; j < 3; ++j) a[q++] += w * ui[i] * ui[j];
        a[10] += w * p[c];
        a[11] += w * p[c] * p[c];
        if (!nut.empty() && !G.empty()) {
            // 亚格子应力 τ_ij = −2 ν_t S_ij
            const Tensor S = symm(G[c]);
            q = 12;
            for (int i = 0; i < 3; ++i)
                for (int j = i; j < 3; ++j) a[q++] += w * (-2.0 * nut[c] * sg[i] * sg[j] * S(ax[i], ax[j]));
        }
    }
    // 壁面切应力（面积平均，所有壁面）
    const auto nuEff = flow_.nuEffFaces();
    double tw[2] = {0, 0};
    for (const auto& pt : m.patches()) {
        if (pt.type != PatchType::Wall) continue;
        for (label f = pt.start; f < pt.end(); ++f) {
            const Vec3 n = m.Sf()[f] / m.magSf()[f];
            Vec3 du = U[m.owner()[f]] - U.bValue(f);
            du -= dot(du, n) * n;
            tw[0] += m.magSf()[f] * nuEff[f] * m.nonOrthDeltaCoeffs()[f] * mag(du);
            tw[1] += m.magSf()[f];
        }
    }
    par::allSumInPlace(tw, 2);
    if (tw[1] > 0) tauAcc_ += dt * tw[0] / tw[1];
    time_ += dt;
}

scalar ChannelStatistics::uTau() const { return time_ > 0 ? std::sqrt(tauAcc_ / time_) : 0.0; }

void ChannelStatistics::write(const std::string& file) const {
    std::vector<double> g = acc_;
    par::allSumInPlace(g.data(), int(g.size()));
    if (!par::master() || time_ <= 0) return;
    const auto parent = std::filesystem::path(file).parent_path();
    if (!parent.empty()) std::filesystem::create_directories(parent);
    std::ofstream os(file);
    os << "y+;u_tau;uAv+;uu++;uv++;uw++;vv++;vw++;ww++;pRMS;uu_sgs++;uv_sgs++;uw_sgs++;vv_sgs++;vw_sgs++;ww_sgs++\n";
    os << std::setprecision(10);
    const scalar ut = uTau(), nu = flow_.nu();
    const scalar iu2 = 1.0 / (ut * ut);
    for (std::size_t l = 0; l < layerY_.size(); ++l) {
        const double* a = g.data() + l * NQ;
        if (a[0] <= 0) continue;
        const double inv = 1.0 / a[0];
        const double m[3] = {a[1] * inv, a[2] * inv, a[3] * inv};
        double R[6];
        int q = 0;
        for (int i = 0; i < 3; ++i)
            for (int j = i; j < 3; ++j, ++q) R[q] = a[4 + q] * inv - m[i] * m[j];
        const double pm = a[10] * inv;
        const double prms = std::sqrt(std::max(a[11] * inv - pm * pm, 0.0));
        os << layerY_[l] * ut / nu << ';' << ut << ';' << m[0] / ut;
        for (int k = 0; k < 6; ++k) os << ';' << R[k] * iu2;
        os << ';' << prms;
        for (int k = 0; k < 6; ++k) os << ';' << a[12 + k] * inv * iu2;
        os << '\n';
    }
}

} // namespace cfd
