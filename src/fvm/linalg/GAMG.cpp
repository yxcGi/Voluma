#include "fvm/linalg/GAMG.h"

#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <map>
#include <stdexcept>
#include <unordered_map>

namespace cfd {

namespace {

void buildAdjacency(GamgLevel& L) {
    L.off.assign(L.n + 1, 0);
    for (std::size_t f = 0; f < L.l.size(); ++f) {
        ++L.off[L.l[f] + 1];
        ++L.off[L.u[f] + 1];
    }
    for (label i = 0; i < L.n; ++i) L.off[i + 1] += L.off[i];
    const std::size_t nnz = L.off[L.n];
    L.adjCell.resize(nnz);
    L.adjFace.resize(nnz);
    L.adjIsL.resize(nnz);
    std::vector<label> fill(L.off.begin(), L.off.end() - 1);
    for (std::size_t f = 0; f < L.l.size(); ++f) {
        label k = fill[L.l[f]]++;
        L.adjCell[k] = L.u[f];
        L.adjFace[k] = label(f);
        L.adjIsL[k] = 1;
        k = fill[L.u[f]]++;
        L.adjCell[k] = L.l[f];
        L.adjFace[k] = label(f);
        L.adjIsL[k] = 0;
    }
}

// 两两配对聚合，返回下一层（只含拓扑与权重）
GamgLevel agglomerate(GamgLevel& L) {
    L.agg.assign(L.n, -1);
    label nc = 0;
    for (label i = 0; i < L.n; ++i) {
        if (L.agg[i] >= 0) continue;
        label best = -1;
        scalar bw = -1;
        for (label k = L.off[i]; k < L.off[i + 1]; ++k) {
            const label j = L.adjCell[k];
            if (L.agg[j] < 0 && j != i && L.weight[L.adjFace[k]] > bw) {
                bw = L.weight[L.adjFace[k]];
                best = j;
            }
        }
        if (best >= 0) {
            L.agg[i] = L.agg[best] = nc++;
            continue;
        }
        // 没有未配对的邻居：并入连接最强的邻居所在聚合
        for (label k = L.off[i]; k < L.off[i + 1]; ++k) {
            if (L.weight[L.adjFace[k]] > bw) {
                bw = L.weight[L.adjFace[k]];
                best = L.adjCell[k];
            }
        }
        L.agg[i] = best >= 0 ? L.agg[best] : nc++;
    }
    GamgLevel C;
    C.n = nc;
    L.faceMap.assign(L.l.size(), -1);
    L.faceFlip.assign(L.l.size(), 0);
    std::unordered_map<std::uint64_t, label> fmap;
    fmap.reserve(L.l.size());
    for (std::size_t f = 0; f < L.l.size(); ++f) {
        const label a = L.agg[L.l[f]], b = L.agg[L.u[f]];
        if (a == b) continue;
        const label lo = std::min(a, b), hi = std::max(a, b);
        const std::uint64_t key = (std::uint64_t(lo) << 32) | std::uint64_t(hi);
        auto it = fmap.find(key);
        label cf;
        if (it == fmap.end()) {
            cf = label(C.l.size());
            fmap.emplace(key, cf);
            C.l.push_back(lo);
            C.u.push_back(hi);
            C.weight.push_back(0.0);
        } else {
            cf = it->second;
        }
        L.faceMap[f] = cf;
        L.faceFlip[f] = a != lo;
        C.weight[cf] += L.weight[f];
    }
    return C;
}

int gamgMergeLevels() {
    const char* e = std::getenv("CFD_GAMG_MERGE");
    return e ? std::max(1, std::atoi(e)) : 2;
}

} // namespace

const GamgHierarchy& gamgHierarchy(const Mesh& m) { return m.cached<GamgHierarchy>("GAMG", [&] { return buildHierarchy(m); }); }

GamgHierarchy buildHierarchy(const Mesh& m) {
    GamgHierarchy H;
    GamgLevel L0;
    L0.n = m.nCells();
    for (label f = 0; f < m.nInternalFaces(); ++f) {
        const label o = m.owner()[f], n = m.neighbour()[f];
        if (o >= m.nCells() || n >= m.nCells() || o == n) continue;
        L0.l.push_back(std::min(o, n));
        L0.u.push_back(std::max(o, n));
        L0.weight.push_back(m.magSf()[f] * m.nonOrthDeltaCoeffs()[f]);
        H.meshFace.push_back(f);
        H.ownerIsL.push_back(o < n);
    }
    buildAdjacency(L0);
    H.levels.push_back(std::move(L0));
    constexpr label nCoarsest = 64;
    // 每层做 mergeLevels 次两两配对（粗化比约 2^mergeLevels），与 OpenFOAM pairGAMGAgglomeration 相同
    const int mergeLevels = gamgMergeLevels();
    while (H.levels.back().n > nCoarsest && H.levels.size() < 40) {
        GamgLevel& F = H.levels.back();
        GamgLevel C = agglomerate(F);
        for (int k = 1; k < mergeLevels && C.n > nCoarsest; ++k) {
            buildAdjacency(C);
            GamgLevel C2 = agglomerate(C);
            for (auto& a : F.agg) a = C.agg[a];
            for (std::size_t f = 0; f < F.faceMap.size(); ++f) {
                const label cf = F.faceMap[f];
                if (cf < 0) continue;
                F.faceMap[f] = C.faceMap[cf];
                F.faceFlip[f] = F.faceFlip[f] != C.faceFlip[cf];
            }
            C = std::move(C2);
        }
        if (C.n >= F.n) {  // 无法继续粗化
            F.agg.clear();
            break;
        }
        buildAdjacency(C);
        H.levels.push_back(std::move(C));
    }
    return H;
}

GamgPreconditioner::GamgPreconditioner(const LduSystem& s) : H_(gamgHierarchy(s.mesh)) {
    const int nL = int(H_.levels.size());
    A_.resize(nL);
    // 第 0 层系数
    {
        const auto& L = H_.levels[0];
        auto& A = A_[0];
        A.diag.assign(s.diag.begin(), s.diag.begin() + L.n);
        A.upper.resize(L.l.size());
        A.lower.resize(L.l.size());
        for (std::size_t k = 0; k < L.l.size(); ++k) {
            const label f = H_.meshFace[k];
            A.upper[k] = H_.ownerIsL[k] ? s.upper[f] : s.lower[f];
            A.lower[k] = H_.ownerIsL[k] ? s.lower[f] : s.upper[f];
        }
    }
    // Galerkin 粗化
    for (int lev = 0; lev + 1 < nL; ++lev) {
        const auto& L = H_.levels[lev];
        const auto& F = A_[lev];
        auto& C = A_[lev + 1];
        const auto& LC = H_.levels[lev + 1];
        C.diag.assign(LC.n, 0.0);
        C.upper.assign(LC.l.size(), 0.0);
        C.lower.assign(LC.l.size(), 0.0);
        for (label i = 0; i < L.n; ++i) C.diag[L.agg[i]] += F.diag[i];
        for (std::size_t f = 0; f < L.l.size(); ++f) {
            const label cf = L.faceMap[f];
            if (cf < 0) {
                C.diag[L.agg[L.l[f]]] += F.upper[f] + F.lower[f];
            } else if (!L.faceFlip[f]) {
                C.upper[cf] += F.upper[f];
                C.lower[cf] += F.lower[f];
            } else {
                C.upper[cf] += F.lower[f];
                C.lower[cf] += F.upper[f];
            }
        }
    }
    // 最粗层稠密 LU
    const auto& LC = H_.levels.back();
    const auto& AC = A_.back();
    nCoarse_ = LC.n;
    const label n = nCoarse_;
    lu_.assign(std::size_t(n) * n, 0.0);
    for (label i = 0; i < n; ++i) lu_[std::size_t(i) * n + i] = AC.diag[i];
    for (std::size_t f = 0; f < LC.l.size(); ++f) {
        lu_[std::size_t(LC.l[f]) * n + LC.u[f]] += AC.upper[f];
        lu_[std::size_t(LC.u[f]) * n + LC.l[f]] += AC.lower[f];
    }
    for (label k = 0; k < n; ++k) {
        const scalar piv = lu_[std::size_t(k) * n + k];
        if (std::abs(piv) < VSMALL) throw std::runtime_error("GAMG: singular coarsest matrix");
        for (label i = k + 1; i < n; ++i) {
            const scalar fct = lu_[std::size_t(i) * n + k] / piv;
            lu_[std::size_t(i) * n + k] = fct;
            if (fct == 0.0) continue;
            for (label j = k + 1; j < n; ++j) lu_[std::size_t(i) * n + j] -= fct * lu_[std::size_t(k) * n + j];
        }
    }
    rbuf_.resize(nL);
    xbuf_.resize(nL);
    tmp_.resize(nL);
    for (int lev = 0; lev < nL; ++lev) {
        rbuf_[lev].resize(H_.levels[lev].n);
        xbuf_[lev].resize(H_.levels[lev].n);
        tmp_[lev].resize(H_.levels[lev].n);
    }
}

void GamgPreconditioner::smooth(int lev, const std::vector<scalar>& r, std::vector<scalar>& x, bool forward) const {
    const auto& L = H_.levels[lev];
    const auto& A = A_[lev];
    auto row = [&](label i) {
        scalar s = r[i];
        for (label k = L.off[i]; k < L.off[i + 1]; ++k) {
            const label f = L.adjFace[k];
            s -= (L.adjIsL[k] ? A.upper[f] : A.lower[f]) * x[L.adjCell[k]];
        }
        x[i] = s / A.diag[i];
    };
    if (forward)
        for (label i = 0; i < L.n; ++i) row(i);
    else
        for (label i = L.n; i-- > 0;) row(i);
}

void GamgPreconditioner::vcycle(int lev, const std::vector<scalar>& r, std::vector<scalar>& x) const {
    const int nL = int(A_.size());
    const auto& L = H_.levels[lev];
    if (lev == nL - 1) {
        const label n = nCoarse_;
        for (label i = 0; i < n; ++i) x[i] = r[i];
        for (label i = 0; i < n; ++i)
            for (label k = 0; k < i; ++k) x[i] -= lu_[std::size_t(i) * n + k] * x[k];
        for (label i = n; i-- > 0;) {
            for (label k = i + 1; k < n; ++k) x[i] -= lu_[std::size_t(i) * n + k] * x[k];
            x[i] /= lu_[std::size_t(i) * n + i];
        }
        return;
    }
    const auto& A = A_[lev];
    constexpr int nSweeps = 2;
    std::fill(x.begin(), x.begin() + L.n, 0.0);
    for (int s = 0; s < nSweeps; ++s) smooth(lev, r, x, true);
    // 残差并限制
    auto& res = tmp_[lev];
    for (label i = 0; i < L.n; ++i) res[i] = r[i] - A.diag[i] * x[i];
    for (std::size_t f = 0; f < L.l.size(); ++f) {
        res[L.l[f]] -= A.upper[f] * x[L.u[f]];
        res[L.u[f]] -= A.lower[f] * x[L.l[f]];
    }
    auto& rc = rbuf_[lev + 1];
    auto& xc = xbuf_[lev + 1];
    std::fill(rc.begin(), rc.end(), 0.0);
    for (label i = 0; i < L.n; ++i) rc[L.agg[i]] += res[i];
    vcycle(lev + 1, rc, xc);
    for (label i = 0; i < L.n; ++i) x[i] += xc[L.agg[i]];
    for (int s = 0; s < nSweeps; ++s) smooth(lev, r, x, false);
}

void GamgPreconditioner::apply(const std::vector<scalar>& r, std::vector<scalar>& w) const {
    auto& x0 = xbuf_[0];
    vcycle(0, r, x0);
    std::copy(x0.begin(), x0.end(), w.begin());
}

} // namespace cfd
