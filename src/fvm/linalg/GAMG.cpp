#include "fvm/linalg/GAMG.h"

#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <iostream>
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


// 并行：进程间面、各层接口与全局最粗编号
void buildInterfaces(const Mesh& m, GamgHierarchy& H) {
    const label nC = m.nCells(), nT = m.nTotalCells();
    const int nL = int(H.levels.size());
    const par::Halo& halo = m.halo();
    // 第 0 层单元 → 各层聚合编号（经 halo 得到对方单元在对方进程上的各层编号）
    std::vector<std::vector<scalar>> id(nL, std::vector<scalar>(nT, -1.0));
    {
        std::vector<label> cur(nC);
        for (label c = 0; c < nC; ++c) cur[c] = c;
        for (int k = 0; k < nL; ++k) {
            if (k > 0)
                for (label c = 0; c < nC; ++c) cur[c] = H.levels[k - 1].agg[cur[c]];
            for (label c = 0; c < nC; ++c) id[k][c] = cur[c];
            halo.exchange(id[k]);
        }
        H.cellToCoarse.assign(cur.begin(), cur.end());
    }
    // 全局最粗编号
    std::vector<double> counts(par::size(), 0.0);
    counts[par::rank()] = H.levels.back().n;
    par::allSumInPlace(counts.data(), par::size());
    for (int r = 0; r < par::size(); ++r) {
        if (r < par::rank()) H.coarseOffset += label(counts[r]);
        H.nGlobalCoarse += label(counts[r]);
    }
    {
        std::vector<scalar> gid(nT, -1.0);
        for (label c = 0; c < nC; ++c) gid[c] = H.coarseOffset + H.cellToCoarse[c];
        halo.exchange(gid);
        for (label f = 0; f < m.nInternalFaces(); ++f) {
            const label o = m.owner()[f], n = m.neighbour()[f];
            const bool oLocal = o < nC, nLocal = n < nC;
            if (oLocal == nLocal) continue;
            GamgHierarchy::ProcFace pf;
            pf.face = f;
            pf.localCell = oLocal ? o : n;
            pf.haloCell = oLocal ? n : o;
            pf.remoteCoarse = label(gid[pf.haloCell]);
            pf.localIsOwner = oLocal;
            H.procFaces.push_back(pf);
        }
    }
    // 各层接口。第 k 层 ghost 顺序：沿网格 halo 的收发列表按首次出现排列——发送方第 j 个发送单元
    // 与接收方第 j 个接收单元是同一个单元，所以双方得到的聚合序列一致，无需额外通信
    const auto& nbrs = halo.neighbours();
    const auto& sendL = halo.sendLists();
    const auto& recvL = halo.recvLists();
    H.iface.resize(nL);
    for (int k = 0; k < nL; ++k) {
        auto& I = H.iface[k];
        const label n = H.levels[k].n;
        std::vector<label> ghostOf(nT - nC, -1);
        if (k == 0) {
            I.nGhost = nT - nC;
            for (label h = nC; h < nT; ++h) ghostOf[h - nC] = h - nC;
        } else {
            std::vector<std::vector<label>> snd(nbrs.size()), rcv(nbrs.size());
            for (std::size_t q = 0; q < nbrs.size(); ++q) {
                std::unordered_map<label, label> seen;
                for (label h : recvL[q]) {
                    const label rid = label(id[k][h]);
                    auto it = seen.find(rid);
                    if (it == seen.end()) {
                        it = seen.emplace(rid, I.nGhost++).first;
                        rcv[q].push_back(n + it->second);
                    }
                    ghostOf[h - nC] = it->second;
                }
                std::unordered_map<label, char> sent;
                for (label c : sendL[q]) {
                    const label lid = label(id[k][c]);
                    if (sent.emplace(lid, 1).second) snd[q].push_back(lid);
                }
            }
            I.halo.setup(nbrs, std::move(snd), std::move(rcv));
        }
        std::unordered_map<std::uint64_t, label> entry;
        I.fromProcFace.resize(H.procFaces.size());
        for (std::size_t p = 0; p < H.procFaces.size(); ++p) {
            const auto& pf = H.procFaces[p];
            const label c = label(id[k][pf.localCell]), g = ghostOf[pf.haloCell - nC];
            const std::uint64_t key = (std::uint64_t(c) << 32) | std::uint64_t(g);
            auto it = entry.find(key);
            if (it == entry.end()) {
                it = entry.emplace(key, label(I.cell.size())).first;
                I.cell.push_back(c);
                I.ghost.push_back(g);
            }
            I.fromProcFace[p] = it->second;
        }
        I.off.assign(n + 1, 0);
        for (label c : I.cell) ++I.off[c + 1];
        for (label i = 0; i < n; ++i) I.off[i + 1] += I.off[i];
        I.idx.resize(I.cell.size());
        std::vector<label> fill(I.off.begin(), I.off.end() - 1);
        for (std::size_t e = 0; e < I.cell.size(); ++e) I.idx[fill[I.cell[e]]++] = label(e);
    }
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
    const bool parallel = par::parallel();
    // 最粗层大小：串行 64；并行时全局最粗矩阵约 256 行
    const label nCoarsest = parallel ? std::max<label>(2, 256 / par::size()) : 64;
    const double nCoarsestGlobal = parallel ? 256.0 : 64.0;
    auto globalCount = [&](label n) { return parallel ? par::allSum(double(n)) : double(n); };
    // 每层做 mergeLevels 次两两配对（粗化比约 2^mergeLevels），与 OpenFOAM pairGAMGAgglomeration 相同；
    // 是否继续粗化按全局单元数判断，保证各进程层数相同
    const int mergeLevels = gamgMergeLevels();
    while (H.levels.size() < 40) {
        GamgLevel& F = H.levels.back();
        const double nF = globalCount(F.n);
        if (nF <= nCoarsestGlobal) break;
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
        if (globalCount(C.n) >= nF) {  // 无法继续粗化
            F.agg.clear();
            break;
        }
        buildAdjacency(C);
        H.levels.push_back(std::move(C));
    }
    if (parallel) buildInterfaces(m, H);
    if (std::getenv("CFD_GAMG_INFO") && par::master()) {
        std::cerr << "GAMG levels (rank 0):";
        for (std::size_t k = 0; k < H.levels.size(); ++k)
            std::cerr << ' ' << H.levels[k].n << (k < H.iface.size() ? "+" + std::to_string(H.iface[k].nGhost) : "");
        std::cerr << '\n';
    }
    return H;
}

GamgPreconditioner::GamgPreconditioner(const LduSystem& s) : mesh_(s.mesh), H_(gamgHierarchy(s.mesh)) {
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
    global_ = par::parallel();
    // 进程间接口系数（本进程行上的系数：owner 行取 upper，neighbour 行取 lower）；分段常数插值下
    // 粗层接口系数就是映射到同一条目的第 0 层系数之和
    if (global_) {
        ifc_.resize(H_.iface.size());
        for (std::size_t lev = 0; lev < H_.iface.size(); ++lev) {
            const auto& I = H_.iface[lev];
            ifc_[lev].assign(I.cell.size(), 0.0);
            for (std::size_t p = 0; p < H_.procFaces.size(); ++p) {
                const auto& pf = H_.procFaces[p];
                ifc_[lev][I.fromProcFace[p]] += pf.localIsOwner ? s.upper[pf.face] : s.lower[pf.face];
            }
        }
    }
    // 最粗层稠密 LU（并行时为全局矩阵：本进程的行 + 进程间耦合，求和归约到所有进程）
    const auto& LC = H_.levels.back();
    const auto& AC = A_.back();
    nCoarse_ = global_ ? H_.nGlobalCoarse : LC.n;
    const label n = nCoarse_, off = global_ ? H_.coarseOffset : 0;
    lu_.assign(std::size_t(n) * n, 0.0);
    for (label i = 0; i < LC.n; ++i) lu_[std::size_t(off + i) * n + off + i] = AC.diag[i];
    for (std::size_t f = 0; f < LC.l.size(); ++f) {
        lu_[std::size_t(off + LC.l[f]) * n + off + LC.u[f]] += AC.upper[f];
        lu_[std::size_t(off + LC.u[f]) * n + off + LC.l[f]] += AC.lower[f];
    }
    if (global_) {
        for (const auto& pf : H_.procFaces) {
            const scalar a = pf.localIsOwner ? s.upper[pf.face] : s.lower[pf.face];
            lu_[std::size_t(off + H_.cellToCoarse[pf.localCell]) * n + pf.remoteCoarse] += a;
        }
        par::allSumInPlace(lu_.data(), int(lu_.size()));
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
        const label nG = global_ && lev < int(H_.iface.size()) ? H_.iface[lev].nGhost : 0;
        rbuf_[lev].resize(H_.levels[lev].n);
        xbuf_[lev].resize(H_.levels[lev].n + nG);
        tmp_[lev].resize(H_.levels[lev].n + nG);
    }
}

void GamgPreconditioner::exchange(int lev, std::vector<scalar>& x) const {
    if (!global_) return;
    if (lev == 0)
        mesh_.halo().exchange(x.data());
    else
        H_.iface[lev].halo.exchange(x.data());
}

void GamgPreconditioner::smooth(int lev, const std::vector<scalar>& r, std::vector<scalar>& x, bool forward) const {
    const auto& L = H_.levels[lev];
    const auto& A = A_[lev];
    const GamgHierarchy::Interface* I = global_ ? &H_.iface[lev] : nullptr;
    const scalar* c = global_ ? ifc_[lev].data() : nullptr;
    auto row = [&](label i) {
        scalar s = r[i];
        for (label k = L.off[i]; k < L.off[i + 1]; ++k) {
            const label f = L.adjFace[k];
            s -= (L.adjIsL[k] ? A.upper[f] : A.lower[f]) * x[L.adjCell[k]];
        }
        if (I)
            for (label k = I->off[i]; k < I->off[i + 1]; ++k) {
                const label e = I->idx[k];
                s -= c[e] * x[L.n + I->ghost[e]];
            }
        x[i] = s / A.diag[i];
    };
    if (forward)
        for (label i = 0; i < L.n; ++i) row(i);
    else
        for (label i = L.n; i-- > 0;) row(i);
}

void GamgPreconditioner::residual(int lev, const std::vector<scalar>& r, const std::vector<scalar>& x,
                                  std::vector<scalar>& res) const {
    const auto& L = H_.levels[lev];
    const auto& A = A_[lev];
    for (label i = 0; i < L.n; ++i) res[i] = r[i] - A.diag[i] * x[i];
    for (std::size_t f = 0; f < L.l.size(); ++f) {
        res[L.l[f]] -= A.upper[f] * x[L.u[f]];
        res[L.u[f]] -= A.lower[f] * x[L.l[f]];
    }
    if (global_) {
        const auto& I = H_.iface[lev];
        const auto& c = ifc_[lev];
        for (std::size_t e = 0; e < I.cell.size(); ++e) res[I.cell[e]] -= c[e] * x[L.n + I.ghost[e]];
    }
}

void GamgPreconditioner::vcycle(int lev, const std::vector<scalar>& r, std::vector<scalar>& x) const {
    const int nL = int(A_.size());
    const auto& L = H_.levels[lev];
    if (lev == nL - 1) {
        const label n = nCoarse_;
        std::vector<scalar>& g = coarse_;
        g.assign(n, 0.0);
        const label off = global_ ? H_.coarseOffset : 0;
        for (label i = 0; i < L.n; ++i) g[off + i] = r[i];
        if (global_) par::allSumInPlace(g.data(), int(n));
        for (label i = 0; i < n; ++i)
            for (label k = 0; k < i; ++k) g[i] -= lu_[std::size_t(i) * n + k] * g[k];
        for (label i = n; i-- > 0;) {
            for (label k = i + 1; k < n; ++k) g[i] -= lu_[std::size_t(i) * n + k] * g[k];
            g[i] /= lu_[std::size_t(i) * n + i];
        }
        for (label i = 0; i < L.n; ++i) x[i] = g[off + i];
        return;
    }
    constexpr int nSweeps = 2;
    std::fill(x.begin(), x.end(), 0.0);
    for (int s = 0; s < nSweeps; ++s) {
        if (s > 0) exchange(lev, x);
        smooth(lev, r, x, true);
    }
    // 残差并限制
    exchange(lev, x);
    auto& res = tmp_[lev];
    residual(lev, r, x, res);
    auto& rc = rbuf_[lev + 1];
    auto& xc = xbuf_[lev + 1];
    std::fill(rc.begin(), rc.end(), 0.0);
    for (label i = 0; i < L.n; ++i) rc[L.agg[i]] += res[i];
    vcycle(lev + 1, rc, xc);
    // 粗网格校正（分段常数插值）
    for (label i = 0; i < L.n; ++i) x[i] += xc[L.agg[i]];
    for (int s = 0; s < nSweeps; ++s) {
        exchange(lev, x);
        smooth(lev, r, x, false);
    }
}

void GamgPreconditioner::apply(const std::vector<scalar>& r, std::vector<scalar>& w) const {
    auto& x0 = xbuf_[0];
    vcycle(0, r, x0);
    std::copy(x0.begin(), x0.begin() + H_.levels[0].n, w.begin());
}

} // namespace cfd
