#include "fvm/io/TecplotWriter.h"

#include <algorithm>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>

namespace cfd {

namespace {

// 本进程的单元 → Tecplot 单元（局部点号，四边形 4 个或六面体 8 个），values 对应复制
struct LocalZone {
    std::vector<Vec3> points;
    std::vector<std::int64_t> elems;  // 每个单元 npe 个
    std::vector<label> cellOfElem;
};

LocalZone buildZone(const Mesh& m, bool& skipped) {
    const VtkPiece& vp = m.vtk();
    LocalZone z;
    z.points = vp.points;
    const bool twoD = m.twoD();
    const int ed = twoD ? m.emptyDir() : 0;
    for (label c = 0; c < m.nCells(); ++c) {
        const std::int64_t b = c == 0 ? 0 : vp.offsets[c - 1], e = vp.offsets[c];
        std::vector<std::int64_t> n(vp.connectivity.begin() + b, vp.connectivity.begin() + e);
        const int t = vp.types[c];
        if (twoD) {
            // 取空方向坐标较小一侧的点，按绕中心角度排序
            if (t == 42) {  // 多面体：点列表可能有重复
                std::sort(n.begin(), n.end());
                n.erase(std::unique(n.begin(), n.end()), n.end());
            }
            scalar zmin = GREAT, zmax = -GREAT;
            for (auto p : n) {
                zmin = std::min(zmin, vp.points[p][ed]);
                zmax = std::max(zmax, vp.points[p][ed]);
            }
            const scalar mid = 0.5 * (zmin + zmax);
            std::vector<std::int64_t> low;
            for (auto p : n)
                if (vp.points[p][ed] < mid) low.push_back(p);
            const int a0 = (ed + 1) % 3, a1 = (ed + 2) % 3;
            Vec3 ctr;
            for (auto p : low) ctr += vp.points[p];
            ctr = ctr / scalar(low.size());
            std::sort(low.begin(), low.end(), [&](std::int64_t p, std::int64_t q) {
                const Vec3 dp = vp.points[p] - ctr, dq = vp.points[q] - ctr;
                return std::atan2(dp[a1], dp[a0]) < std::atan2(dq[a1], dq[a0]);
            });
            if (low.size() == 3) {
                z.elems.insert(z.elems.end(), {low[0], low[1], low[2], low[2]});
                z.cellOfElem.push_back(c);
            } else if (low.size() == 4) {
                z.elems.insert(z.elems.end(), low.begin(), low.end());
                z.cellOfElem.push_back(c);
            } else {
                // 多边形：扇形剖分为退化四边形
                for (std::size_t k = 1; k + 1 < low.size(); ++k) {
                    z.elems.insert(z.elems.end(), {low[0], low[k], low[k + 1], low[k + 1]});
                    z.cellOfElem.push_back(c);
                }
            }
        } else {
            std::vector<std::int64_t> h;
            switch (t) {
            case 12: h = n; break;                                                  // hex
            case 13: h = {n[0], n[1], n[2], n[2], n[3], n[4], n[5], n[5]}; break;  // wedge
            case 10: h = {n[0], n[1], n[2], n[2], n[3], n[3], n[3], n[3]}; break;  // tet
            case 14: h = {n[0], n[1], n[2], n[3], n[4], n[4], n[4], n[4]}; break;  // pyramid
            default: skipped = true; continue;
            }
            z.elems.insert(z.elems.end(), h.begin(), h.end());
            z.cellOfElem.push_back(c);
        }
    }
    return z;
}

} // namespace

void TecplotWriter::write(const std::string& file, const std::string& title) const {
    const Mesh& m = *mesh_;
    bool skipped = false;
    const LocalZone z = buildZone(m, skipped);
    const bool twoD = m.twoD();
    const int npe = twoD ? 4 : 8;
    int ncTotal = 0;
    for (const auto& it : items_) ncTotal += it.nc;
    // 单元值（按 Tecplot 单元复制）
    std::vector<double> vals(z.cellOfElem.size() * std::size_t(ncTotal));
    {
        int off = 0;
        for (const auto& it : items_) {
            const double* d = it.data();
            for (std::size_t e = 0; e < z.cellOfElem.size(); ++e)
                for (int k = 0; k < it.nc; ++k)
                    vals[e * ncTotal + off + k] = d[std::size_t(z.cellOfElem[e]) * it.nc + k];
            off += it.nc;
        }
    }
    skipped = par::allMax(skipped ? 1.0 : 0.0) > 0.5;

    if (!par::master()) {
        par::Packer p;
        p.putVec(z.points);
        p.putVec(z.elems);
        p.putVec(vals);
        par::sendBytes(p.buf, 0, 71);
        return;
    }
    std::vector<std::vector<Vec3>> pts(par::size());
    std::vector<std::vector<std::int64_t>> els(par::size());
    std::vector<std::vector<double>> vs(par::size());
    pts[0] = z.points;
    els[0] = z.elems;
    vs[0] = vals;
    for (int r = 1; r < par::size(); ++r) {
        const auto buf = par::recvBytes(r, 71);
        par::Unpacker u(buf);
        pts[r] = u.getVec<Vec3>();
        els[r] = u.getVec<std::int64_t>();
        vs[r] = u.getVec<double>();
    }
    std::size_t nP = 0, nE = 0;
    for (int r = 0; r < par::size(); ++r) {
        nP += pts[r].size();
        nE += els[r].size() / npe;
    }
    const auto parent = std::filesystem::path(file).parent_path();
    if (!parent.empty()) std::filesystem::create_directories(parent);
    std::ofstream os(file);
    os << std::setprecision(10);
    int dims[3] = {0, 1, 2};
    int nd = 3;
    if (twoD) {
        nd = 2;
        dims[0] = (m.emptyDir() + 1) % 3;
        dims[1] = (m.emptyDir() + 2) % 3;
        if (dims[0] > dims[1]) std::swap(dims[0], dims[1]);
    }
    static const char* axis[3] = {"X", "Y", "Z"};
    os << "Title=\"" << title << "\"\nVARIABLES=";
    for (int k = 0; k < nd; ++k) os << (k ? "," : "") << '"' << axis[dims[k]] << '"';
    static const char* cmp[3] = {"x", "y", "z"};
    for (const auto& it : items_)
        for (int k = 0; k < it.nc; ++k)
            os << ",\"" << it.name << (it.nc == 3 ? std::string("_") + cmp[k] : (it.nc > 1 ? std::to_string(k) : "")) << '"';
    os << "\nZONE T=\"" << title << "\",N=" << nP << ",E=" << nE
       << ",ZONETYPE=" << (twoD ? "FEQUADRILATERAL" : "FEBRICK") << "\nDATAPACKING=BLOCK\n";
    os << "VARLOCATION=([1-" << nd << "]=NODAL";
    if (ncTotal > 0) os << ", [" << nd + 1 << "-" << nd + ncTotal << "]=CELLCENTERED";
    os << ")\n";
    for (int k = 0; k < nd; ++k)
        for (int r = 0; r < par::size(); ++r)
            for (const auto& p : pts[r]) os << p[dims[k]] << '\n';
    for (int k = 0; k < ncTotal; ++k)
        for (int r = 0; r < par::size(); ++r)
            for (std::size_t e = 0; e < vs[r].size() / std::max(ncTotal, 1); ++e) os << vs[r][e * ncTotal + k] << '\n';
    std::size_t base = 0;
    for (int r = 0; r < par::size(); ++r) {
        for (std::size_t e = 0; e < els[r].size(); e += npe) {
            for (int k = 0; k < npe; ++k) os << (k ? " " : "") << base + els[r][e + k] + 1;
            os << '\n';
        }
        base += pts[r].size();
    }
    if (skipped) std::cout << "Warning: Tecplot output skipped polyhedral cells (use VTK for those)\n";
}

} // namespace cfd
