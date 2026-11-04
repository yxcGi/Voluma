#include "fvm/post/CutPlane.h"

#include <algorithm>
#include <cmath>
#include <cstring>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <stdexcept>

namespace cfd {

CutPlane::CutPlane(MeshPtr mesh, const Vec3& point, const Vec3& normal) : mesh_(std::move(mesh)), x0_(point) {
    const scalar mn = mag(normal);
    if (!(mn > 0)) throw std::runtime_error("plane: normal must be non-zero");
    n_ = normal / mn;
    const Mesh& m = *mesh_;
    const VtkPiece& vp = m.vtk();
    for (label c = 0; c < m.nCells(); ++c) {
        const std::int64_t b = c == 0 ? 0 : vp.offsets[c - 1], e = vp.offsets[c];
        scalar smin = GREAT, smax = -GREAT;
        for (std::int64_t k = b; k < e; ++k) {
            const scalar s = dot(vp.points[vp.connectivity[k]] - x0_, n_);
            smin = std::min(smin, s);
            smax = std::max(smax, s);
        }
        const scalar eps = 1e-8 * std::cbrt(m.V()[c]);
        if (smin < eps && smax > eps) {
            cells_.push_back(c);
            const Vec3& C = m.C()[c];
            points_.push_back(C - dot(C - x0_, n_) * n_);
        }
    }
}

glabel CutPlane::nCells() const { return par::allSum(glabel(cells_.size())); }

void CutPlane::writeTable(const std::string& file, const double* d, int nc, const std::string& name) const {
    std::vector<double> v(cells_.size() * nc);
    for (std::size_t k = 0; k < cells_.size(); ++k)
        for (int i = 0; i < nc; ++i) v[k * nc + i] = d[std::size_t(cells_[k]) * nc + i];
    writeTable(file, v, nc, name);
}

void CutPlane::writeTable(const std::string& file, const std::vector<double>& v, int nc, const std::string& name) const {
    const Mesh& m = *mesh_;
    // 记录：全局单元号, x, y, z, 值...
    const int nr = 4 + nc;
    std::vector<double> rec(cells_.size() * nr);
    for (std::size_t k = 0; k < cells_.size(); ++k) {
        double* r = rec.data() + k * nr;
        r[0] = double(m.cellGlobal()[cells_[k]]);
        r[1] = points_[k].x;
        r[2] = points_[k].y;
        r[3] = points_[k].z;
        for (int i = 0; i < nc; ++i) r[4 + i] = v[k * nc + i];
    }
    rec = par::gatherDoublesToMaster(rec);
    if (par::master()) {
        const std::size_t nRows = rec.size() / nr;
        std::vector<std::size_t> order(nRows);
        for (std::size_t k = 0; k < nRows; ++k) order[k] = k;
        std::sort(order.begin(), order.end(), [&](std::size_t a, std::size_t b) { return rec[a * nr] < rec[b * nr]; });
        const auto parent = std::filesystem::path(file).parent_path();
        if (!parent.empty()) std::filesystem::create_directories(parent);
        std::ofstream os(file);
        os << "x,y,z";
        for (int c = 0; c < nc; ++c)
            os << ',' << name << (nc == 1 ? "" : (nc == 3 ? std::string("_") + "xyz"[c] : std::to_string(c)));
        os << '\n' << std::setprecision(10);
        for (std::size_t k : order) {
            const double* r = rec.data() + k * nr;
            os << r[1] << ',' << r[2] << ',' << r[3];
            for (int c = 0; c < nc; ++c) os << ',' << r[4 + c];
            os << '\n';
        }
    }
}

} // namespace cfd
