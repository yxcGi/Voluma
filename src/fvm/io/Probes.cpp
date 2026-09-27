#include "fvm/io/Probes.h"

#include <filesystem>
#include <iomanip>

namespace cfd {

Probes::Probes(MeshPtr mesh, std::vector<Vec3> points, std::string dir, std::string name)
    : mesh_(std::move(mesh)), points_(std::move(points)), dir_(std::move(dir)), name_(std::move(name)) {
    const Mesh& m = *mesh_;
    const std::size_t np = points_.size();
    cell_.assign(np, -1);
    centres_.assign(np, Vec3{});
    for (std::size_t k = 0; k < np; ++k) {
        scalar best = GREAT;
        label bc = -1;
        for (label c = 0; c < m.nCells(); ++c) {
            const scalar d = magSqr(m.C()[c] - points_[k]);
            // 距离相同时取全局编号小者，保证与分区无关
            if (d < best || (d == best && bc >= 0 && m.cellGlobal()[c] < m.cellGlobal()[bc])) {
                best = d;
                bc = c;
            }
        }
        const scalar gBest = par::allMin(best);
        const double myId = (bc >= 0 && best == gBest) ? double(m.cellGlobal()[bc]) : 1e300;
        const double gId = par::allMin(myId);
        if (bc >= 0 && double(m.cellGlobal()[bc]) == gId) cell_[k] = bc;
        double xyz[3] = {0, 0, 0};
        if (cell_[k] >= 0) {
            xyz[0] = m.C()[bc].x;
            xyz[1] = m.C()[bc].y;
            xyz[2] = m.C()[bc].z;
        }
        par::allSumInPlace(xyz, 3);
        centres_[k] = {xyz[0], xyz[1], xyz[2]};
    }
    if (par::master()) std::filesystem::create_directories(dir_);
}

std::vector<double> Probes::values(const double* d, int nc) const {
    const std::size_t np = points_.size();
    std::vector<double> v(np * nc, 0.0);
    for (std::size_t k = 0; k < np; ++k)
        if (cell_[k] >= 0)
            for (int c = 0; c < nc; ++c) v[k * nc + c] = d[std::size_t(cell_[k]) * nc + c];
    par::allSumInPlace(v.data(), int(v.size()));
    return v;
}

void Probes::writeTable(const std::string& file, const double* data, int nc, const std::string& name) const {
    const auto v = values(data, nc);
    if (!par::master()) return;
    std::ofstream os(file);
    os << "x,y,z";
    for (int c = 0; c < nc; ++c) os << ',' << name << (nc == 1 ? "" : std::string("_") + "xyz"[c]);
    os << '\n' << std::setprecision(10);
    for (std::size_t k = 0; k < points_.size(); ++k) {
        os << centres_[k].x << ',' << centres_[k].y << ',' << centres_[k].z;
        for (int c = 0; c < nc; ++c) os << ',' << v[k * nc + c];
        os << '\n';
    }
}

void Probes::sample(scalar time) {
    const std::size_t np = points_.size();
    for (auto& it : items_) {
        std::vector<double> v(np * it.nc, 0.0);
        const double* d = it.data();
        for (std::size_t k = 0; k < np; ++k)
            if (cell_[k] >= 0)
                for (int c = 0; c < it.nc; ++c) v[k * it.nc + c] = d[std::size_t(cell_[k]) * it.nc + c];
        par::allSumInPlace(v.data(), int(v.size()));  // 只有拥有者非零，求和精确
        if (!par::master()) continue;
        if (!it.os) {
            it.os = std::make_shared<std::ofstream>(dir_ + "/" + name_ + "_" + it.name + ".csv");
            auto& os = *it.os;
            os << "# probe locations (cell centres):\n";
            for (std::size_t k = 0; k < np; ++k)
                os << "# " << k << ' ' << centres_[k].x << ' ' << centres_[k].y << ' ' << centres_[k].z << '\n';
            os << "time";
            for (std::size_t k = 0; k < np; ++k)
                for (int c = 0; c < it.nc; ++c)
                    os << ',' << k << (it.nc == 1 ? "" : std::string("_") + "xyz"[c]);
            os << '\n';
            os << std::setprecision(10);
        }
        auto& os = *it.os;
        os << time;
        for (double x : v) os << ',' << x;
        os << '\n';
        os.flush();
    }
}

} // namespace cfd
