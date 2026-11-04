#include "fvm/post/FieldAverage.h"

#include "fvm/io/Restart.h"

#include <filesystem>
#include <fstream>

namespace cfd {

void FieldAverage::add(const std::string& name, const std::vector<scalar>& src, bool prime2) {
    Item it;
    it.name = name;
    it.nc = 1;
    it.prime2 = prime2;
    it.src = [&src] { return src.data(); };
    const label nC = mesh_->nCells(), nT = mesh_->nTotalCells();
    it.sum.assign(nC, 0.0);
    it.meanS.assign(nT, 0.0);
    if (prime2) {
        it.sum2.assign(nC, 0.0);
        it.p2S.assign(nT, 0.0);
    }
    items_.push_back(std::move(it));
}

void FieldAverage::add(const std::string& name, const std::vector<Vec3>& src, bool prime2) {
    Item it;
    it.name = name;
    it.nc = 3;
    it.prime2 = prime2;
    it.src = [&src] { return reinterpret_cast<const double*>(src.data()); };
    const label nC = mesh_->nCells(), nT = mesh_->nTotalCells();
    it.sum.assign(std::size_t(nC) * 3, 0.0);
    it.meanV.assign(nT, Vec3{});
    if (prime2) {
        it.sum2.assign(std::size_t(nC) * 9, 0.0);
        it.p2V.assign(nT, Tensor{});
    }
    items_.push_back(std::move(it));
}

void FieldAverage::reset() {
    for (auto& it : items_) {
        std::fill(it.sum.begin(), it.sum.end(), 0.0);
        std::fill(it.sum2.begin(), it.sum2.end(), 0.0);
    }
    time_ = 0;
    n_ = 0;
    update();
}

void FieldAverage::sample(scalar w) {
    if (!(w > 0)) return;
    const label nC = mesh_->nCells();
    for (auto& it : items_) {
        const double* d = it.src();
        if (it.nc == 1) {
            for (label c = 0; c < nC; ++c) {
                it.sum[c] += w * d[c];
                if (it.prime2) it.sum2[c] += w * d[c] * d[c];
            }
        } else {
            for (label c = 0; c < nC; ++c) {
                const double* v = d + 3 * c;
                for (int i = 0; i < 3; ++i) it.sum[3 * c + i] += w * v[i];
                if (it.prime2)
                    for (int i = 0; i < 3; ++i)
                        for (int j = 0; j < 3; ++j) it.sum2[9 * c + 3 * i + j] += w * v[i] * v[j];
            }
        }
    }
    time_ += w;
    ++n_;
    update();
}

void FieldAverage::update() {
    const label nC = mesh_->nCells();
    const double inv = time_ > 0 ? 1.0 / time_ : 0.0;
    for (auto& it : items_) {
        if (it.nc == 1) {
            for (label c = 0; c < nC; ++c) {
                it.meanS[c] = it.sum[c] * inv;
                if (it.prime2) it.p2S[c] = std::max(it.sum2[c] * inv - it.meanS[c] * it.meanS[c], 0.0);
            }
            mesh_->halo().exchange(it.meanS);
            if (it.prime2) mesh_->halo().exchange(it.p2S);
        } else {
            for (label c = 0; c < nC; ++c) {
                Vec3& m = it.meanV[c];
                for (int i = 0; i < 3; ++i) m[i] = it.sum[3 * c + i] * inv;
                if (it.prime2) {
                    Tensor& t = it.p2V[c];
                    for (int i = 0; i < 3; ++i)
                        for (int j = 0; j < 3; ++j) t(i, j) = it.sum2[9 * c + 3 * i + j] * inv - m[i] * m[j];
                    for (int i = 0; i < 3; ++i) t(i, i) = std::max(t(i, i), 0.0);
                }
            }
            mesh_->halo().exchange(it.meanV);
            if (it.prime2) mesh_->halo().exchange(it.p2V);
        }
    }
}

void FieldAverage::addTo(VtkWriter& w) const {
    for (const auto& it : items_) {
        if (it.nc == 1) {
            w.add(it.name + "Mean", it.meanS);
            if (it.prime2) w.add(it.name + "Prime2Mean", it.p2S);
        } else {
            w.add(it.name + "Mean", it.meanV);
            if (it.prime2) w.add(it.name + "Prime2Mean", it.p2V);
        }
    }
}

void FieldAverage::addTo(TecplotWriter& w) const {
    for (const auto& it : items_) {
        if (it.nc == 1) {
            w.add(it.name + "Mean", it.meanS, 1);
            if (it.prime2) w.add(it.name + "Prime2Mean", it.p2S, 1);
        } else {
            w.add(it.name + "Mean", it.meanV, 3);
            if (it.prime2) w.add(it.name + "Prime2Mean", it.p2V, 9);
        }
    }
}

std::vector<FieldAverage::View> FieldAverage::views() const {
    std::vector<View> v;
    for (const auto& it : items_) {
        if (it.nc == 1) {
            v.push_back({it.name + "Mean", it.meanS.data(), 1});
            if (it.prime2) v.push_back({it.name + "Prime2Mean", it.p2S.data(), 1});
        } else {
            v.push_back({it.name + "Mean", reinterpret_cast<const double*>(it.meanV.data()), 3});
            if (it.prime2) v.push_back({it.name + "Prime2Mean", reinterpret_cast<const double*>(it.p2V.data()), 9});
        }
    }
    return v;
}

void FieldAverage::writeRestart(const std::string& dir) const {
    const Mesh& m = *mesh_;
    if (par::master()) std::filesystem::create_directories(dir);
    par::barrier();
    for (const auto& it : items_) {
        restart::writeCellArray(dir + "/avg_" + it.name + ".bin", m, it.sum.data(), it.nc);
        if (it.prime2) restart::writeCellArray(dir + "/avg_" + it.name + "2.bin", m, it.sum2.data(), it.nc2());
    }
    if (par::master()) {
        std::ofstream os(dir + "/fieldAverage.txt");
        os.precision(17);
        os << "time " << time_ << "\nsamples " << n_ << '\n';
    }
}

bool FieldAverage::readRestart(const std::string& dir) {
    double hdr[2] = {0, 0};
    if (par::master()) {
        std::ifstream is(dir + "/fieldAverage.txt");
        std::string key;
        double v;
        while (is >> key >> v) {
            if (key == "time") hdr[0] = v;
            if (key == "samples") hdr[1] = v;
        }
    }
    par::broadcast(hdr, 2);
    if (!(hdr[0] > 0)) return false;
    const Mesh& m = *mesh_;
    for (auto& it : items_) {
        if (!restart::readCellArray(dir + "/avg_" + it.name + ".bin", m, it.sum.data(), it.nc)) return false;
        if (it.prime2 && !restart::readCellArray(dir + "/avg_" + it.name + "2.bin", m, it.sum2.data(), it.nc2())) return false;
    }
    time_ = hdr[0];
    n_ = int(hdr[1]);
    update();
    return true;
}

} // namespace cfd
