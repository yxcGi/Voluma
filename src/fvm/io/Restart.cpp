#include "fvm/io/Restart.h"

#include <cstring>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <sstream>
#include <stdexcept>

namespace cfd::restart {

namespace {
constexpr std::uint64_t MAGIC = 0x4346445253543031ull;  // "CFDRST01"

void writeGlobal(const std::string& file, const std::vector<double>& g, glabel n, int nc) {
    if (!par::master()) return;
    std::ofstream os(file, std::ios::binary);
    if (!os) throw std::runtime_error("cannot write " + file);
    const std::uint64_t h[3] = {MAGIC, std::uint64_t(n), std::uint64_t(nc)};
    os.write(reinterpret_cast<const char*>(h), sizeof(h));
    os.write(reinterpret_cast<const char*>(g.data()), std::streamsize(g.size() * sizeof(double)));
}

// 0 号读入全局数组；所有进程返回是否存在
bool readGlobal(const std::string& file, glabel n, int nc, std::vector<double>& g) {
    double ok = 0;
    if (par::master()) {
        std::ifstream is(file, std::ios::binary);
        if (is) {
            std::uint64_t h[3];
            is.read(reinterpret_cast<char*>(h), sizeof(h));
            if (h[0] != MAGIC || h[1] != std::uint64_t(n) || h[2] != std::uint64_t(nc))
                throw std::runtime_error("restart file " + file + " does not match the mesh");
            g.resize(std::size_t(n) * nc);
            is.read(reinterpret_cast<char*>(g.data()), std::streamsize(g.size() * sizeof(double)));
            ok = 1;
        }
    }
    par::broadcast(&ok, 1);
    return ok > 0.5;
}
} // namespace

void writeCellArray(const std::string& file, const Mesh& m, const double* data, int nc) {
    std::vector<double> g;
    const auto& ord = m.cellOrdering();
    // 按分量打包成 nc 个 double 的块
    std::vector<double> local(data, data + std::size_t(m.nCells()) * nc);
    if (nc == 1)
        g = ord.gather(local.data());
    else {
        struct B3 { double v[3]; };
        struct B9 { double v[9]; };
        if (nc == 3) {
            auto r = ord.gather(reinterpret_cast<const B3*>(local.data()));
            g.assign(reinterpret_cast<double*>(r.data()), reinterpret_cast<double*>(r.data()) + r.size() * 3);
        } else if (nc == 9) {
            auto r = ord.gather(reinterpret_cast<const B9*>(local.data()));
            g.assign(reinterpret_cast<double*>(r.data()), reinterpret_cast<double*>(r.data()) + r.size() * 9);
        } else {
            throw std::runtime_error("writeCellArray: unsupported component count");
        }
    }
    writeGlobal(file, g, m.nGlobalCells(), nc);
}

bool readCellArray(const std::string& file, const Mesh& m, double* data, int nc) {
    std::vector<double> g;
    if (!readGlobal(file, m.nGlobalCells(), nc, g)) return false;
    const auto& ord = m.cellOrdering();
    std::vector<double> tmp(std::size_t(m.nTotalCells()) * nc, 0.0);
    if (nc == 1) {
        ord.scatter(g, tmp.data());
    } else if (nc == 3) {
        struct B3 { double v[3]; };
        std::vector<B3> gb(par::master() ? g.size() / 3 : 0);
        if (par::master()) std::memcpy(static_cast<void*>(gb.data()), g.data(), g.size() * sizeof(double));
        ord.scatter(gb, reinterpret_cast<B3*>(tmp.data()));
    } else {
        throw std::runtime_error("readCellArray: unsupported component count");
    }
    std::memcpy(data, tmp.data(), sizeof(double) * std::size_t(m.nCells()) * nc);
    return true;
}

void writeFaceArray(const std::string& file, const Mesh& m, const double* data, int nc) {
    const auto g = m.gatherFaces(data, nc);
    writeGlobal(file, g, m.nGlobalFaces(), nc);
}

bool readFaceArray(const std::string& file, const Mesh& m, double* data, int nc) {
    std::vector<double> g;
    if (!readGlobal(file, m.nGlobalFaces(), nc, g)) return false;
    m.scatterFaces(g, data, nc);
    return true;
}

void write(const std::string& dir, IncompressibleFlow& flow, const std::map<std::string, double>& extra) {
    const Mesh& m = flow.mesh();
    if (par::master()) std::filesystem::create_directories(dir);
    par::barrier();
    auto& U = flow.U();
    auto vd = [](std::vector<Vec3>& v) { return reinterpret_cast<const double*>(v.data()); };
    writeCellArray(dir + "/U.bin", m, vd(U.internal()), 3);
    writeCellArray(dir + "/p.bin", m, flow.p().internal().data(), 1);
    writeFaceArray(dir + "/phi.bin", m, flow.phi().data(), 1);
    const int nOld = U.nOldTimes();
    if (nOld >= 1) {
        writeCellArray(dir + "/U_0.bin", m, vd(U.oldRef()), 3);
        writeFaceArray(dir + "/phi_0.bin", m, flow.phiOld().data(), 1);
    }
    if (nOld >= 2) {
        writeCellArray(dir + "/U_00.bin", m, vd(U.oldOldRef()), 3);
        writeFaceArray(dir + "/phi_00.bin", m, flow.phiOldOld().data(), 1);
    }
    if (!flow.nut().empty()) writeCellArray(dir + "/nut.bin", m, flow.nut().data(), 1);
    if (par::master()) {
        std::ofstream os(dir + "/state.txt");
        os << std::setprecision(17);
        const auto& ts = flow.time();
        os << "time " << ts.time << "\ndt " << ts.dt << "\ndt0 " << ts.dt0 << "\ntimeIndex " << ts.timeIndex
           << "\nnOldTimes " << nOld << "\ngradP " << flow.forcing.gradP << '\n';
        for (auto& [k, v] : extra) os << "extra." << k << ' ' << v << '\n';
    }
    par::barrier();
}

bool read(const std::string& dir, IncompressibleFlow& flow, std::map<std::string, double>* extra) {
    const Mesh& m = flow.mesh();
    // 状态文件（0 号读取后广播）
    std::vector<char> buf;
    if (par::master()) {
        std::ifstream is(dir + "/state.txt");
        if (is) {
            std::stringstream ss;
            ss << is.rdbuf();
            const std::string s = ss.str();
            buf.assign(s.begin(), s.end());
        }
    }
    par::broadcastBytes(buf);
    if (buf.empty()) return false;
    std::map<std::string, double> kv;
    {
        std::istringstream is(std::string(buf.begin(), buf.end()));
        std::string k;
        double v;
        while (is >> k >> v) kv[k] = v;
    }
    auto& U = flow.U();
    auto vd = [](std::vector<Vec3>& v) { return reinterpret_cast<double*>(v.data()); };
    readCellArray(dir + "/U.bin", m, vd(U.internal()), 3);
    readCellArray(dir + "/p.bin", m, flow.p().internal().data(), 1);
    flow.phi().resize(m.nFaces());
    readFaceArray(dir + "/phi.bin", m, flow.phi().data(), 1);
    const int nOld = int(kv["nOldTimes"]);
    if (nOld >= 1) {
        U.oldRef().assign(m.nTotalCells(), Vec3{});
        readCellArray(dir + "/U_0.bin", m, vd(U.oldRef()), 3);
        m.halo().exchange(U.oldRef());
        flow.phiOld().assign(m.nFaces(), 0.0);
        readFaceArray(dir + "/phi_0.bin", m, flow.phiOld().data(), 1);
    }
    if (nOld >= 2) {
        U.oldOldRef().assign(m.nTotalCells(), Vec3{});
        readCellArray(dir + "/U_00.bin", m, vd(U.oldOldRef()), 3);
        m.halo().exchange(U.oldOldRef());
        flow.phiOldOld().assign(m.nFaces(), 0.0);
        readFaceArray(dir + "/phi_00.bin", m, flow.phiOldOld().data(), 1);
    }
    U.setNOldTimes(nOld);
    if (!flow.nut().empty()) {
        readCellArray(dir + "/nut.bin", m, flow.nut().data(), 1);
        m.halo().exchange(flow.nut());
    }
    auto& ts = flow.time();
    ts.time = kv["time"];
    ts.dt = kv["dt"];
    ts.dt0 = kv["dt0"];
    ts.timeIndex = int(kv["timeIndex"]);
    flow.forcing.gradP = kv["gradP"];
    if (extra)
        for (auto& [k, v] : kv)
            if (k.rfind("extra.", 0) == 0) (*extra)[k.substr(6)] = v;
    flow.markRestarted();
    return true;
}

} // namespace cfd::restart
