// meshCheck：导入并检查网格（类似 OpenFOAM checkMesh），可转换为 polyMesh 与 VTK。
//
//   meshCheck <网格>  [--scale s] [--depth d] [--patch-type name=wall ...] [--polymesh 输出目录] [--vtk 输出目录]
//
// <网格> 可为 Gmsh .msh、Fluent .msh/.cas 或 polyMesh 目录（自动识别）。
// 检查项：单元闭合性 |ΣS_f|/Σ|S_f|、面朝向（owner→neighbour）、单元体积、非正交角、各 patch 面数与面积。
#include "fvm/io/VtkWriter.h"
#include "fvm/mesh/Mesh.h"
#include "fvm/mesh/MeshReaders.h"

#include <cmath>
#include <iostream>
#include <string>

using namespace cfd;

namespace {

Vec3 faceArea(const RawMesh& m, glabel f, Vec3& c) {
    const glabel b = m.faceOffsets[f], e = m.faceOffsets[f + 1];
    Vec3 avg;
    for (glabel k = b; k < e; ++k) avg += m.points[m.facePoints[k]];
    avg /= scalar(e - b);
    Vec3 n;
    for (glabel k = b; k < e; ++k) {
        const Vec3& p = m.points[m.facePoints[k]];
        const Vec3& q = m.points[m.facePoints[k + 1 < e ? k + 1 : b]];
        n += cross(p - avg, q - avg);
    }
    c = avg;
    return 0.5 * n;
}

// 在原始网格上检查闭合性与朝向（只在 0 号进程）
int checkRaw(const RawMesh& m) {
    int bad = 0;
    const glabel nc = m.nCells;
    std::vector<Vec3> sum(nc), cc(nc);
    std::vector<scalar> area(nc, 0);
    std::vector<int> nf(nc, 0);
    std::vector<Vec3> Sf(m.nFaces()), Cf(m.nFaces());
    for (glabel f = 0; f < m.nFaces(); ++f) {
        Sf[f] = faceArea(m, f, Cf[f]);
        const glabel o = m.owner[f];
        sum[o] += Sf[f];
        area[o] += mag(Sf[f]);
        cc[o] += Cf[f];
        ++nf[o];
        if (f < m.nInternalFaces()) {
            const glabel n = m.neighbour[f];
            if (n <= o) ++bad;
            sum[n] -= Sf[f];
            area[n] += mag(Sf[f]);
            cc[n] += Cf[f];
            ++nf[n];
        }
    }
    scalar maxOpen = 0;
    for (glabel c = 0; c < nc; ++c) {
        maxOpen = std::max(maxOpen, mag(sum[c]) / std::max(area[c], VSMALL));
        cc[c] /= scalar(std::max(nf[c], 1));
    }
    glabel wrong = 0;
    for (glabel f = 0; f < m.nFaces(); ++f) {
        const Vec3 d = f < m.nInternalFaces() ? cc[m.neighbour[f]] - cc[m.owner[f]] : Cf[f] - cc[m.owner[f]];
        if (dot(d, Sf[f]) <= 0) ++wrong;
    }
    std::cout << "  cells " << nc << ", faces " << m.nFaces() << " (internal " << m.nInternalFaces() << "), points "
              << m.points.size() << '\n';
    std::cout << "  max cell openness |sum Sf|/sum|Sf| = " << maxOpen << (maxOpen > 1e-6 ? "  ***" : "") << '\n';
    std::cout << "  faces pointing into their owner: " << wrong << (wrong ? "  ***" : "") << '\n';
    if (bad) std::cout << "  internal faces with neighbour <= owner: " << bad << "  ***\n";
    for (const auto& p : m.patches) {
        scalar a = 0;
        for (glabel f = p.start; f < p.start + p.size; ++f) a += mag(Sf[f]);
        std::cout << "  patch " << p.name << " (" << toString(p.type) << "): " << p.size << " faces, area " << a << '\n';
    }
    return (maxOpen > 1e-6 || wrong || bad) ? 1 : 0;
}

} // namespace

int main(int argc, char** argv) {
    par::Environment env(argc, argv);
    if (argc < 2) {
        std::cout << "usage: meshCheck <mesh file or polyMesh dir> [--scale s] [--depth d] [--patch-type name=type]"
                     " [--polymesh dir] [--vtk dir]\n";
        return 1;
    }
    MeshImportOptions opt;
    std::string polyOut, vtkOut;
    for (int i = 2; i + 1 < argc; i += 2) {
        const std::string k = argv[i], v = argv[i + 1];
        if (k == "--scale") opt.scale = std::stod(v);
        else if (k == "--depth") opt.depth2D = std::stod(v);
        else if (k == "--polymesh") polyOut = v;
        else if (k == "--vtk") vtkOut = v;
        else if (k == "--patch-type") {
            const auto e = v.find('=');
            opt.patchTypes[v.substr(0, e)] = patchTypeFromString(v.substr(e + 1));
        } else {
            std::cerr << "unknown option " << k << '\n';
            return 1;
        }
    }
    int status = 0;
    try {
        RawMesh raw;
        if (par::master()) {
            raw = readMeshFile(argv[1], opt);
            std::cout << "Raw mesh " << argv[1] << ":\n";
            status = checkRaw(raw);
            if (!polyOut.empty()) {
                writePolyMesh(raw, polyOut);
                std::cout << "  written polyMesh to " << polyOut << '\n';
            }
        }
        auto m = Mesh::build(par::master() ? &raw : nullptr);
        m->printSummary();
        scalar vmin = GREAT, vmax = 0;
        for (label c = 0; c < m->nCells(); ++c) {
            vmin = std::min(vmin, m->V()[c]);
            vmax = std::max(vmax, m->V()[c]);
        }
        vmin = par::allMin(vmin);
        vmax = par::allMax(vmax);
        if (par::master()) std::cout << "  cell volume min " << vmin << " max " << vmax << (vmin <= 0 ? "  ***" : "") << '\n';
        if (vmin <= 0) status = 1;
        if (!vtkOut.empty()) {
            VtkWriter w(m, vtkOut, "mesh");
            std::vector<scalar> V(m->V().begin(), m->V().begin() + m->nCells());
            w.add("cellVolume", V);
            w.write(0);
        }
    } catch (const std::exception& e) {
        std::cerr << "[rank " << par::rank() << "] error: " << e.what() << '\n';
        par::abort(1);
    }
    if (par::master()) std::cout << (status ? "Mesh check FAILED" : "Mesh OK") << '\n';
    return status;
}
