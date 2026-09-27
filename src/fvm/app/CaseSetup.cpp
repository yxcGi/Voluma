#include "fvm/app/CaseSetup.h"

#include "fvm/mesh/AirfoilMesh.h"
#include "fvm/mesh/MeshReaders.h"
#include "fvm/mesh/RawMesh.h"

#include <algorithm>
#include <filesystem>
#include <stdexcept>

namespace cfd {

std::string caseDirectory(const std::string& caseFile) {
    const auto p = std::filesystem::path(caseFile).parent_path();
    return p.empty() ? std::string(".") : p.string();
}

std::string resolvePath(const std::string& caseDir, const std::string& p) {
    if (p.empty() || std::filesystem::path(p).is_absolute()) return p;
    return (std::filesystem::path(caseDir) / p).string();
}

MeshPtr buildMesh(const Json& d, const std::string& caseDir) {
    RawMesh raw;
    if (par::master()) {
        if (d.has("polyMesh")) {
            raw = readPolyMesh(resolvePath(caseDir, d["polyMesh"].string()));
        } else if (d.has("file")) {
            // 外部网格：Gmsh .msh / Fluent .msh/.cas / polyMesh 目录，按内容自动识别
            MeshImportOptions o;
            o.scale = d.get("scale", 1.0);
            o.depth2D = d.get("depth", o.depth2D);
            if (d.has("patchTypes"))
                for (auto& [name, t] : d["patchTypes"].items()) o.patchTypes[name] = patchTypeFromString(t.string());
            const std::string f = resolvePath(caseDir, d["file"].string());
            const std::string fmt = d.get("format", "auto");
            raw = fmt == "gmsh" ? readGmsh(f, o) : fmt == "fluent" ? readFluent(f, o) : readMeshFile(f, o);
            if (d.has("writePolyMesh")) writePolyMesh(raw, resolvePath(caseDir, d["writePolyMesh"].string()));
        } else if (d.has("box")) {
            const Json& b = d["box"];
            BoxSpec s;
            const Json& n = b["n"];
            if (!n.isArray() || n.size() != 3) throw std::runtime_error("mesh.box.n must be [nx, ny, nz]");
            for (int k = 0; k < 3; ++k) s.n[k] = int(n[k].number());
            s.lo = b.get("lo", Vec3{0, 0, 0});
            s.hi = b.get("hi", Vec3{1, 1, 1});
            if (b.has("periodic"))
                for (int k = 0; k < 3; ++k) s.periodic[k] = b["periodic"][k].boolean();
            s.twoD = b.get("twoD", false);
            static const char* dirNames[3] = {"x", "y", "z"};
            if (b.has("stretch"))
                for (int k = 0; k < 3; ++k)
                    if (b["stretch"].has(dirNames[k])) s.stretch[k] = tanhStretch(b["stretch"][dirNames[k]].number());
            static const char* sides[6] = {"xMin", "xMax", "yMin", "yMax", "zMin", "zMax"};
            if (b.has("walls"))
                for (auto& w : b["walls"].array())
                    for (int k = 0; k < 6; ++k)
                        if (w.string() == sides[k]) s.types[k] = PatchType::Wall;
            if (b.has("names"))
                for (int k = 0; k < 6; ++k)
                    if (b["names"].has(sides[k])) s.names[k] = b["names"][sides[k]].string();
            raw = generateBox(s);
        } else if (d.has("airfoil")) {
            const Json& a = d["airfoil"];
            AirfoilMeshSpec s;
            s.naca = a.get("naca", s.naca);
            s.closedTE = a.get("closedTE", s.closedTE);
            s.chord = a.get("chord", s.chord);
            s.alphaDeg = a.get("alpha", s.alphaDeg);
            s.upstream = a.get("upstream", s.upstream);
            s.downstream = a.get("downstream", s.downstream);
            s.halfHeight = a.get("halfHeight", s.halfHeight);
            s.nAirfoil = a.get("nAirfoil", s.nAirfoil);
            s.nWake = a.get("nWake", s.nWake);
            s.nNormal = a.get("nNormal", s.nNormal);
            s.firstCell = a.get("firstCell", s.firstCell);
            s.smoothIter = a.get("smoothIter", s.smoothIter);
            s.depth = a.get("depth", s.depth);
            raw = generateAirfoilCMesh(s);
            if (a.has("writePolyMesh")) writePolyMesh(raw, resolvePath(caseDir, a["writePolyMesh"].string()));
        } else {
            throw std::runtime_error("mesh: need \"polyMesh\", \"file\", \"box\" or \"airfoil\"");
        }
    } else {
        // 其他进程也要把字典标记为已读，避免误报
        for (const char* k : {"file", "format", "scale", "depth", "patchTypes", "writePolyMesh"}) (void)d[k];
        (void)d["polyMesh"];
        (void)d["box"];
        (void)d["airfoil"];
    }
    return Mesh::build(par::master() ? &raw : nullptr);
}

template <class T> void setBoundaryCondition(VolField<T>& f, const std::string& patch, const Json& d) {
    const Mesh& m = f.mesh();
    if (m.isEmptyPatch(patch)) return;
    const std::string type = d.get("type", "");
    if (type == "fixedValue") {
        const T v = d["value"].template as<T>();
        if (d.has("ramp")) {
            // 缓启动：value · s(t/duration)，s 为 smoothstep（6y⁵−15y⁴+10y³，同 OpenLB PolynomialStartScale）或 linear
            const Json& r = d["ramp"];
            const scalar dur = r["duration"].number();
            const bool smooth = r.get("shape", "smoothstep") == "smoothstep";
            if (!(dur > 0)) throw std::runtime_error(f.name() + " on patch " + patch + ": ramp.duration must be > 0");
            f.template set<FixedValueBC<T>>(patch, typename FixedValueBC<T>::Fn([v, dur, smooth](const Vec3&, scalar t) {
                const scalar y = std::clamp(t / dur, scalar(0), scalar(1));
                const scalar s = smooth ? y * y * y * (10 + y * (6 * y - 15)) : y;
                return T(v * s);
            }));
        } else {
            f.template set<FixedValueBC<T>>(patch, v);
        }
    } else if (type == "noSlip") {
        f.template set<FixedValueBC<T>>(patch, Traits<T>::zero());
    } else if (type == "zeroGradient") {
        f.template set<ZeroGradientBC<T>>(patch);
    } else if (type == "fixedGradient") {
        f.template set<FixedGradientBC<T>>(patch, d["gradient"].template as<T>());
    } else if (type == "robin") {
        f.template set<RobinBC<T>>(patch, d["a"].number(), d["b"].number(), d["c"].template as<T>());
    } else if (type == "symmetry" || type == "slip") {
        f.template set<SymmetryBC<T>>(patch);
    } else {
        throw std::runtime_error(f.name() + " on patch " + patch + ": unknown boundary type '" + type + "'");
    }
}

template <class T> void setBoundaryConditions(VolField<T>& f, const Json& boundary, const Json& defaultWall) {
    const Mesh& m = f.mesh();
    for (const auto& p : m.patches()) {
        if (p.type == PatchType::Empty || p.type == PatchType::Cyclic) continue;
        const Json& pd = boundary[p.name];
        if (pd.has(f.name())) {
            setBoundaryCondition(f, p.name, pd[f.name()]);
        } else if (p.type == PatchType::Wall && !defaultWall.isNull()) {
            setBoundaryCondition(f, p.name, defaultWall);
        } else if (p.type == PatchType::Symmetry) {
            f.template set<SymmetryBC<T>>(p.name);
        } else {
            throw std::runtime_error("no boundary condition for " + f.name() + " on patch " + p.name);
        }
    }
    // 其他进程可能没有某些 patch 的面，但 patch 列表全局一致，因此所有进程都会读取字典
}

template void setBoundaryCondition(VolField<scalar>&, const std::string&, const Json&);
template void setBoundaryCondition(VolField<Vec3>&, const std::string&, const Json&);
template void setBoundaryConditions(VolField<scalar>&, const Json&, const Json&);
template void setBoundaryConditions(VolField<Vec3>&, const Json&, const Json&);

SolverControls parseSolverControls(const Json& d, SolverControls c) {
    if (d.isNull()) return c;
    c.solver = d.get("solver", c.solver);
    c.preconditioner = d.get("preconditioner", c.preconditioner);
    c.tolerance = d.get("tolerance", c.tolerance);
    c.relTol = d.get("relTol", c.relTol);
    c.maxIter = d.get("maxIter", c.maxIter);
    c.minIter = d.get("minIter", c.minIter);
    c.nSweeps = d.get("nSweeps", c.nSweeps);
    return c;
}

} // namespace cfd
