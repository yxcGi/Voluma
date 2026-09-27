#include "fvm/app/CaseSetup.h"

#include "fvm/mesh/RawMesh.h"

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
        } else {
            throw std::runtime_error("mesh: need \"polyMesh\" or \"box\"");
        }
    } else {
        // 其他进程也要把字典标记为已读，避免误报
        (void)d["polyMesh"];
        (void)d["box"];
    }
    return Mesh::build(par::master() ? &raw : nullptr);
}

template <class T> void setBoundaryCondition(VolField<T>& f, const std::string& patch, const Json& d) {
    const Mesh& m = f.mesh();
    if (m.isEmptyPatch(patch)) return;
    const std::string type = d.get("type", "");
    if (type == "fixedValue") {
        f.template set<FixedValueBC<T>>(patch, d["value"].template as<T>());
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
