#include "fvm/app/CaseSetup.h"

#include "fvm/field/Functions.h"
#include "fvm/mesh/AirfoilMesh.h"
#include "fvm/mesh/CylinderMesh.h"
#include "fvm/mesh/PeriodicHill.h"
#include "fvm/mesh/PipeMesh.h"
#include "fvm/mesh/MeshReaders.h"
#include "fvm/mesh/RawMesh.h"

#include <algorithm>
#include <filesystem>
#include <functional>
#include <stdexcept>
#include <type_traits>

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
            s.wakeGrowth = a.get("wakeGrowth", s.wakeGrowth);
            s.smoothIter = a.get("smoothIter", s.smoothIter);
            s.depth = a.get("depth", s.depth);
            raw = generateAirfoilCMesh(s);
            if (a.has("writePolyMesh")) writePolyMesh(raw, resolvePath(caseDir, a["writePolyMesh"].string()));
        } else if (d.has("periodicHill")) {
            const Json& a = d["periodicHill"];
            PeriodicHillSpec s;
            if (a.has("n")) {
                const Json& n = a["n"];
                if (!n.isArray() || n.size() != 3) throw std::runtime_error("mesh.periodicHill.n must be [nx, ny, nz]");
                for (int k = 0; k < 3; ++k) s.n[k] = int(n[k].number());
            }
            s.H = a.get("H", s.H);
            s.length = a.get("length", s.length);
            s.span = a.get("span", s.span);
            s.height = a.get("height", s.height);
            s.stretch = a.get("stretch", s.stretch);
            raw = generatePeriodicHillMesh(s);
            if (a.has("writePolyMesh")) writePolyMesh(raw, resolvePath(caseDir, a["writePolyMesh"].string()));
        } else if (d.has("pipe")) {
            const Json& a = d["pipe"];
            PipeMeshSpec s;
            s.radius = a.get("radius", s.radius);
            s.length = a.get("length", s.length);
            s.nCore = int(a.get("nCore", double(s.nCore)));
            s.nRing = int(a.get("nRing", double(s.nRing)));
            s.nAxial = int(a.get("nAxial", double(s.nAxial)));
            s.coreFrac = a.get("coreFrac", s.coreFrac);
            s.stretch = a.get("stretch", s.stretch);
            s.periodic = a.get("periodic", s.periodic);
            raw = generatePipeMesh(s);
            if (a.has("writePolyMesh")) writePolyMesh(raw, resolvePath(caseDir, a["writePolyMesh"].string()));
        } else if (d.has("cylinder")) {
            const Json& a = d["cylinder"];
            CylinderMeshSpec s;
            s.D = a.get("D", s.D);
            s.upstream = a.get("upstream", s.upstream);
            s.downstream = a.get("downstream", s.downstream);
            s.halfHeight = a.get("halfHeight", s.halfHeight);
            s.nTheta = a.get("nTheta", s.nTheta);
            s.nRadial = a.get("nRadial", s.nRadial);
            s.nSpan = a.get("nSpan", s.nSpan);
            s.span = a.get("span", s.span);
            s.firstCell = a.get("firstCell", s.firstCell);
            s.twoD = a.get("twoD", s.twoD);
            if (s.twoD) s.nSpan = 1;
            raw = generateCylinderOMesh(s);
            if (a.has("writePolyMesh")) writePolyMesh(raw, resolvePath(caseDir, a["writePolyMesh"].string()));
        } else {
            throw std::runtime_error(
                "mesh: need \"polyMesh\", \"file\", \"box\", \"airfoil\", \"periodicHill\", \"pipe\" or \"cylinder\"");
        }
    } else {
        // 其他进程也要把字典标记为已读，避免误报
        for (const char* k : {"file", "format", "scale", "depth", "patchTypes", "writePolyMesh"}) (void)d[k];
        (void)d["polyMesh"];
        (void)d["box"];
        (void)d["airfoil"];
        (void)d["periodicHill"];
        (void)d["pipe"];
        (void)d["cylinder"];
    }
    return Mesh::build(par::master() ? &raw : nullptr);
}

namespace {

Vec3 unitOrThrow(const Json& j, const char* what) {
    const Vec3 v = j.vec3();
    const scalar m = mag(v);
    if (!(m > 0)) throw std::runtime_error(j.where() + ": " + what + " must be non-zero");
    return v / m;
}

// patch 的面积加权平均单位法向与面心（全局）
void patchGeometry(const Mesh& m, label p, Vec3& nMean, Vec3& centre) {
    const Patch& pt = m.patches()[p];
    double a[7] = {0, 0, 0, 0, 0, 0, 0};
    for (label f = pt.start; f < pt.end(); ++f) {
        const scalar s = m.magSf()[f];
        for (int k = 0; k < 3; ++k) {
            a[k] += m.Sf()[f][k];
            a[3 + k] += m.Cf()[f][k] * s;
        }
        a[6] += s;
    }
    par::allSumInPlace(a, 7);
    nMean = {a[0], a[1], a[2]};
    const scalar mn = mag(nMean);
    nMean = mn > 0 ? nMean / mn : Vec3{0, 0, 1};
    centre = a[6] > 0 ? Vec3{a[3], a[4], a[5]} / a[6] : Vec3{};
}

// 标量函数 s(t)：常数或时间表
std::function<scalar(scalar)> scalarFunction(const Json& j, const std::string& caseDir) {
    TimeTable<scalar> tt = TimeTable<scalar>::parse(j, caseDir);
    return [tt](scalar t) { return tt(t); };
}

// 定值类的值：恰好给出 value / timeTable / profile / parabolic 之一，ramp、noise 为修饰
template <class T> struct ValueSpec {
    bool uniform = true;
    T value{};
    typename FixedValueBC<T>::Fn fn;
};

template <class T>
ValueSpec<T> parseValue(const Json& d, const std::vector<std::string>& keys, const VolField<T>& f, const std::string& patch,
                        const BCContext* ctx) {
    using Fn = typename FixedValueBC<T>::Fn;
    const std::string caseDir = ctx ? ctx->caseDir : ".";
    const std::string where = f.name() + " on patch " + patch;
    ValueSpec<T> spec;
    Fn base;
    int n = 0;
    for (const auto& k : keys)
        if (d.has(k)) {
            spec.value = d[k].template as<T>();
            ++n;
        }
    if (d.has("timeTable")) {
        TimeTable<T> tt = TimeTable<T>::parse(d["timeTable"], caseDir);
        base = [tt](const Vec3&, const Vec3&, scalar t) { return tt(t); };
        ++n;
    }
    if (d.has("profile")) {
        ProfileTable<T> pr = ProfileTable<T>::parse(d["profile"], caseDir);
        base = [pr](const Vec3& x, const Vec3&, scalar) { return pr(x); };
        ++n;
    }
    if (d.has("parabolic")) {
        if constexpr (!std::is_same_v<T, Vec3>) {
            throw std::runtime_error(where + ": parabolic profile is only for velocity");
        } else {
            // 抛物线入流：圆管 u = Umax (1 − r²/R²)，Umax = 2 Ubulk；平面槽道 u = Umax (1 − s²/h²)，Umax = 1.5 Ubulk
            const Json& pb = d["parabolic"];
            const Mesh& m = f.mesh();
            const label p = m.findPatch(patch);
            Vec3 nMean, centre;
            patchGeometry(m, p, nMean, centre);
            const Vec3 c = pb.get("center", centre);
            const Vec3 dir = pb.has("direction") ? unitOrThrow(pb["direction"], "parabolic.direction") : -nMean;
            if (pb.has("radius")) {
                const scalar R = pb["radius"].number();
                const scalar Umax = pb.has("Umax") ? pb["Umax"].number() : 2.0 * pb["Ubulk"].number();
                base = [=](const Vec3& x, const Vec3&, scalar) {
                    Vec3 r = x - c;
                    r -= dot(r, nMean) * nMean;
                    return dir * (Umax * std::max(1.0 - magSqr(r) / (R * R), 0.0));
                };
            } else if (pb.has("halfWidth")) {
                const scalar h = pb["halfWidth"].number();
                const Vec3 axis = unitOrThrow(pb["axis"], "parabolic.axis");
                const scalar Umax = pb.has("Umax") ? pb["Umax"].number() : 1.5 * pb["Ubulk"].number();
                base = [=](const Vec3& x, const Vec3&, scalar) {
                    const scalar s = dot(x - c, axis) / h;
                    return dir * (Umax * std::max(1.0 - s * s, 0.0));
                };
            } else {
                throw std::runtime_error(where + ": parabolic needs \"radius\" (pipe) or \"halfWidth\" + \"axis\" (channel)");
            }
            ++n;
        }
    }
    if (n != 1) {
        std::string ks;
        for (const auto& k : keys) ks += (ks.empty() ? "" : "/") + k;
        throw std::runtime_error(where + ": give exactly one of " + ks + ", timeTable, profile, parabolic");
    }
    // 修饰：缓启动 ramp、随机扰动 noise
    if (d.has("ramp")) {
        // value · s(t/duration)，s 为 smoothstep（6y⁵−15y⁴+10y³，同 OpenLB PolynomialStartScale）或 linear
        const Json& r = d["ramp"];
        const scalar dur = r["duration"].number();
        const bool smooth = r.get("shape", "smoothstep") == "smoothstep";
        if (!(dur > 0)) throw std::runtime_error(where + ": ramp.duration must be > 0");
        const T v = spec.value;
        Fn inner = base ? base : Fn([v](const Vec3&, const Vec3&, scalar) { return v; });
        base = [inner, dur, smooth](const Vec3& x, const Vec3& nf, scalar t) {
            const scalar y = std::clamp(t / dur, scalar(0), scalar(1));
            const scalar s = smooth ? y * y * y * (10 + y * (6 * y - 15)) : y;
            return T(inner(x, nf, t) * s);
        };
    }
    if (d.has("noise")) {
        if constexpr (!std::is_same_v<T, Vec3>) {
            throw std::runtime_error(where + ": noise is only for velocity");
        } else {
            // 每面、每步独立的随机扰动，幅值为 noise·|U|（由面心坐标与时间散列，与进程数无关）
            const scalar amp = d["noise"].number();
            const T v = spec.value;
            Fn inner = base ? base : Fn([v](const Vec3&, const Vec3&, scalar) { return v; });
            base = [inner, amp](const Vec3& x, const Vec3& nf, scalar t) {
                Vec3 u = inner(x, nf, t);
                const scalar mu = mag(u);
                for (int k = 0; k < 3; ++k) u[k] += amp * mu * coordHash(x, k, t);
                return u;
            };
        }
    }
    if (base) {
        spec.uniform = false;
        spec.fn = base;
    }
    return spec;
}

template <class T> void setFixedValue(VolField<T>& f, const std::string& patch, const ValueSpec<T>& v) {
    if (v.uniform)
        f.template set<FixedValueBC<T>>(patch, v.value);
    else
        f.template set<FixedValueBC<T>>(patch, v.fn);
}

} // namespace

template <class T> void setBoundaryCondition(VolField<T>& f, const std::string& patch, const Json& d, const BCContext* ctx) {
    const Mesh& m = f.mesh();
    if (m.isEmptyPatch(patch)) return;
    const std::string type = d.get("type", "");
    const std::string where = f.name() + " on patch " + patch;
    const std::string caseDir = ctx ? ctx->caseDir : ".";
    constexpr bool isVec = std::is_same_v<T, Vec3>;
    constexpr bool isScalar = std::is_same_v<T, scalar>;
    if (type == "fixedValue") {
        setFixedValue(f, patch, parseValue(d, {"value"}, f, patch, ctx));
    } else if (type == "noSlip") {
        f.template set<FixedValueBC<T>>(patch, Traits<T>::zero());
    } else if (type == "movingWall") {
        if constexpr (!isVec) throw std::runtime_error(where + ": movingWall is only for velocity");
        setFixedValue(f, patch, parseValue(d, {"velocity", "value"}, f, patch, ctx));
    } else if (type == "rotatingWall") {
        if constexpr (!isVec) {
            throw std::runtime_error(where + ": rotatingWall is only for velocity");
        } else {
            // U = ω (axis × (x − origin))，去掉法向分量（壁面不完全圆时）
            const Vec3 origin = d.get("origin", Vec3{0, 0, 0});
            const Vec3 axis = unitOrThrow(d["axis"], "rotatingWall.axis");
            auto omega = scalarFunction(d["omega"], caseDir);
            f.template set<FixedValueBC<Vec3>>(patch, typename FixedValueBC<Vec3>::Fn([=](const Vec3& x, const Vec3& n, scalar t) {
                const Vec3 u = omega(t) * cross(axis, x - origin);
                return u - dot(u, n) * n;
            }));
        }
    } else if (type == "flowRateInlet") {
        if constexpr (!isVec) {
            throw std::runtime_error(where + ": flowRateInlet is only for velocity");
        } else {
            std::function<scalar(scalar)> Q;
            if (d.has("volumetricFlowRate")) {
                Q = scalarFunction(d["volumetricFlowRate"], caseDir);
            } else if (d.has("massFlowRate")) {
                auto mdot = scalarFunction(d["massFlowRate"], caseDir);
                const scalar rho = d.get("rho", 1.0);
                Q = [mdot, rho](scalar t) { return mdot(t) / rho; };
            } else {
                throw std::runtime_error(where + ": flowRateInlet needs volumetricFlowRate or massFlowRate");
            }
            f.template set<FlowRateInletBC>(patch, Q);
        }
    } else if (type == "zeroGradient") {
        f.template set<ZeroGradientBC<T>>(patch);
    } else if (type == "fixedGradient") {
        f.template set<FixedGradientBC<T>>(patch, d["gradient"].template as<T>());
    } else if (type == "fixedFluxPressure") {
        if constexpr (!isScalar) throw std::runtime_error(where + ": fixedFluxPressure is only for pressure");
        else f.template set<FixedFluxPressureBC>(patch);
    } else if (type == "robin") {
        f.template set<RobinBC<T>>(patch, d["a"].number(), d["b"].number(), d["c"].template as<T>());
    } else if (type == "symmetry" || type == "slip") {
        f.template set<SymmetryBC<T>>(patch);
    } else if (type == "inletOutlet") {
        if (!ctx) throw std::runtime_error(where + ": inletOutlet needs the flow context");
        Json zero = Json::parse(isVec ? "{\"inletValue\": [0, 0, 0]}" : "{\"inletValue\": 0}");
        const bool given = d.has("inletValue") || d.has("value") || d.has("timeTable") || d.has("profile");
        const ValueSpec<T> v = parseValue(given ? d : zero, {"inletValue", "value"}, f, patch, ctx);
        if (v.uniform)
            f.template set<InletOutletBC<T>>(patch, ctx, v.value);
        else
            f.template set<InletOutletBC<T>>(patch, ctx, v.fn);
    } else if (type == "freestream") {
        const ValueSpec<T> v = parseValue(d, {"freestreamValue", "value"}, f, patch, ctx);
        if constexpr (isVec) {
            if (!v.uniform) throw std::runtime_error(where + ": freestream velocity must be uniform");
            f.template set<FreestreamVelocityBC>(patch, v.value);
        } else {
            if (!ctx) throw std::runtime_error(where + ": freestream needs the flow context");
            if (v.uniform)
                f.template set<InletOutletBC<T>>(patch, ctx, v.value);
            else
                f.template set<InletOutletBC<T>>(patch, ctx, v.fn);
        }
    } else if (type == "freestreamPressure") {
        if constexpr (!isScalar) {
            throw std::runtime_error(where + ": freestreamPressure is only for pressure");
        } else {
            const scalar pInf = d.has("freestreamValue") ? d["freestreamValue"].number() : d.get("value", 0.0);
            if (d.has("U")) {
                const Vec3 U = d["U"].vec3();
                f.template set<FreestreamPressureBC>(patch, ctx, pInf, &U);
            } else {
                if (!ctx) throw std::runtime_error(where + ": freestreamPressure needs \"U\" or the flow context");
                f.template set<FreestreamPressureBC>(patch, ctx, pInf, nullptr);
            }
        }
    } else if (type == "totalPressure") {
        if constexpr (!isScalar) {
            throw std::runtime_error(where + ": totalPressure is only for pressure");
        } else {
            if (!ctx) throw std::runtime_error(where + ": totalPressure needs the flow context");
            f.template set<TotalPressureBC>(patch, ctx, scalarFunction(d["p0"], caseDir));
        }
    } else if (type == "pressureInletOutletVelocity") {
        if constexpr (!isVec) {
            throw std::runtime_error(where + ": pressureInletOutletVelocity is only for velocity");
        } else {
            if (!ctx) throw std::runtime_error(where + ": pressureInletOutletVelocity needs the flow context");
            f.template set<PressureInletOutletVelocityBC>(patch, ctx, d.get("tangentialVelocity", Vec3{0, 0, 0}));
        }
    } else if (type == "advective") {
        if (!ctx) throw std::runtime_error(where + ": advective needs the flow context");
        const bool haveInf = d.has("fieldInf");
        const T fieldInf = haveInf ? d["fieldInf"].template as<T>() : Traits<T>::zero();
        const scalar lInf = d.get("lInf", 0.0);
        if (haveInf && !(lInf > 0)) throw std::runtime_error(where + ": advective with fieldInf needs lInf > 0");
        f.template set<AdvectiveBC<T>>(patch, ctx, haveInf, fieldInf, lInf);
    } else {
        throw std::runtime_error(where + ": unknown boundary type '" + type + "'");
    }
}

template <class T>
void setBoundaryConditions(VolField<T>& f, const Json& boundary, const Json& defaultWall, const BCContext* ctx) {
    const Mesh& m = f.mesh();
    for (const auto& p : m.patches()) {
        if (p.type == PatchType::Empty || p.type == PatchType::Cyclic) continue;
        const Json& pd = boundary[p.name];
        if (pd.has(f.name())) {
            setBoundaryCondition(f, p.name, pd[f.name()], ctx);
        } else if (p.type == PatchType::Wall && !defaultWall.isNull()) {
            setBoundaryCondition(f, p.name, defaultWall, ctx);
        } else if (p.type == PatchType::Symmetry) {
            f.template set<SymmetryBC<T>>(p.name);
        } else {
            throw std::runtime_error("no boundary condition for " + f.name() + " on patch " + p.name);
        }
    }
    // 其他进程可能没有某些 patch 的面，但 patch 列表全局一致，因此所有进程都会读取字典
}

template void setBoundaryCondition(VolField<scalar>&, const std::string&, const Json&, const BCContext*);
template void setBoundaryCondition(VolField<Vec3>&, const std::string&, const Json&, const BCContext*);
template void setBoundaryConditions(VolField<scalar>&, const Json&, const Json&, const BCContext*);
template void setBoundaryConditions(VolField<Vec3>&, const Json&, const Json&, const BCContext*);

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
