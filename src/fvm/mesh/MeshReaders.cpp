#include "fvm/mesh/MeshReaders.h"

#include "fvm/mesh/Extrude2D.h"

#include <algorithm>
#include <array>
#include <cctype>
#include <cmath>
#include <cstring>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <unordered_map>

namespace cfd {

namespace {

std::string lower(std::string s) {
    for (auto& c : s) c = char(std::tolower(static_cast<unsigned char>(c)));
    return s;
}

// 面节点排序后的键（最多 4 个角点的单元面；多边形面用 vector 键）
struct FaceKey {
    std::array<glabel, 4> v{-1, -1, -1, -1};
    bool operator==(const FaceKey& o) const { return v == o.v; }
};
struct FaceKeyHash {
    std::size_t operator()(const FaceKey& k) const {
        std::size_t h = 1469598103934665603ull;
        for (glabel x : k.v) h = (h ^ std::size_t(x + 1)) * 1099511628211ull;
        return h;
    }
};
FaceKey makeKey(const std::vector<glabel>& n) {
    if (n.size() > 4) throw std::runtime_error("mesh import: element face with more than 4 nodes");
    FaceKey k;
    std::copy(n.begin(), n.end(), k.v.begin());
    std::sort(k.v.begin(), k.v.begin() + n.size());
    return k;
}

Vec3 polygonNormal(const std::vector<Vec3>& P, const std::vector<glabel>& f, Vec3& centre) {
    Vec3 c;
    for (glabel i : f) c += P[i];
    c /= scalar(f.size());
    Vec3 n;
    for (std::size_t k = 0; k < f.size(); ++k) n += cross(P[f[k]] - c, P[f[(k + 1) % f.size()]] - c);
    centre = c;
    return 0.5 * n;
}

// 类型默认值：名称含 wall → wall，含 symmetry → symmetry；可被 opt.patchTypes 覆盖
PatchType defaultPatchType(const std::string& name, PatchType t, const MeshImportOptions& opt) {
    auto it = opt.patchTypes.find(name);
    if (it != opt.patchTypes.end()) return it->second;
    return t;
}

// 二维多边形网格（任意朝向）+ 边界边 → 拉伸一层
RawMesh extrudeImported(Mesh2D m2, std::vector<RawPatch> patches,
                        const std::map<std::pair<glabel, glabel>, int>& edgePatch, const MeshImportOptions& opt) {
    for (auto& c : m2.cells) {
        scalar a = 0;
        for (std::size_t k = 0; k < c.size(); ++k) {
            const Vec3& p = m2.points[c[k]];
            const Vec3& q = m2.points[c[(k + 1) % c.size()]];
            a += p.x * q.y - q.x * p.y;
        }
        if (a < 0) std::reverse(c.begin(), c.end());
    }
    std::vector<BoundaryPatchSpec> specs;
    for (auto& p : patches) specs.push_back({p.name, defaultPatchType(p.name, p.type, opt)});
    const int defaultIdx = int(specs.size());
    specs.push_back({"defaultFaces", defaultPatchType("defaultFaces", PatchType::Wall, opt)});
    glabel nDefault = 0;
    RawMesh r = extrude2D(
        m2, specs,
        [&](glabel a, glabel b) {
            auto it = edgePatch.find({std::min(a, b), std::max(a, b)});
            if (it != edgePatch.end()) return it->second;
            ++nDefault;
            return defaultIdx;
        },
        opt.depth2D);
    if (nDefault) std::cout << "mesh import: " << nDefault << " boundary edges without a group -> patch defaultFaces\n";
    // 去掉空 patch
    std::vector<RawPatch> kept;
    for (auto& p : r.patches)
        if (p.size > 0) kept.push_back(p);
    r.patches = kept;
    return r;
}

// 按边列表重建二维单元多边形（Fluent 二维网格）
std::vector<std::vector<glabel>> polygonsFromEdges(glabel nCells, const std::vector<std::array<glabel, 3>>& cellEdges) {
    // cellEdges: (cell, a, b)
    std::vector<std::vector<std::pair<glabel, glabel>>> e(nCells);
    for (auto& ce : cellEdges) e[ce[0]].push_back({ce[1], ce[2]});
    std::vector<std::vector<glabel>> cells(nCells);
    for (glabel c = 0; c < nCells; ++c) {
        auto& E = e[c];
        if (E.size() < 3) throw std::runtime_error("mesh import: 2D cell " + std::to_string(c) + " has < 3 edges");
        std::vector<bool> used(E.size(), false);
        std::vector<glabel> poly{E[0].first, E[0].second};
        used[0] = true;
        for (std::size_t k = 1; k < E.size(); ++k) {
            const glabel last = poly.back();
            bool found = false;
            for (std::size_t j = 0; j < E.size() && !found; ++j) {
                if (used[j]) continue;
                if (E[j].first == last) poly.push_back(E[j].second), found = true;
                else if (E[j].second == last) poly.push_back(E[j].first), found = true;
                if (found) used[j] = true;
            }
            if (!found) throw std::runtime_error("mesh import: 2D cell " + std::to_string(c) + " edges do not form a loop");
        }
        if (poly.back() != poly.front()) throw std::runtime_error("mesh import: 2D cell " + std::to_string(c) + " is not closed");
        poly.pop_back();
        cells[c] = poly;
    }
    return cells;
}

} // namespace

// ============================================================ 公共装配
RawMesh assembleFaceMesh(std::vector<Vec3> points, glabel nCells, std::vector<FaceRecord> faces,
                         const std::vector<RawPatch>& patches) {
    // 单元近似中心：各面中心的平均（凸单元在内部，足以判定朝向）
    std::vector<Vec3> cc(nCells);
    std::vector<int> nf(nCells, 0);
    std::vector<Vec3> fc(faces.size()), fn(faces.size());
    for (std::size_t f = 0; f < faces.size(); ++f) {
        auto& F = faces[f];
        if (F.c0 < 0 && F.c1 >= 0) std::swap(F.c0, F.c1);
        if (F.c0 < 0 || F.c0 >= nCells || F.c1 >= nCells) throw std::runtime_error("mesh import: face with invalid cells");
        fn[f] = polygonNormal(points, F.nodes, fc[f]);
        cc[F.c0] += fc[f];
        ++nf[F.c0];
        if (F.c1 >= 0) {
            cc[F.c1] += fc[f];
            ++nf[F.c1];
        }
    }
    for (glabel c = 0; c < nCells; ++c) {
        if (nf[c] < 4) throw std::runtime_error("mesh import: cell " + std::to_string(c) + " has fewer than 4 faces");
        cc[c] /= scalar(nf[c]);
    }
    std::vector<std::size_t> internal;
    std::vector<std::vector<std::size_t>> bnd(patches.size());
    for (std::size_t f = 0; f < faces.size(); ++f) {
        auto& F = faces[f];
        if (F.c1 >= 0) {
            if (F.c1 < F.c0) std::swap(F.c0, F.c1);
            if (dot(fn[f], cc[F.c1] - cc[F.c0]) < 0) std::reverse(F.nodes.begin(), F.nodes.end());
            internal.push_back(f);
        } else {
            if (F.patch < 0 || F.patch >= int(patches.size())) throw std::runtime_error("mesh import: boundary face without patch");
            if (dot(fn[f], fc[f] - cc[F.c0]) < 0) std::reverse(F.nodes.begin(), F.nodes.end());
            bnd[F.patch].push_back(f);
        }
    }
    std::sort(internal.begin(), internal.end(), [&](std::size_t a, std::size_t b) {
        return std::tie(faces[a].c0, faces[a].c1) < std::tie(faces[b].c0, faces[b].c1);
    });
    RawMesh r;
    r.points = std::move(points);
    r.nCells = nCells;
    for (auto f : internal) {
        r.addFace(faces[f].nodes);
        r.owner.push_back(faces[f].c0);
        r.neighbour.push_back(faces[f].c1);
    }
    for (std::size_t p = 0; p < patches.size(); ++p) {
        if (bnd[p].empty()) continue;
        RawPatch rp = patches[p];
        rp.start = r.nFaces();
        for (auto f : bnd[p]) {
            r.addFace(faces[f].nodes);
            r.owner.push_back(faces[f].c0);
        }
        rp.size = r.nFaces() - rp.start;
        r.patches.push_back(rp);
    }
    return r;
}

// ============================================================ Gmsh
namespace {

struct GmshType {
    int dim;
    int corners;  // 角点数（高阶单元只取前 corners 个节点）
};
bool gmshType(int t, GmshType& g) {
    switch (t) {
    case 15: g = {0, 1}; return true;
    case 1: case 8: case 26: case 27: case 28: g = {1, 2}; return true;
    case 2: case 9: case 20: case 21: case 22: case 23: case 24: case 25: g = {2, 3}; return true;
    case 3: case 10: case 16: case 36: case 37: case 38: g = {2, 4}; return true;
    case 4: case 11: case 29: case 30: case 31: g = {3, 4}; return true;
    case 5: case 12: case 17: case 92: case 93: g = {3, 8}; return true;
    case 6: case 13: case 18: case 90: case 91: g = {3, 6}; return true;
    case 7: case 14: case 19: case 118: case 119: g = {3, 5}; return true;
    default: return false;
    }
}

// 单元面（Gmsh 节点编号，朝向随后按几何修正）
const std::vector<std::vector<int>>& elementFaces(int corners) {
    static const std::vector<std::vector<int>> tet{{0, 2, 1}, {0, 1, 3}, {0, 3, 2}, {1, 2, 3}};
    static const std::vector<std::vector<int>> pyr{{0, 3, 2, 1}, {0, 1, 4}, {1, 2, 4}, {2, 3, 4}, {3, 0, 4}};
    static const std::vector<std::vector<int>> prism{{0, 2, 1}, {3, 4, 5}, {0, 1, 4, 3}, {1, 2, 5, 4}, {2, 0, 3, 5}};
    static const std::vector<std::vector<int>> hex{{0, 3, 2, 1}, {4, 5, 6, 7}, {0, 1, 5, 4},
                                                   {1, 2, 6, 5}, {2, 3, 7, 6}, {3, 0, 4, 7}};
    switch (corners) {
    case 4: return tet;
    case 5: return pyr;
    case 6: return prism;
    default: return hex;
    }
}

struct GmshElement {
    int dim;
    int phys;  // 物理组编号（无则 -1）
    std::vector<glabel> nodes;  // 节点编号（Gmsh tag）
};

} // namespace

RawMesh readGmsh(const std::string& file, const MeshImportOptions& opt) {
    std::ifstream in(file);
    if (!in) throw std::runtime_error("cannot open " + file);
    double version = 0;
    std::map<std::pair<int, int>, std::string> physNames;          // (dim, tag) → name
    std::map<std::pair<int, int>, std::vector<int>> entityPhys;    // v4：(dim, entity) → 物理组
    std::unordered_map<glabel, glabel> nodeIndex;
    std::vector<Vec3> points;
    std::vector<GmshElement> elems;
    std::string line;
    auto endSection = [&](const std::string& name) {
        while (std::getline(in, line))
            if (line.rfind("$End" + name, 0) == 0) return;
        throw std::runtime_error(file + ": missing $End" + name);
    };
    while (std::getline(in, line)) {
        if (!line.empty() && line.back() == '\r') line.pop_back();
        if (line == "$MeshFormat") {
            int ft = 0, ds = 0;
            in >> version >> ft >> ds;
            if (ft != 0) throw std::runtime_error(file + ": binary Gmsh files are not supported; save as ASCII (Mesh.Binary=0)");
            if (version < 2 || version >= 5) throw std::runtime_error(file + ": unsupported Gmsh version " + std::to_string(version));
            endSection("MeshFormat");
        } else if (line == "$PhysicalNames") {
            int n;
            in >> n;
            for (int i = 0; i < n; ++i) {
                int d, t;
                in >> d >> t;
                std::getline(in, line);
                const auto a = line.find('"'), b = line.rfind('"');
                physNames[{d, t}] = (a != std::string::npos && b > a) ? line.substr(a + 1, b - a - 1) : std::to_string(t);
            }
            endSection("PhysicalNames");
        } else if (line == "$Entities" && version >= 4) {
            std::size_t n[4];
            in >> n[0] >> n[1] >> n[2] >> n[3];
            for (int d = 0; d < 4; ++d)
                for (std::size_t i = 0; i < n[d]; ++i) {
                    int tag;
                    double x;
                    in >> tag;
                    for (int k = 0; k < (d == 0 ? 3 : 6); ++k) in >> x;
                    std::size_t np;
                    in >> np;
                    std::vector<int> ph(np);
                    for (auto& p : ph) in >> p;
                    for (auto& p : ph) p = std::abs(p);
                    entityPhys[{d, tag}] = ph;
                    if (d > 0) {
                        std::size_t nb;
                        in >> nb;
                        for (std::size_t k = 0; k < nb; ++k) in >> tag;
                    }
                }
            endSection("Entities");
        } else if (line == "$Nodes") {
            if (version >= 4) {
                std::size_t nBlocks, nNodes, minT, maxT;
                in >> nBlocks >> nNodes >> minT >> maxT;
                points.reserve(nNodes);
                for (std::size_t b = 0; b < nBlocks; ++b) {
                    int d, tag, param;
                    std::size_t n;
                    in >> d >> tag >> param >> n;
                    std::vector<glabel> tags(n);
                    for (auto& t : tags) in >> t;
                    for (std::size_t i = 0; i < n; ++i) {
                        Vec3 p;
                        in >> p.x >> p.y >> p.z;
                        double u;
                        for (int k = 0; k < (param ? d : 0); ++k) in >> u;
                        nodeIndex[tags[i]] = glabel(points.size());
                        points.push_back(p * opt.scale);
                    }
                }
            } else {
                std::size_t n;
                in >> n;
                points.reserve(n);
                for (std::size_t i = 0; i < n; ++i) {
                    glabel t;
                    Vec3 p;
                    in >> t >> p.x >> p.y >> p.z;
                    nodeIndex[t] = glabel(points.size());
                    points.push_back(p * opt.scale);
                }
            }
            endSection("Nodes");
        } else if (line == "$Elements") {
            if (version >= 4) {
                std::size_t nBlocks, nEl, minT, maxT;
                in >> nBlocks >> nEl >> minT >> maxT;
                elems.reserve(nEl);
                std::getline(in, line);
                for (std::size_t b = 0; b < nBlocks; ++b) {
                    int d, tag, type;
                    std::size_t n;
                    in >> d >> tag >> type >> n;
                    std::getline(in, line);
                    GmshType gt;
                    const bool known = gmshType(type, gt);
                    if (!known) throw std::runtime_error(file + ": unsupported Gmsh element type " + std::to_string(type));
                    auto it = entityPhys.find({d, tag});
                    const int phys = (it != entityPhys.end() && !it->second.empty()) ? it->second[0] : -1;
                    for (std::size_t i = 0; i < n; ++i) {
                        std::getline(in, line);
                        std::istringstream ls(line);
                        glabel id, v;
                        ls >> id;
                        GmshElement e{gt.dim, phys, {}};
                        for (int k = 0; k < gt.corners && ls >> v; ++k) e.nodes.push_back(v);
                        elems.push_back(std::move(e));
                    }
                }
            } else {
                std::size_t n;
                in >> n;
                std::getline(in, line);
                elems.reserve(n);
                for (std::size_t i = 0; i < n; ++i) {
                    std::getline(in, line);
                    std::istringstream ls(line);
                    glabel id;
                    int type, nt;
                    ls >> id >> type >> nt;
                    std::vector<int> tags(nt);
                    for (auto& t : tags) ls >> t;
                    GmshType gt;
                    if (!gmshType(type, gt)) throw std::runtime_error(file + ": unsupported Gmsh element type " + std::to_string(type));
                    GmshElement e{gt.dim, nt > 0 ? tags[0] : -1, {}};
                    if (e.phys == 0) e.phys = -1;
                    glabel v;
                    for (int k = 0; k < gt.corners && ls >> v; ++k) e.nodes.push_back(v);
                    elems.push_back(std::move(e));
                }
            }
            endSection("Elements");
        } else if (line == "$Periodic") {
            std::cout << "mesh import: Gmsh $Periodic section ignored (periodic boundaries are imported as ordinary patches)\n";
        }
    }
    if (points.empty() || elems.empty()) throw std::runtime_error(file + ": no nodes or elements");
    auto node = [&](glabel t) {
        auto it = nodeIndex.find(t);
        if (it == nodeIndex.end()) throw std::runtime_error(file + ": element refers to unknown node " + std::to_string(t));
        return it->second;
    };
    int meshDim = 0;
    for (auto& e : elems) meshDim = std::max(meshDim, e.dim);
    if (meshDim < 2) throw std::runtime_error(file + ": no 2D or 3D elements");

    // patches：维度 meshDim−1 的物理组
    std::map<int, int> physToPatch;
    std::vector<RawPatch> patches;
    for (auto& e : elems)
        if (e.dim == meshDim - 1 && e.phys >= 0 && !physToPatch.count(e.phys)) {
            RawPatch p;
            auto it = physNames.find({meshDim - 1, e.phys});
            p.name = it != physNames.end() ? it->second : "patch" + std::to_string(e.phys);
            p.type = lower(p.name).find("wall") != std::string::npos ? PatchType::Wall : PatchType::Patch;
            if (lower(p.name).find("symmetry") != std::string::npos) p.type = PatchType::Symmetry;
            physToPatch[e.phys] = int(patches.size());
            patches.push_back(p);
        }
    // 同名物理组合并
    for (std::size_t a = 0; a < patches.size(); ++a)
        for (std::size_t b = a + 1; b < patches.size(); ++b)
            if (patches[a].name == patches[b].name) throw std::runtime_error(file + ": duplicate physical name " + patches[a].name);

    if (meshDim == 2) {
        Mesh2D m2;
        m2.points = points;
        const scalar z0 = points[0].z;
        for (auto& p : points)
            if (std::abs(p.z - z0) > 1e-9 * (1 + std::abs(z0)))
                throw std::runtime_error(file + ": 2D Gmsh mesh must lie in a plane z = const");
        std::map<std::pair<glabel, glabel>, int> edgePatch;
        for (auto& e : elems) {
            if (e.dim == 2) {
                std::vector<glabel> c;
                for (auto t : e.nodes) c.push_back(node(t));
                m2.cells.push_back(c);
            } else if (e.dim == 1 && e.phys >= 0) {
                const glabel a = node(e.nodes[0]), b = node(e.nodes[1]);
                edgePatch[{std::min(a, b), std::max(a, b)}] = physToPatch[e.phys];
            }
        }
        for (auto& p : patches) p.type = defaultPatchType(p.name, p.type, opt);
        return extrudeImported(std::move(m2), patches, edgePatch, opt);
    }

    // 三维：由单元生成面并配对
    std::unordered_map<FaceKey, int, FaceKeyHash> bFacePatch;
    for (auto& e : elems)
        if (e.dim == 2 && e.phys >= 0) {
            std::vector<glabel> f;
            for (auto t : e.nodes) f.push_back(node(t));
            bFacePatch[makeKey(f)] = physToPatch[e.phys];
        }
    const int defaultIdx = int(patches.size());
    patches.push_back({"defaultFaces", PatchType::Wall, 0, 0, ""});
    for (auto& p : patches) p.type = defaultPatchType(p.name, p.type, opt);

    std::vector<FaceRecord> faces;
    std::unordered_map<FaceKey, std::size_t, FaceKeyHash> faceOf;
    glabel nCells = 0;
    for (auto& e : elems) {
        if (e.dim != 3) continue;
        const glabel c = nCells++;
        std::vector<glabel> cn;
        for (auto t : e.nodes) cn.push_back(node(t));
        for (auto& fd : elementFaces(int(cn.size()))) {
            std::vector<glabel> f;
            for (int k : fd) f.push_back(cn[k]);
            const FaceKey k = makeKey(f);
            auto it = faceOf.find(k);
            if (it == faceOf.end()) {
                faceOf.emplace(k, faces.size());
                faces.push_back({f, c, -1, -1});
            } else {
                auto& F = faces[it->second];
                if (F.c1 >= 0) throw std::runtime_error(file + ": face shared by more than two cells");
                F.c1 = c;
            }
        }
    }
    glabel nDefault = 0;
    for (auto& F : faces)
        if (F.c1 < 0) {
            auto it = bFacePatch.find(makeKey(F.nodes));
            F.patch = it != bFacePatch.end() ? it->second : (++nDefault, defaultIdx);
        }
    if (nDefault) std::cout << "mesh import: " << nDefault << " boundary faces without a physical group -> patch defaultFaces\n";
    // 去掉未被使用的点
    std::vector<glabel> used(points.size(), -1);
    std::vector<Vec3> pts;
    for (auto& F : faces)
        for (auto& v : F.nodes) {
            if (used[v] < 0) {
                used[v] = glabel(pts.size());
                pts.push_back(points[v]);
            }
            v = used[v];
        }
    return assembleFaceMesh(std::move(pts), nCells, std::move(faces), patches);
}

// ============================================================ Fluent
namespace {

class FluentReader {
public:
    explicit FluentReader(const std::string& file) : file_(file) {
        std::ifstream in(file, std::ios::binary);
        if (!in) throw std::runtime_error("cannot open " + file);
        std::stringstream ss;
        ss << in.rdbuf();
        s_ = ss.str();
    }

    int dim = 3;
    std::vector<Vec3> points;
    glabel nCells = 0;
    struct Zone {
        int id;
        int bcType;
        std::vector<FaceRecord> faces;  // 节点/单元为 0 基
    };
    std::vector<Zone> faceZones;
    std::map<int, std::pair<std::string, std::string>> zoneInfo;  // id → (type, name)

    void parse() {
        while (true) {
            skipWs();
            if (p_ >= s_.size()) break;
            if (s_[p_] != '(') {
                ++p_;
                continue;
            }
            ++p_;
            const int idx = int(readInt(10));
            section(idx);
        }
    }

private:
    std::string file_, s_;
    std::size_t p_ = 0;

    [[noreturn]] void fail(const std::string& m) const {
        throw std::runtime_error(file_ + ": " + m + " (at byte " + std::to_string(p_) + ")");
    }
    void skipWs() {
        while (p_ < s_.size() && std::isspace(static_cast<unsigned char>(s_[p_]))) ++p_;
    }
    void expect(char c) {
        skipWs();
        if (p_ >= s_.size() || s_[p_] != c) fail(std::string("expected '") + c + "'");
        ++p_;
    }
    bool peek(char c) {
        skipWs();
        return p_ < s_.size() && s_[p_] == c;
    }
    long long readInt(int base) {
        skipWs();
        const char* b = s_.c_str() + p_;
        char* e;
        const long long v = std::strtoll(b, &e, base);
        if (e == b) fail("expected integer");
        p_ += std::size_t(e - b);
        return v;
    }
    double readDouble() {
        skipWs();
        const char* b = s_.c_str() + p_;
        char* e;
        const double v = std::strtod(b, &e);
        if (e == b) fail("expected number");
        p_ += std::size_t(e - b);
        return v;
    }
    std::string readWord() {
        skipWs();
        const std::size_t b = p_;
        if (p_ < s_.size() && s_[p_] == '"') {
            const auto e = s_.find('"', p_ + 1);
            p_ = e + 1;
            return s_.substr(b + 1, e - b - 1);
        }
        while (p_ < s_.size() && !std::isspace(static_cast<unsigned char>(s_[p_])) && s_[p_] != '(' && s_[p_] != ')') ++p_;
        return s_.substr(b, p_ - b);
    }
    // 跳过括号平衡的 ASCII 内容（当前位置在 '(' 之后，结束于匹配的 ')' 之后）
    void skipBalanced() {
        int depth = 1;
        while (p_ < s_.size() && depth > 0) {
            const char c = s_[p_++];
            if (c == '"') {
                const auto e = s_.find('"', p_);
                p_ = e == std::string::npos ? s_.size() : e + 1;
            } else if (c == '(') {
                ++depth;
            } else if (c == ')') {
                --depth;
            }
        }
    }
    // 二进制段结束：") End of Binary Section xxxx)"
    void endBinary() {
        const auto e = s_.find("End of Binary Section", p_);
        if (e == std::string::npos) fail("missing 'End of Binary Section'");
        p_ = s_.find(')', e);
        if (p_ == std::string::npos) fail("unterminated binary section");
        ++p_;
    }
    template <class T> T bin() {
        T v;
        if (p_ + sizeof(T) > s_.size()) fail("unexpected end of binary data");
        std::memcpy(&v, s_.data() + p_, sizeof(T));
        p_ += sizeof(T);
        return v;
    }
    // 头部 (a b c d e)，十六进制
    std::vector<long long> header() {
        expect('(');
        std::vector<long long> h;
        while (!peek(')')) h.push_back(readInt(16));
        ++p_;
        return h;
    }

    void section(int idx) {
        const int kind = idx % 1000;
        const int bin = idx / 1000;  // 0 ASCII；2 单精度/32 位；3 双精度/64 位
        switch (kind) {
        case 0: skipBalanced(); return;  // 注释
        case 1: skipBalanced(); return;  // 头
        case 2:
            dim = int(readInt(10));
            expect(')');
            return;
        case 10: nodes(bin); return;
        case 12: cells(bin); return;
        case 13: faces(bin); return;
        case 39:
        case 45: zone(); return;
        case 18:
            std::cout << "mesh import: Fluent periodic face pairs ignored (periodic zones are imported as ordinary patches)\n";
            break;
        case 58:
        case 59:
        case 61:
        case 62:
        case 63:
            fail("non-conformal / hanging-node Fluent meshes (section " + std::to_string(idx) + ") are not supported");
        default: break;
        }
        if (bin) endBinary();
        else skipBalanced();
    }

    void nodes(int bin) {
        const auto h = header();
        if (h.size() < 4) fail("bad node header");
        if (h[0] == 0) {  // 声明
            points.resize(std::size_t(h[2]));
            expect(')');
            return;
        }
        const glabel first = glabel(h[1]) - 1, last = glabel(h[2]) - 1;
        const int nd = h.size() >= 5 ? int(h[4]) : dim;
        if (glabel(points.size()) <= last) points.resize(std::size_t(last) + 1);
        expect('(');
        for (glabel i = first; i <= last; ++i) {
            Vec3 x;
            for (int k = 0; k < nd; ++k) x[k] = bin == 0 ? readDouble() : (bin == 2 ? double(this->bin<float>()) : this->bin<double>());
            points[i] = x;
        }
        if (bin) endBinary();
        else {
            expect(')');
            expect(')');
        }
    }

    void cells(int bin) {
        const auto h = header();
        if (h.size() < 3) fail("bad cell header");
        if (h[0] == 0) {
            nCells = std::max(nCells, glabel(h[2]));
            expect(')');
            return;
        }
        nCells = std::max(nCells, glabel(h[2]));
        if (peek('(')) {  // 混合单元类型列表：不需要
            ++p_;
            if (bin) endBinary();
            else {
                skipBalanced();
                expect(')');
            }
        } else {
            expect(')');
        }
    }

    void faces(int bin) {
        const auto h = header();
        if (h.size() < 5) {
            if (h.size() >= 3 && h[0] == 0) {
                expect(')');
                return;
            }
            fail("bad face header");
        }
        if (h[0] == 0) {
            expect(')');
            return;
        }
        Zone z{int(h[0]), int(h[3]), {}};
        const long long n = h[2] - h[1] + 1;
        const int faceType = int(h[4]);
        z.faces.reserve(std::size_t(n));
        expect('(');
        auto rd = [&]() -> long long {
            if (bin == 0) return readInt(16);
            return bin == 3 ? this->bin<std::int64_t>() : this->bin<std::int32_t>();
        };
        for (long long i = 0; i < n; ++i) {
            long long nn = faceType;
            if (faceType == 0 || faceType == 5) {
                nn = rd();
                if (faceType == 0 && nn == 5) nn = rd();
            }
            FaceRecord F;
            for (long long k = 0; k < nn; ++k) F.nodes.push_back(glabel(rd()) - 1);
            F.c0 = glabel(rd()) - 1;
            F.c1 = glabel(rd()) - 1;
            z.faces.push_back(std::move(F));
        }
        if (bin) endBinary();
        else {
            expect(')');
            expect(')');
        }
        faceZones.push_back(std::move(z));
    }

    void zone() {
        expect('(');
        const int id = int(readInt(10));
        const std::string type = readWord();
        const std::string name = readWord();
        skipBalanced();  // 头部剩余
        zoneInfo[id] = {type, name};
        skipWs();
        if (peek('(')) {
            ++p_;
            skipBalanced();
        }
        expect(')');
    }
};

PatchType fluentPatchType(const std::string& type, int bc) {
    const std::string t = lower(type);
    if (t == "wall" || (t.empty() && bc == 3)) return PatchType::Wall;
    if (t == "symmetry" || t == "axis" || (t.empty() && (bc == 7 || bc == 37))) return PatchType::Symmetry;
    return PatchType::Patch;
}

} // namespace

RawMesh readFluent(const std::string& file, const MeshImportOptions& opt) {
    FluentReader fr(file);
    fr.parse();
    if (fr.points.empty() || fr.faceZones.empty()) throw std::runtime_error(file + ": no nodes or faces found");
    for (auto& p : fr.points) p *= opt.scale;
    glabel nc = fr.nCells;
    for (auto& z : fr.faceZones)
        for (auto& F : z.faces) nc = std::max({nc, F.c0 + 1, F.c1 + 1});

    // 边界 zone → patch
    std::vector<RawPatch> patches;
    std::map<int, int> zoneToPatch;
    glabel nBaffle = 0;
    for (auto& z : fr.faceZones) {
        bool hasBoundary = false;
        for (auto& F : z.faces)
            if (F.c0 < 0 || F.c1 < 0) hasBoundary = true;
            else if (z.bcType != 2) ++nBaffle;
        if (!hasBoundary) continue;
        auto it = fr.zoneInfo.find(z.id);
        RawPatch p;
        p.name = it != fr.zoneInfo.end() ? it->second.second : "zone" + std::to_string(z.id);
        p.type = defaultPatchType(p.name, fluentPatchType(it != fr.zoneInfo.end() ? it->second.first : "", z.bcType), opt);
        zoneToPatch[z.id] = int(patches.size());
        patches.push_back(p);
    }
    if (nBaffle) std::cout << "mesh import: " << nBaffle << " two-sided faces in boundary zones treated as internal faces\n";

    if (fr.dim == 2) {
        Mesh2D m2;
        m2.points = fr.points;
        std::vector<std::array<glabel, 3>> cellEdges;
        std::map<std::pair<glabel, glabel>, int> edgePatch;
        for (auto& z : fr.faceZones)
            for (auto& F : z.faces) {
                if (F.nodes.size() != 2) throw std::runtime_error(file + ": 2D face with " + std::to_string(F.nodes.size()) + " nodes");
                const glabel a = F.nodes[0], b = F.nodes[1];
                if (F.c0 >= 0) cellEdges.push_back({F.c0, a, b});
                if (F.c1 >= 0) cellEdges.push_back({F.c1, a, b});
                if (F.c0 < 0 || F.c1 < 0) edgePatch[{std::min(a, b), std::max(a, b)}] = zoneToPatch.at(z.id);
            }
        m2.cells = polygonsFromEdges(nc, cellEdges);
        return extrudeImported(std::move(m2), patches, edgePatch, opt);
    }

    std::vector<FaceRecord> faces;
    for (auto& z : fr.faceZones)
        for (auto& F : z.faces) {
            if (F.c0 < 0 || F.c1 < 0) F.patch = zoneToPatch.at(z.id);
            faces.push_back(std::move(F));
        }
    return assembleFaceMesh(std::move(fr.points), nc, std::move(faces), patches);
}

RawMesh readMeshFile(const std::string& path, const MeshImportOptions& opt) {
    namespace fs = std::filesystem;
    if (fs::is_directory(path)) return readPolyMesh(path);
    std::ifstream in(path);
    if (!in) throw std::runtime_error("cannot open " + path);
    std::string first;
    std::getline(in, first);
    if (first.rfind("$MeshFormat", 0) == 0) return readGmsh(path, opt);
    const std::string ext = lower(fs::path(path).extension().string());
    if (ext == ".msh" || ext == ".cas" || first.find('(') != std::string::npos) return readFluent(path, opt);
    throw std::runtime_error(path + ": unknown mesh format");
}

} // namespace cfd
