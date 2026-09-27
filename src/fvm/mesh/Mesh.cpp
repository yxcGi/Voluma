#include "fvm/mesh/Mesh.h"

#include <algorithm>
#include <cmath>
#include <iostream>
#include <limits>
#include <map>
#include <set>
#include <array>
#include <stdexcept>
#include <unordered_map>
#include <cstring>

namespace cfd {

namespace {

constexpr std::uint8_t VTK_TETRA = 10, VTK_HEXAHEDRON = 12, VTK_WEDGE = 13, VTK_PYRAMID = 14, VTK_POLYHEDRON = 42;

// ------------------------------------------------ 几何（与 OpenFOAM 相同的三角形/棱锥分解）
void faceGeometry(const RawMesh& m, glabel f, Vec3& Sf, Vec3& Cf) {
    const glabel b = m.faceOffsets[f], e = m.faceOffsets[f + 1];
    const glabel n = e - b;
    const auto& P = m.points;
    const auto* fp = m.facePoints.data() + b;
    if (n == 3) {
        Cf = (P[fp[0]] + P[fp[1]] + P[fp[2]]) / 3.0;
        Sf = 0.5 * cross(P[fp[1]] - P[fp[0]], P[fp[2]] - P[fp[0]]);
        return;
    }
    Vec3 pAvg;
    for (glabel k = 0; k < n; ++k) pAvg += P[fp[k]];
    pAvg /= scalar(n);
    Vec3 sumN, sumAc;
    scalar sumA = 0;
    for (glabel k = 0; k < n; ++k) {
        const Vec3& p0 = P[fp[k]];
        const Vec3& p1 = P[fp[(k + 1) % n]];
        const Vec3 c = p0 + p1 + pAvg;
        const Vec3 nk = cross(p1 - p0, pAvg - p0);
        const scalar a = mag(nk);
        sumN += nk;
        sumA += a;
        sumAc += a * c;
    }
    Cf = sumA < VSMALL ? pAvg : sumAc / (3.0 * sumA);
    Sf = 0.5 * sumN;
}

// ------------------------------------------------ 递归坐标二分（确定性）
void bisect(const std::vector<Vec3>& C, std::vector<glabel>::iterator b, std::vector<glabel>::iterator e, int first,
            int nParts, std::vector<int>& part) {
    if (nParts == 1) {
        for (auto it = b; it != e; ++it) part[*it] = first;
        return;
    }
    Vec3 lo{GREAT, GREAT, GREAT}, hi{-GREAT, -GREAT, -GREAT};
    for (auto it = b; it != e; ++it)
        for (int a = 0; a < 3; ++a) {
            lo[a] = std::min(lo[a], C[*it][a]);
            hi[a] = std::max(hi[a], C[*it][a]);
        }
    int axis = 0;
    for (int a = 1; a < 3; ++a)
        if (hi[a] - lo[a] > hi[axis] - lo[axis]) axis = a;
    const int nLeft = nParts / 2;
    const glabel n = glabel(e - b);
    const glabel nl = (n * nLeft + nParts / 2) / nParts;
    auto mid = b + nl;
    std::nth_element(b, mid, e, [&](glabel x, glabel y) {
        if (C[x][axis] != C[y][axis]) return C[x][axis] < C[y][axis];
        return x < y;
    });
    bisect(C, b, mid, first, nLeft, part);
    bisect(C, mid, e, first + nLeft, nParts - nLeft, part);
}

std::vector<int> decomposeRCB(const std::vector<Vec3>& C, int nParts) {
    if (glabel(nParts) > glabel(C.size())) throw std::runtime_error("more processes than cells");
    std::vector<int> part(C.size(), 0);
    std::vector<glabel> ids(C.size());
    for (std::size_t i = 0; i < ids.size(); ++i) ids[i] = glabel(i);
    bisect(C, ids.begin(), ids.end(), 0, nParts, part);
    return part;
}

// ------------------------------------------------ 全局计算网格（仅 0 号进程）
struct CompPatch {
    std::string name;
    PatchType type;
    glabel start, size;
};

struct GlobalComp {
    glabel nCells = 0;
    std::vector<Vec3> C;
    std::vector<scalar> V;
    std::vector<glabel> own, nei;
    std::vector<Vec3> Sf, Cf, d;
    glabel nInternal = 0;
    std::vector<CompPatch> patches;
    std::vector<std::string> emptyNames;
    bool twoD = false;
    int emptyDir = -1;
    // VTK 单元形状
    std::vector<std::uint8_t> shape;
    std::vector<glabel> shapeOff{0}, shapePts;
    std::vector<glabel> polyOff{0}, poly;
};

GlobalComp buildGlobal(const RawMesh& m) {
    GlobalComp g;
    const glabel nF = m.nFaces(), nI = m.nInternalFaces(), nC = m.nCells;
    g.nCells = nC;
    std::vector<Vec3> Sf(nF), Cf(nF);
    for (glabel f = 0; f < nF; ++f) faceGeometry(m, f, Sf[f], Cf[f]);

    // 单元-面 CSR
    std::vector<glabel> cOff(nC + 1, 0), cList;
    for (glabel f = 0; f < nF; ++f) {
        ++cOff[m.owner[f] + 1];
        if (f < nI) ++cOff[m.neighbour[f] + 1];
    }
    for (glabel c = 0; c < nC; ++c) cOff[c + 1] += cOff[c];
    cList.resize(cOff[nC]);
    {
        std::vector<glabel> fill(cOff.begin(), cOff.end() - 1);
        for (glabel f = 0; f < nF; ++f) {
            cList[fill[m.owner[f]]++] = f;
            if (f < nI) cList[fill[m.neighbour[f]]++] = f;
        }
    }

    // 单元几何：棱锥分解
    g.C.assign(nC, Vec3{});
    g.V.assign(nC, 0.0);
    for (glabel c = 0; c < nC; ++c) {
        Vec3 cEst;
        const glabel nfc = cOff[c + 1] - cOff[c];
        for (glabel k = cOff[c]; k < cOff[c + 1]; ++k) cEst += Cf[cList[k]];
        cEst /= scalar(nfc);
        Vec3 cSum;
        scalar vSum = 0;
        for (glabel k = cOff[c]; k < cOff[c + 1]; ++k) {
            const glabel f = cList[k];
            const Vec3 S = (m.owner[f] == c) ? Sf[f] : -Sf[f];
            const scalar pv = std::max(dot(S, Cf[f] - cEst), VSMALL);
            cSum += pv * (0.75 * Cf[f] + 0.25 * cEst);
            vSum += pv;
        }
        g.C[c] = cSum / vSum;
        g.V[c] = vSum / 3.0;
    }

    // VTK 形状
    g.shape.resize(nC);
    for (glabel c = 0; c < nC; ++c) {
        std::vector<std::vector<glabel>> fs;  // 外法向顺序
        for (glabel k = cOff[c]; k < cOff[c + 1]; ++k) {
            const glabel f = cList[k];
            std::vector<glabel> pts(m.facePoints.begin() + m.faceOffsets[f], m.facePoints.begin() + m.faceOffsets[f + 1]);
            if (m.owner[f] != c) std::reverse(pts.begin(), pts.end());
            fs.push_back(std::move(pts));
        }
        std::vector<glabel> uniq;
        int nTri = 0, nQuad = 0;
        for (auto& p : fs) {
            uniq.insert(uniq.end(), p.begin(), p.end());
            nTri += p.size() == 3;
            nQuad += p.size() == 4;
        }
        std::sort(uniq.begin(), uniq.end());
        uniq.erase(std::unique(uniq.begin(), uniq.end()), uniq.end());
        const int nfc = int(fs.size()), np = int(uniq.size());
        auto inSet = [](const std::vector<glabel>& v, glabel p) { return std::find(v.begin(), v.end(), p) != v.end(); };
        // 从 base 面上的点 b 沿不在 base 上的棱找到对面的点
        auto across = [&](const std::vector<glabel>& base, glabel b) -> glabel {
            for (auto& f : fs) {
                const int n = int(f.size());
                for (int i = 0; i < n; ++i) {
                    const glabel a = f[i], a2 = f[(i + 1) % n];
                    if (a == b && !inSet(base, a2)) return a2;
                    if (a2 == b && !inSet(base, a)) return a;
                }
            }
            return -1;
        };
        std::vector<glabel> pts;
        std::uint8_t type = VTK_POLYHEDRON;
        if (nfc == 6 && nQuad == 6 && np == 8) {
            std::vector<glabel> base(fs[0].rbegin(), fs[0].rend());
            pts = base;
            for (auto b : base) pts.push_back(across(base, b));
            type = VTK_HEXAHEDRON;
        } else if (nfc == 4 && nTri == 4 && np == 4) {
            std::vector<glabel> base(fs[0].rbegin(), fs[0].rend());
            pts = base;
            for (auto p : uniq)
                if (!inSet(base, p)) pts.push_back(p);
            type = VTK_TETRA;
        } else if (nfc == 5 && nTri == 2 && nQuad == 3 && np == 6) {
            std::vector<glabel> base;
            for (auto& f : fs)
                if (f.size() == 3) {
                    base = f;
                    break;
                }
            pts = base;
            for (auto b : base) pts.push_back(across(base, b));
            type = VTK_WEDGE;
        } else if (nfc == 5 && nTri == 4 && nQuad == 1 && np == 5) {
            std::vector<glabel> base;
            for (auto& f : fs)
                if (f.size() == 4) base.assign(f.rbegin(), f.rend());
            pts = base;
            for (auto p : uniq)
                if (!inSet(base, p)) pts.push_back(p);
            type = VTK_PYRAMID;
        }
        if (type != VTK_POLYHEDRON && std::find(pts.begin(), pts.end(), glabel(-1)) != pts.end())
            type = VTK_POLYHEDRON;
        if (type == VTK_POLYHEDRON) {
            pts = uniq;
            g.poly.push_back(nfc);
            for (auto& f : fs) {
                g.poly.push_back(glabel(f.size()));
                g.poly.insert(g.poly.end(), f.begin(), f.end());
            }
        }
        g.polyOff.push_back(glabel(g.poly.size()));
        g.shape[c] = type;
        g.shapePts.insert(g.shapePts.end(), pts.begin(), pts.end());
        g.shapeOff.push_back(glabel(g.shapePts.size()));
    }

    // 计算面：内部面
    for (glabel f = 0; f < nI; ++f) {
        g.own.push_back(m.owner[f]);
        g.nei.push_back(m.neighbour[f]);
        g.Sf.push_back(Sf[f]);
        g.Cf.push_back(Cf[f]);
        g.d.push_back(g.C[m.neighbour[f]] - g.C[m.owner[f]]);
    }
    // 周期面并入内部面
    std::vector<bool> done(m.patches.size(), false);
    for (std::size_t ia = 0; ia < m.patches.size(); ++ia) {
        const auto& A = m.patches[ia];
        if (A.type != PatchType::Cyclic || done[ia]) continue;
        std::size_t ib = m.patches.size();
        for (std::size_t k = 0; k < m.patches.size(); ++k)
            if (m.patches[k].name == A.neighbourPatch) ib = k;
        if (ib == m.patches.size()) throw std::runtime_error("cyclic patch " + A.name + ": neighbourPatch not found");
        const auto& B = m.patches[ib];
        if (A.size != B.size) throw std::runtime_error("cyclic patches " + A.name + "/" + B.name + " differ in size");
        done[ia] = done[ib] = true;
        Vec3 ca, cb;
        scalar minA = GREAT;
        for (glabel k = 0; k < A.size; ++k) {
            ca += Cf[A.start + k];
            cb += Cf[B.start + k];
            minA = std::min({minA, mag(Sf[A.start + k]), mag(Sf[B.start + k])});
        }
        const Vec3 sep = (cb - ca) / scalar(A.size);
        const scalar h = 0.25 * std::sqrt(minA);
        auto key = [&](const Vec3& p) {
            return std::array<glabel, 3>{glabel(std::floor(p.x / h)), glabel(std::floor(p.y / h)), glabel(std::floor(p.z / h))};
        };
        std::map<std::array<glabel, 3>, std::vector<glabel>> grid;
        for (glabel k = 0; k < B.size; ++k) grid[key(Cf[B.start + k])].push_back(B.start + k);
        for (glabel k = 0; k < A.size; ++k) {
            const glabel fa = A.start + k;
            const Vec3 target = Cf[fa] + sep;
            const auto kk = key(target);
            glabel best = -1;
            scalar bestD = GREAT;
            for (int dx = -1; dx <= 1; ++dx)
                for (int dy = -1; dy <= 1; ++dy)
                    for (int dz = -1; dz <= 1; ++dz) {
                        auto it = grid.find({kk[0] + dx, kk[1] + dy, kk[2] + dz});
                        if (it == grid.end()) continue;
                        for (glabel fb : it->second) {
                            const scalar dd = mag(Cf[fb] - target);
                            if (dd < bestD) {
                                bestD = dd;
                                best = fb;
                            }
                        }
                    }
            if (best < 0 || bestD > 1e-3 * std::sqrt(mag(Sf[fa])))
                throw std::runtime_error("cyclic patch " + A.name + ": no matching face (non-translational?)");
            const glabel o = m.owner[fa], n = m.owner[best];
            g.own.push_back(o);
            g.nei.push_back(n);
            g.Sf.push_back(Sf[fa]);
            g.Cf.push_back(Cf[fa]);
            g.d.push_back(g.C[n] + (Cf[fa] - Cf[best]) - g.C[o]);
        }
    }
    g.nInternal = glabel(g.own.size());
    // 边界面
    Vec3 emptyN;
    for (const auto& p : m.patches) {
        if (p.type == PatchType::Cyclic) continue;
        if (p.type == PatchType::Empty) {
            g.emptyNames.push_back(p.name);
            g.twoD = true;
            for (glabel k = 0; k < p.size; ++k) {
                const Vec3 s = Sf[p.start + k];
                emptyN += Vec3{std::abs(s.x), std::abs(s.y), std::abs(s.z)};
            }
            continue;
        }
        CompPatch cp{p.name, p.type, glabel(g.own.size()), p.size};
        for (glabel k = 0; k < p.size; ++k) {
            const glabel f = p.start + k;
            g.own.push_back(m.owner[f]);
            g.Sf.push_back(Sf[f]);
            g.Cf.push_back(Cf[f]);
            g.d.push_back(Cf[f] - g.C[m.owner[f]]);
        }
        g.patches.push_back(cp);
    }
    if (g.twoD) {
        g.emptyDir = 0;
        for (int a = 1; a < 3; ++a)
            if (emptyN[a] > emptyN[g.emptyDir]) g.emptyDir = a;
    }
    return g;
}

std::vector<char> packRank(const GlobalComp& g, const RawMesh& m, const std::vector<int>& part, int r,
                           const std::vector<glabel>& owned, const std::vector<glabel>& faces,
                           const std::vector<glabel>& ghosts, const std::vector<int>& nbrs,
                           const std::vector<std::vector<glabel>>& sendG, const std::vector<std::vector<glabel>>& recvG) {
    par::Packer pk;
    std::unordered_map<glabel, label> loc;
    loc.reserve(owned.size() + ghosts.size());
    for (std::size_t i = 0; i < owned.size(); ++i) loc[owned[i]] = label(i);
    for (std::size_t i = 0; i < ghosts.size(); ++i) loc[ghosts[i]] = label(owned.size() + i);

    pk.put<glabel>(g.nCells);
    pk.put<label>(label(owned.size()));
    pk.put<label>(label(ghosts.size()));
    std::vector<glabel> cg(owned);
    cg.insert(cg.end(), ghosts.begin(), ghosts.end());
    std::vector<Vec3> C;
    std::vector<scalar> V;
    for (auto c : cg) {
        C.push_back(g.C[c]);
        V.push_back(g.V[c]);
    }
    pk.putVec(cg);
    pk.putVec(C);
    pk.putVec(V);

    std::vector<label> own, nei;
    std::vector<Vec3> Sf, Cf, d;
    label nInt = 0;
    for (auto f : faces) {
        own.push_back(loc.at(g.own[f]));
        if (f < g.nInternal) {
            nei.push_back(loc.at(g.nei[f]));
            ++nInt;
        }
        Sf.push_back(g.Sf[f]);
        Cf.push_back(g.Cf[f]);
        d.push_back(g.d[f]);
    }
    pk.put<label>(nInt);
    pk.putVec(faces);
    pk.put<glabel>(glabel(g.own.size()));
    pk.putVec(own);
    pk.putVec(nei);
    pk.putVec(Sf);
    pk.putVec(Cf);
    pk.putVec(d);

    // patches
    pk.put<std::uint64_t>(g.patches.size());
    for (const auto& p : g.patches) {
        label start = -1, size = 0;
        for (std::size_t k = 0; k < faces.size(); ++k)
            if (faces[k] >= p.start && faces[k] < p.start + p.size) {
                if (start < 0) start = label(k);
                ++size;
            }
        if (start < 0) start = label(faces.size());
        pk.putString(p.name);
        pk.put<int>(int(p.type));
        pk.put<label>(start);
        pk.put<label>(size);
        pk.put<glabel>(p.size);
    }
    pk.put<std::uint64_t>(g.emptyNames.size());
    for (auto& s : g.emptyNames) pk.putString(s);
    pk.put<int>(g.twoD ? 1 : 0);
    pk.put<int>(g.emptyDir);

    // halo
    pk.putVec(nbrs);
    for (std::size_t n = 0; n < nbrs.size(); ++n) {
        std::vector<label> s, rr;
        for (auto c : sendG[n]) s.push_back(loc.at(c));
        for (auto c : recvG[n]) rr.push_back(loc.at(c));
        pk.putVec(s);
        pk.putVec(rr);
    }

    // VTK piece
    std::vector<glabel> gp;
    for (auto c : owned) gp.insert(gp.end(), g.shapePts.begin() + g.shapeOff[c], g.shapePts.begin() + g.shapeOff[c + 1]);
    std::sort(gp.begin(), gp.end());
    gp.erase(std::unique(gp.begin(), gp.end()), gp.end());
    std::unordered_map<glabel, std::int64_t> ploc;
    for (std::size_t i = 0; i < gp.size(); ++i) ploc[gp[i]] = std::int64_t(i);
    std::vector<Vec3> pts;
    for (auto p : gp) pts.push_back(m.points[p]);
    std::vector<std::int64_t> conn, off, fcs, foff;
    std::vector<std::uint8_t> types;
    bool hasPoly = false;
    for (auto c : owned) {
        for (glabel k = g.shapeOff[c]; k < g.shapeOff[c + 1]; ++k) conn.push_back(ploc.at(g.shapePts[k]));
        off.push_back(std::int64_t(conn.size()));
        types.push_back(g.shape[c]);
        if (g.shape[c] == VTK_POLYHEDRON) {
            hasPoly = true;
            glabel k = g.polyOff[c];
            const glabel nfc = g.poly[k++];
            fcs.push_back(nfc);
            for (glabel f = 0; f < nfc; ++f) {
                const glabel n = g.poly[k++];
                fcs.push_back(n);
                for (glabel i = 0; i < n; ++i) fcs.push_back(ploc.at(g.poly[k++]));
            }
            foff.push_back(std::int64_t(fcs.size()));
        } else {
            foff.push_back(-1);
        }
    }
    pk.putVec(pts);
    pk.putVec(conn);
    pk.putVec(off);
    pk.putVec(types);
    pk.put<int>(hasPoly ? 1 : 0);
    pk.putVec(fcs);
    pk.putVec(foff);
    (void)part;
    (void)r;
    return std::move(pk.buf);
}

} // namespace

// =================================================================== build
std::shared_ptr<Mesh> Mesh::build(const RawMesh* raw) {
    const int np = par::size();
    auto mesh = std::make_shared<Mesh>();
    std::vector<char> myPack;
    // 0 号进程上：各进程条目在全局序列中的位置
    std::vector<glabel> cellPos;
    std::vector<std::vector<glabel>> patchPos;

    if (par::master()) {
        if (!raw) throw std::runtime_error("Mesh::build: raw mesh required on master");
        const GlobalComp g = buildGlobal(*raw);
        const std::vector<int> part = decomposeRCB(g.C, np);
        std::vector<std::vector<glabel>> owned(np), faces(np), ghosts(np);
        for (glabel c = 0; c < g.nCells; ++c) owned[part[c]].push_back(c);
        for (glabel f = 0; f < glabel(g.own.size()); ++f) {
            const int po = part[g.own[f]];
            faces[po].push_back(f);
            if (f < g.nInternal) {
                const int pn = part[g.nei[f]];
                if (pn != po) {
                    faces[pn].push_back(f);
                    ghosts[po].push_back(g.nei[f]);
                    ghosts[pn].push_back(g.own[f]);
                }
            }
        }
        for (int r = 0; r < np; ++r) {
            std::sort(faces[r].begin(), faces[r].end());
            std::sort(ghosts[r].begin(), ghosts[r].end());
            ghosts[r].erase(std::unique(ghosts[r].begin(), ghosts[r].end()), ghosts[r].end());
        }
        // recv[r][q]：r 从 q 接收的幽灵（全局号，升序）；send[q][r] 与之相同
        std::vector<std::map<int, std::vector<glabel>>> recv(np), send(np);
        for (int r = 0; r < np; ++r)
            for (auto c : ghosts[r]) recv[r][part[c]].push_back(c);
        for (int r = 0; r < np; ++r)
            for (auto& [q, list] : recv[r]) send[q][r] = list;
        for (int r = 0; r < np; ++r) {
            std::set<int> nb;
            for (auto& kv : recv[r]) nb.insert(kv.first);
            for (auto& kv : send[r]) nb.insert(kv.first);
            std::vector<int> nbrs(nb.begin(), nb.end());
            std::vector<std::vector<glabel>> sG, rG;
            for (int q : nbrs) {
                sG.push_back(send[r].count(q) ? send[r][q] : std::vector<glabel>{});
                rG.push_back(recv[r].count(q) ? recv[r][q] : std::vector<glabel>{});
            }
            auto pack = packRank(g, *raw, part, r, owned[r], faces[r], ghosts[r], nbrs, sG, rG);
            if (r == 0)
                myPack = std::move(pack);
            else
                par::sendBytes(pack, r, 100 + r);
        }
        for (int r = 0; r < np; ++r) cellPos.insert(cellPos.end(), owned[r].begin(), owned[r].end());
        mesh->faceListsOnMaster_ = faces;
        patchPos.resize(g.patches.size());
        for (std::size_t p = 0; p < g.patches.size(); ++p)
            for (int r = 0; r < np; ++r)
                for (auto f : faces[r])
                    if (f >= g.patches[p].start && f < g.patches[p].start + g.patches[p].size)
                        patchPos[p].push_back(f - g.patches[p].start);
    } else {
        myPack = par::recvBytes(0, 100 + par::rank());
    }

    // ------------------------------------------------ 解包
    Mesh& M = *mesh;
    par::Unpacker u(myPack);
    M.nGlobalCells_ = u.get<glabel>();
    M.nCells_ = u.get<label>();
    M.nGhost_ = u.get<label>();
    M.cellGlobal_ = u.getVec<glabel>();
    M.C_ = u.getVec<Vec3>();
    M.V_ = u.getVec<scalar>();
    M.nInternalFaces_ = u.get<label>();
    M.faceGlobal_ = u.getVec<glabel>();
    M.nGlobalFaces_ = u.get<glabel>();
    M.owner_ = u.getVec<label>();
    M.neighbour_ = u.getVec<label>();
    M.Sf_ = u.getVec<Vec3>();
    M.Cf_ = u.getVec<Vec3>();
    M.d_ = u.getVec<Vec3>();
    const auto nP = u.get<std::uint64_t>();
    for (std::uint64_t p = 0; p < nP; ++p) {
        Patch pt;
        pt.name = u.getString();
        pt.type = PatchType(u.get<int>());
        pt.start = u.get<label>();
        pt.size = u.get<label>();
        pt.globalSize = u.get<glabel>();
        M.patches_.push_back(std::move(pt));
    }
    const auto nE = u.get<std::uint64_t>();
    for (std::uint64_t k = 0; k < nE; ++k) M.emptyPatches_.push_back(u.getString());
    M.twoD_ = u.get<int>() != 0;
    M.emptyDir_ = u.get<int>();
    {
        auto nbrs = u.getVec<int>();
        std::vector<std::vector<label>> s(nbrs.size()), r(nbrs.size());
        for (std::size_t n = 0; n < nbrs.size(); ++n) {
            s[n] = u.getVec<label>();
            r[n] = u.getVec<label>();
        }
        M.halo_.setup(std::move(nbrs), std::move(s), std::move(r));
    }
    M.vtk_.points = u.getVec<Vec3>();
    M.vtk_.connectivity = u.getVec<std::int64_t>();
    M.vtk_.offsets = u.getVec<std::int64_t>();
    M.vtk_.types = u.getVec<std::uint8_t>();
    M.vtk_.hasPolyhedra = u.get<int>() != 0;
    M.vtk_.faces = u.getVec<std::int64_t>();
    M.vtk_.faceOffsets = u.getVec<std::int64_t>();

    M.cellOrdering_.setup(M.nCells_, M.nGlobalCells_, std::move(cellPos));
    for (std::size_t p = 0; p < M.patches_.size(); ++p)
        M.patches_[p].ordering.setup(M.patches_[p].size, M.patches_[p].globalSize,
                                     par::master() ? std::move(patchPos[p]) : std::vector<glabel>{});
    M.finishLocal();
    return mesh;
}

void Mesh::finishLocal() {
    const label nF = nFaces();
    magSf_.resize(nF);
    w_.resize(nF);
    dc_.resize(nF);
    ndc_.resize(nF);
    corr_.resize(nF);
    for (label f = 0; f < nF; ++f) {
        magSf_[f] = mag(Sf_[f]);
        const Vec3 n = Sf_[f] / magSf_[f];
        const Vec3& dd = d_[f];
        if (f < nInternalFaces_) {
            const Vec3& CP = C_[owner_[f]];
            const scalar sOwn = std::abs(dot(Sf_[f], Cf_[f] - CP));
            const scalar sNei = std::abs(dot(Sf_[f], CP + dd - Cf_[f]));
            w_[f] = sNei / (sOwn + sNei);
        } else {
            w_[f] = 1.0;
        }
        dc_[f] = 1.0 / mag(dd);
        ndc_[f] = 1.0 / std::max(dot(n, dd), 0.05 * mag(dd));
        corr_[f] = n - dd * ndc_[f];
    }
    // 自有单元的面 CSR（按局部面顺序）
    cfOff_.assign(nCells_ + 1, 0);
    for (label f = 0; f < nF; ++f) {
        if (owner_[f] < nCells_) ++cfOff_[owner_[f] + 1];
        if (f < nInternalFaces_ && neighbour_[f] < nCells_) ++cfOff_[neighbour_[f] + 1];
    }
    for (label c = 0; c < nCells_; ++c) cfOff_[c + 1] += cfOff_[c];
    cfList_.resize(cfOff_[nCells_]);
    std::vector<label> fill(cfOff_.begin(), cfOff_.end() - 1);
    for (label f = 0; f < nF; ++f) {
        if (owner_[f] < nCells_) cfList_[fill[owner_[f]]++] = f;
        if (f < nInternalFaces_ && neighbour_[f] < nCells_) cfList_[fill[neighbour_[f]]++] = f;
    }
    solved_.clear();
    for (int c = 0; c < 3; ++c)
        if (!twoD_ || c != emptyDir_) solved_.push_back(c);
    par::SumAcc v;
    for (label c = 0; c < nCells_; ++c) v.add(V_[c]);
    totalVolume_ = v.allReduce();
}

std::vector<double> Mesh::gatherFaces(const double* local, int nc) const {
    std::vector<double> out;
    if (par::master()) out.assign(std::size_t(nGlobalFaces_) * nc, 0.0);
    if (par::master()) {
        for (std::size_t k = 0; k < faceGlobal_.size(); ++k)
            for (int c = 0; c < nc; ++c) out[faceGlobal_[k] * nc + c] = local[k * nc + c];
        for (int r = 1; r < par::size(); ++r) {
            auto buf = par::recvBytes(r, 300);
            const double* v = reinterpret_cast<const double*>(buf.data());
            const auto& fl = faceListsOnMaster_[r];
            for (std::size_t k = 0; k < fl.size(); ++k)
                for (int c = 0; c < nc; ++c) out[fl[k] * nc + c] = v[k * nc + c];
        }
    } else {
        std::vector<char> buf(sizeof(double) * faceGlobal_.size() * nc);
        std::memcpy(buf.data(), local, buf.size());
        par::sendBytes(buf, 0, 300);
    }
    return out;
}

void Mesh::scatterFaces(const std::vector<double>& global, double* local, int nc) const {
    if (par::master()) {
        for (int r = 1; r < par::size(); ++r) {
            const auto& fl = faceListsOnMaster_[r];
            std::vector<char> buf(sizeof(double) * fl.size() * nc);
            double* v = reinterpret_cast<double*>(buf.data());
            for (std::size_t k = 0; k < fl.size(); ++k)
                for (int c = 0; c < nc; ++c) v[k * nc + c] = global[fl[k] * nc + c];
            par::sendBytes(buf, r, 301);
        }
        for (std::size_t k = 0; k < faceGlobal_.size(); ++k)
            for (int c = 0; c < nc; ++c) local[k * nc + c] = global[faceGlobal_[k] * nc + c];
    } else {
        auto buf = par::recvBytes(0, 301);
        std::memcpy(local, buf.data(), buf.size());
    }
}

label Mesh::findPatch(const std::string& name) const {
    for (std::size_t p = 0; p < patches_.size(); ++p)
        if (patches_[p].name == name) return label(p);
    return -1;
}

bool Mesh::isEmptyPatch(const std::string& name) const {
    return std::find(emptyPatches_.begin(), emptyPatches_.end(), name) != emptyPatches_.end();
}

void Mesh::printSummary() const {
    const glabel nf = par::allSum(glabel(nInternalFaces_));
    const glabel ng = par::allSum(glabel(nGhost_));
    const double minC = par::allMin(nCells_), maxC = par::allMax(nCells_);
    scalar maxNonOrth = 0;
    for (label f = 0; f < nInternalFaces_; ++f) {
        const scalar c = dot(Sf_[f], d_[f]) / (magSf_[f] * mag(d_[f]));
        maxNonOrth = std::max(maxNonOrth, std::acos(std::clamp(c, -1.0, 1.0)) * 180.0 / PI);
    }
    maxNonOrth = par::allMax(maxNonOrth);
    std::cout << "Mesh: " << nGlobalCells_ << " cells, " << (twoD_ ? "2D" : "3D") << ", total volume "
              << totalVolume_ << ", max non-orthogonality " << maxNonOrth << " deg\n";
    std::cout << "  patches:";
    for (const auto& p : patches_) std::cout << ' ' << p.name << '(' << toString(p.type) << ',' << p.globalSize << ')';
    for (const auto& e : emptyPatches_) std::cout << ' ' << e << "(empty)";
    std::cout << '\n';
    if (par::parallel())
        std::cout << "  " << par::size() << " processes, cells/process " << minC << ".." << maxC
                  << ", inter-process+internal faces " << nf << ", ghost cells " << ng << '\n';
}

} // namespace cfd
