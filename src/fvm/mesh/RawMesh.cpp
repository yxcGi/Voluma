#include "fvm/mesh/RawMesh.h"

#include <algorithm>
#include <cstdlib>
#include <cstring>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <sstream>
#include <stdexcept>

namespace cfd {

PatchType patchTypeFromString(const std::string& s) {
    if (s == "wall") return PatchType::Wall;
    if (s == "symmetry" || s == "symmetryPlane") return PatchType::Symmetry;
    if (s == "empty") return PatchType::Empty;
    if (s == "cyclic") return PatchType::Cyclic;
    return PatchType::Patch;
}

std::string toString(PatchType t) {
    switch (t) {
    case PatchType::Wall: return "wall";
    case PatchType::Symmetry: return "symmetryPlane";
    case PatchType::Empty: return "empty";
    case PatchType::Cyclic: return "cyclic";
    default: return "patch";
    }
}

namespace {

// ------------------------------------------------ OpenFOAM ASCII 文件分词
class FoamTokens {
public:
    explicit FoamTokens(const std::string& path) {
        std::ifstream in(path, std::ios::binary);
        if (!in) throw std::runtime_error("cannot open " + path);
        std::stringstream ss;
        ss << in.rdbuf();
        text_ = ss.str();
        stripComments();
        skipHeader();
    }
    bool eof() {
        skipWs();
        return pos_ >= text_.size();
    }
    std::string next() {
        skipWs();
        if (pos_ >= text_.size()) throw std::runtime_error("unexpected end of file");
        const char c = text_[pos_];
        if (c == '(' || c == ')' || c == '{' || c == '}' || c == ';') {
            ++pos_;
            return std::string(1, c);
        }
        const std::size_t b = pos_;
        while (pos_ < text_.size() && !std::isspace(static_cast<unsigned char>(text_[pos_])) &&
               std::string("(){};").find(text_[pos_]) == std::string::npos)
            ++pos_;
        return text_.substr(b, pos_ - b);
    }
    void expect(const std::string& t) {
        const auto s = next();
        if (s != t) throw std::runtime_error("expected '" + t + "' but got '" + s + "'");
    }
    glabel nextInt() { return std::strtoll(next().c_str(), nullptr, 10); }
    double nextDouble() { return std::strtod(next().c_str(), nullptr); }
    // 读列表前缀：N ( ；也支持 N{v} 的均匀列表（返回 -1 表示均匀）
    glabel listStart() {
        const glabel n = nextInt();
        skipWs();
        if (pos_ < text_.size() && text_[pos_] == '{') throw std::runtime_error("uniform lists are not supported");
        expect("(");
        return n;
    }
    const std::string& headerClass() const { return class_; }

private:
    void skipWs() {
        while (pos_ < text_.size() && std::isspace(static_cast<unsigned char>(text_[pos_]))) ++pos_;
    }
    void stripComments() {
        std::string out;
        out.reserve(text_.size());
        for (std::size_t i = 0; i < text_.size();) {
            if (text_.compare(i, 2, "//") == 0) {
                while (i < text_.size() && text_[i] != '\n') ++i;
            } else if (text_.compare(i, 2, "/*") == 0) {
                const auto e = text_.find("*/", i + 2);
                i = (e == std::string::npos) ? text_.size() : e + 2;
            } else {
                out.push_back(text_[i++]);
            }
        }
        text_.swap(out);
    }
    void skipHeader() {
        const auto p = text_.find("FoamFile");
        if (p == std::string::npos) return;
        const auto b = text_.find('{', p);
        const auto e = text_.find('}', b);
        const std::string hdr = text_.substr(b + 1, e - b - 1);
        std::istringstream hs(hdr);
        std::string w;
        while (hs >> w)
            if (w == "class") {
                hs >> class_;
                if (!class_.empty() && class_.back() == ';') class_.pop_back();
            }
        pos_ = e + 1;
    }
    std::string text_;
    std::size_t pos_ = 0;
    std::string class_;
};

// ------------------------------------------------ OpenFOAM 二进制列表文件
// 头部 format binary; arch "LSB;label=32;scalar=64";  列表：N ( <N×元素字节> )
struct FoamBinaryFile {
    std::string text;
    std::size_t pos = 0;
    int labelBytes = 4, scalarBytes = 8;
    bool binary = false;
    std::string cls;

    explicit FoamBinaryFile(const std::string& path) {
        std::ifstream in(path, std::ios::binary);
        if (!in) throw std::runtime_error("cannot open " + path);
        std::stringstream ss;
        ss << in.rdbuf();
        text = ss.str();
        const auto p = text.find("FoamFile");
        if (p == std::string::npos) return;
        const auto b = text.find('{', p), e = text.find('}', b);
        const std::string hdr = text.substr(b + 1, e - b - 1);
        std::istringstream hs(hdr);
        std::string w;
        while (hs >> w) {
            if (w == "format") {
                hs >> w;
                binary = w.rfind("binary", 0) == 0;
            } else if (w == "class") {
                hs >> cls;
                if (!cls.empty() && cls.back() == ';') cls.pop_back();
            }
        }
        const auto a = hdr.find("arch");
        if (a != std::string::npos) {
            const auto l = hdr.find("label=", a), sc = hdr.find("scalar=", a);
            if (l != std::string::npos) labelBytes = std::atoi(hdr.c_str() + l + 6) / 8;
            if (sc != std::string::npos) scalarBytes = std::atoi(hdr.c_str() + sc + 7) / 8;
            if (hdr.find("MSB", a) != std::string::npos) throw std::runtime_error(path + ": big-endian binary not supported");
        }
        pos = e + 1;
    }
    void skipWsComments() {
        while (pos < text.size()) {
            if (std::isspace(static_cast<unsigned char>(text[pos]))) ++pos;
            else if (text.compare(pos, 2, "//") == 0) pos = text.find('\n', pos);
            else if (text.compare(pos, 2, "/*") == 0) pos = text.find("*/", pos) + 2;
            else break;
        }
    }
    // 返回列表元素起始字节与数目
    std::pair<const char*, glabel> list(std::size_t elemBytes) {
        skipWsComments();
        char* e;
        const glabel n = std::strtoll(text.c_str() + pos, &e, 10);
        pos = std::size_t(e - text.c_str());
        skipWsComments();
        if (text[pos] != '(') throw std::runtime_error("binary list: expected '('");
        ++pos;
        const char* data = text.data() + pos;
        pos += std::size_t(n) * elemBytes;
        if (pos >= text.size() || text[pos] != ')') throw std::runtime_error("binary list: size mismatch");
        ++pos;
        return {data, n};
    }
    std::vector<glabel> labels() {
        auto [d, n] = list(labelBytes);
        std::vector<glabel> v(n);
        for (glabel i = 0; i < n; ++i) {
            if (labelBytes == 4) {
                std::int32_t x;
                std::memcpy(&x, d + 4 * i, 4);
                v[i] = x;
            } else {
                std::int64_t x;
                std::memcpy(&x, d + 8 * i, 8);
                v[i] = glabel(x);
            }
        }
        return v;
    }
    std::vector<Vec3> vectors() {
        auto [d, n] = list(3 * scalarBytes);
        std::vector<Vec3> v(n);
        for (glabel i = 0; i < n; ++i)
            for (int k = 0; k < 3; ++k) {
                if (scalarBytes == 8) {
                    double x;
                    std::memcpy(&x, d + 24 * i + 8 * k, 8);
                    v[i][k] = x;
                } else {
                    float x;
                    std::memcpy(&x, d + 12 * i + 4 * k, 4);
                    v[i][k] = x;
                }
            }
        return v;
    }
};

void writeHeader(std::ostream& os, const std::string& cls, const std::string& obj) {
    os << "FoamFile\n{\n    format      ascii;\n    class       " << cls << ";\n    location    \"constant/polyMesh\";\n"
       << "    object      " << obj << ";\n}\n\n";
}

} // namespace

RawMesh readPolyMesh(const std::string& dir) {
    RawMesh m;
    if (!std::filesystem::exists(dir + "/points") && std::filesystem::exists(dir + "/points.gz"))
        throw std::runtime_error(dir + ": compressed polyMesh (*.gz) is not supported; gunzip the files first");
    if (FoamBinaryFile(dir + "/points").binary) {
        m.points = FoamBinaryFile(dir + "/points").vectors();
        FoamBinaryFile ff(dir + "/faces");
        if (ff.cls != "faceCompactList") throw std::runtime_error(dir + "/faces: binary faces must be faceCompactList");
        m.faceOffsets = ff.labels();
        m.facePoints = ff.labels();
        m.owner = FoamBinaryFile(dir + "/owner").labels();
        m.neighbour = FoamBinaryFile(dir + "/neighbour").labels();
    } else {
    {
        FoamTokens t(dir + "/points");
        const glabel n = t.listStart();
        m.points.resize(n);
        for (glabel i = 0; i < n; ++i) {
            t.expect("(");
            m.points[i].x = t.nextDouble();
            m.points[i].y = t.nextDouble();
            m.points[i].z = t.nextDouble();
            t.expect(")");
        }
    }
    {
        FoamTokens t(dir + "/faces");
        if (t.headerClass() == "faceCompactList") {
            const glabel n1 = t.listStart();
            std::vector<glabel> off(n1);
            for (auto& o : off) o = t.nextInt();
            t.expect(")");
            const glabel n2 = t.listStart();
            m.facePoints.resize(n2);
            for (auto& p : m.facePoints) p = t.nextInt();
            m.faceOffsets = off;
        } else {
            const glabel n = t.listStart();
            m.faceOffsets.reserve(n + 1);
            for (glabel f = 0; f < n; ++f) {
                // 形如 4(0 1 2 3)
                std::string tok = t.next();
                const glabel np = std::strtoll(tok.c_str(), nullptr, 10);
                t.expect("(");
                for (glabel k = 0; k < np; ++k) m.facePoints.push_back(t.nextInt());
                t.expect(")");
                m.faceOffsets.push_back(glabel(m.facePoints.size()));
            }
        }
    }
    {
        FoamTokens t(dir + "/owner");
        const glabel n = t.listStart();
        m.owner.resize(n);
        for (auto& o : m.owner) o = t.nextInt();
    }
    {
        FoamTokens t(dir + "/neighbour");
        const glabel n = t.listStart();
        m.neighbour.resize(n);
        for (auto& o : m.neighbour) o = t.nextInt();
    }
    }
    {
        FoamTokens t(dir + "/boundary");
        const glabel n = t.listStart();
        for (glabel p = 0; p < n; ++p) {
            RawPatch rp;
            rp.name = t.next();
            t.expect("{");
            int depth = 1;
            while (depth > 0) {
                std::string key = t.next();
                if (key == "}") {
                    --depth;
                    continue;
                }
                if (key == "{") {
                    ++depth;
                    continue;
                }
                if (key == ";") continue;
                // 读取值直到 ';'
                std::vector<std::string> vals;
                for (std::string v = t.next(); v != ";"; v = t.next()) {
                    if (v == "{") {  // 嵌套字典：跳过
                        int d = 1;
                        while (d > 0) {
                            const auto w = t.next();
                            if (w == "{") ++d;
                            if (w == "}") --d;
                        }
                        break;
                    }
                    vals.push_back(v);
                }
                if (vals.empty()) continue;
                if (key == "type") rp.type = patchTypeFromString(vals[0]);
                else if (key == "nFaces") rp.size = std::strtoll(vals[0].c_str(), nullptr, 10);
                else if (key == "startFace") rp.start = std::strtoll(vals[0].c_str(), nullptr, 10);
                else if (key == "neighbourPatch") rp.neighbourPatch = vals[0];
            }
            m.patches.push_back(rp);
        }
    }
    glabel nc = 0;
    for (auto o : m.owner) nc = std::max(nc, o + 1);
    for (auto o : m.neighbour) nc = std::max(nc, o + 1);
    m.nCells = nc;
    if (m.faceOffsets.size() != m.owner.size() + 1)
        throw std::runtime_error("polyMesh: faces and owner sizes differ in " + dir);
    return m;
}

void writePolyMesh(const RawMesh& m, const std::string& dir) {
    std::filesystem::create_directories(dir);
    {
        std::ofstream os(dir + "/points");
        os << std::setprecision(17);
        writeHeader(os, "vectorField", "points");
        os << m.points.size() << "\n(\n";
        for (const auto& p : m.points) os << '(' << p.x << ' ' << p.y << ' ' << p.z << ")\n";
        os << ")\n";
    }
    {
        std::ofstream os(dir + "/faces");
        writeHeader(os, "faceList", "faces");
        os << m.nFaces() << "\n(\n";
        for (glabel f = 0; f < m.nFaces(); ++f) {
            os << (m.faceOffsets[f + 1] - m.faceOffsets[f]) << '(';
            for (glabel k = m.faceOffsets[f]; k < m.faceOffsets[f + 1]; ++k)
                os << m.facePoints[k] << (k + 1 < m.faceOffsets[f + 1] ? " " : "");
            os << ")\n";
        }
        os << ")\n";
    }
    {
        std::ofstream os(dir + "/owner");
        writeHeader(os, "labelList", "owner");
        os << m.owner.size() << "\n(\n";
        for (auto o : m.owner) os << o << '\n';
        os << ")\n";
    }
    {
        std::ofstream os(dir + "/neighbour");
        writeHeader(os, "labelList", "neighbour");
        os << m.neighbour.size() << "\n(\n";
        for (auto o : m.neighbour) os << o << '\n';
        os << ")\n";
    }
    {
        std::ofstream os(dir + "/boundary");
        writeHeader(os, "polyBoundaryMesh", "boundary");
        os << m.patches.size() << "\n(\n";
        for (const auto& p : m.patches) {
            os << "    " << p.name << "\n    {\n        type            " << toString(p.type) << ";\n";
            if (p.type == PatchType::Cyclic) os << "        neighbourPatch  " << p.neighbourPatch << ";\n";
            os << "        nFaces          " << p.size << ";\n        startFace       " << p.start << ";\n    }\n";
        }
        os << ")\n";
    }
}

// ------------------------------------------------------------ 结构网格生成
std::function<scalar(scalar)> tanhStretch(scalar beta) {
    return [beta](scalar t) { return 0.5 * (1.0 + std::tanh(beta * (2.0 * t - 1.0)) / std::tanh(beta)); };
}

RawMesh generateBox(const BoxSpec& s) {
    const int nx = s.n[0], ny = s.n[1], nz = s.n[2];
    if (nx < 1 || ny < 1 || nz < 1) throw std::invalid_argument("generateBox: bad resolution");
    if (s.twoD && nz != 1) throw std::invalid_argument("generateBox: twoD requires n[2]==1");
    RawMesh m;
    auto P = [&](int i, int j, int k) -> glabel { return i + glabel(nx + 1) * (j + glabel(ny + 1) * k); };
    auto C = [&](int i, int j, int k) -> glabel { return i + glabel(nx) * (j + glabel(ny) * k); };

    auto coord = [&](int d, int i, int n) {
        scalar t = scalar(i) / n;
        if (s.stretch[d]) t = s.stretch[d](t);
        return s.lo[d] + (s.hi[d] - s.lo[d]) * t;
    };
    m.points.resize(std::size_t(nx + 1) * (ny + 1) * (nz + 1));
    for (int k = 0; k <= nz; ++k)
        for (int j = 0; j <= ny; ++j)
            for (int i = 0; i <= nx; ++i) {
                Vec3 p{coord(0, i, nx), coord(1, j, ny), coord(2, k, nz)};
                if (s.map) p = s.map(p);
                m.points[P(i, j, k)] = p;
            }
    m.nCells = glabel(nx) * ny * nz;

    // 内部面：按 owner 单元顺序，依次 +x +y +z
    for (int k = 0; k < nz; ++k)
        for (int j = 0; j < ny; ++j)
            for (int i = 0; i < nx; ++i) {
                const glabel c = C(i, j, k);
                if (i + 1 < nx) {
                    m.addFace({P(i + 1, j, k), P(i + 1, j + 1, k), P(i + 1, j + 1, k + 1), P(i + 1, j, k + 1)});
                    m.owner.push_back(c);
                    m.neighbour.push_back(C(i + 1, j, k));
                }
                if (j + 1 < ny) {
                    m.addFace({P(i, j + 1, k), P(i, j + 1, k + 1), P(i + 1, j + 1, k + 1), P(i + 1, j + 1, k)});
                    m.owner.push_back(c);
                    m.neighbour.push_back(C(i, j + 1, k));
                }
                if (k + 1 < nz) {
                    m.addFace({P(i, j, k + 1), P(i + 1, j, k + 1), P(i + 1, j + 1, k + 1), P(i, j + 1, k + 1)});
                    m.owner.push_back(c);
                    m.neighbour.push_back(C(i, j, k + 1));
                }
            }

    auto beginPatch = [&](const std::string& name, PatchType t) {
        RawPatch rp;
        rp.name = name;
        rp.type = t;
        rp.start = m.nFaces();
        m.patches.push_back(rp);
    };
    auto endPatch = [&]() { m.patches.back().size = m.nFaces() - m.patches.back().start; };

    // xMin / xMax
    for (int side = 0; side < 2; ++side) {
        const bool cyc = s.periodic[0];
        beginPatch(s.names[side], cyc ? PatchType::Cyclic : s.types[side]);
        if (cyc) m.patches.back().neighbourPatch = s.names[1 - side];
        const int i = side == 0 ? 0 : nx;
        for (int k = 0; k < nz; ++k)
            for (int j = 0; j < ny; ++j) {
                if (side == 0)
                    m.addFace({P(i, j, k), P(i, j, k + 1), P(i, j + 1, k + 1), P(i, j + 1, k)});
                else
                    m.addFace({P(i, j, k), P(i, j + 1, k), P(i, j + 1, k + 1), P(i, j, k + 1)});
                m.owner.push_back(C(side == 0 ? 0 : nx - 1, j, k));
            }
        endPatch();
    }
    // yMin / yMax
    for (int side = 0; side < 2; ++side) {
        const bool cyc = s.periodic[1];
        beginPatch(s.names[2 + side], cyc ? PatchType::Cyclic : s.types[2 + side]);
        if (cyc) m.patches.back().neighbourPatch = s.names[3 - side];
        const int j = side == 0 ? 0 : ny;
        for (int k = 0; k < nz; ++k)
            for (int i = 0; i < nx; ++i) {
                if (side == 0)
                    m.addFace({P(i, j, k), P(i + 1, j, k), P(i + 1, j, k + 1), P(i, j, k + 1)});
                else
                    m.addFace({P(i, j, k), P(i, j, k + 1), P(i + 1, j, k + 1), P(i + 1, j, k)});
                m.owner.push_back(C(i, side == 0 ? 0 : ny - 1, k));
            }
        endPatch();
    }
    // zMin / zMax（二维时合并为一个 empty 面组 frontAndBack）
    if (s.twoD) beginPatch("frontAndBack", PatchType::Empty);
    for (int side = 0; side < 2; ++side) {
        const bool cyc = s.periodic[2] && !s.twoD;
        if (!s.twoD) {
            beginPatch(s.names[4 + side], cyc ? PatchType::Cyclic : s.types[4 + side]);
            if (cyc) m.patches.back().neighbourPatch = s.names[5 - side];
        }
        const int k = side == 0 ? 0 : nz;
        for (int j = 0; j < ny; ++j)
            for (int i = 0; i < nx; ++i) {
                if (side == 0)
                    m.addFace({P(i, j, k), P(i, j + 1, k), P(i + 1, j + 1, k), P(i + 1, j, k)});
                else
                    m.addFace({P(i, j, k), P(i + 1, j, k), P(i + 1, j + 1, k), P(i, j + 1, k)});
                m.owner.push_back(C(i, j, side == 0 ? 0 : nz - 1));
            }
        if (!s.twoD) endPatch();
    }
    if (s.twoD) endPatch();
    return m;
}

} // namespace cfd
