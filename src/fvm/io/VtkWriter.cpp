#include "fvm/io/VtkWriter.h"

#include <cstdio>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <sstream>

namespace cfd {

VtkWriter::VtkWriter(MeshPtr mesh, std::string dir, std::string base)
    : mesh_(std::move(mesh)), dir_(std::move(dir)), base_(std::move(base)) {
    if (par::master()) std::filesystem::create_directories(dir_ + "/VTK");
    par::barrier();
}

namespace {
struct Appender {
    std::vector<char> data;
    std::ostringstream xml;
    template <class T> void array(const std::string& type, const std::string& name, int nc, const T* v, std::size_t n) {
        xml << "<DataArray type=\"" << type << "\"";
        if (!name.empty()) xml << " Name=\"" << name << "\"";
        if (nc > 1) xml << " NumberOfComponents=\"" << nc << "\"";
        xml << " format=\"appended\" offset=\"" << data.size() << "\"/>\n";
        const std::uint64_t nb = sizeof(T) * n;
        const char* h = reinterpret_cast<const char*>(&nb);
        data.insert(data.end(), h, h + sizeof(nb));
        const char* p = reinterpret_cast<const char*>(v);
        data.insert(data.end(), p, p + nb);
    }
};
const char* vtkType(int nc) { return nc == 1 ? "1" : (nc == 3 ? "3" : "9"); }
} // namespace

void VtkWriter::write(scalar time) {
    const Mesh& m = *mesh_;
    const VtkPiece& vp = m.vtk();
    const int idx = int(frames_.size());
    std::ostringstream frame;
    frame << base_ << '_' << std::setw(5) << std::setfill('0') << idx;
    const std::string fdir = dir_ + "/VTK/" + frame.str();
    if (par::master()) std::filesystem::create_directories(fdir);
    par::barrier();

    Appender a;
    const label nC = m.nCells();
    a.xml << "<Piece NumberOfPoints=\"" << vp.points.size() << "\" NumberOfCells=\"" << nC << "\">\n<Points>\n";
    a.array("Float64", "", 3, reinterpret_cast<const double*>(vp.points.data()), vp.points.size() * 3);
    a.xml << "</Points>\n<Cells>\n";
    a.array("Int64", "connectivity", 1, vp.connectivity.data(), vp.connectivity.size());
    a.array("Int64", "offsets", 1, vp.offsets.data(), vp.offsets.size());
    a.array("UInt8", "types", 1, vp.types.data(), vp.types.size());
    if (vp.hasPolyhedra) {
        a.array("Int64", "faces", 1, vp.faces.data(), vp.faces.size());
        a.array("Int64", "faceoffsets", 1, vp.faceOffsets.data(), vp.faceOffsets.size());
    }
    a.xml << "</Cells>\n<CellData>\n";
    for (const auto& it : items_) a.array("Float64", it.name, it.nc, it.data(), std::size_t(nC) * it.nc);
    std::vector<double> rk(nC, double(par::rank()));
    a.array("Float64", "procId", 1, rk.data(), rk.size());
    a.xml << "</CellData>\n</Piece>\n";

    {
        std::ofstream os(fdir + "/" + std::to_string(par::rank()) + ".vtu", std::ios::binary);
        os << "<?xml version=\"1.0\"?>\n<VTKFile type=\"UnstructuredGrid\" version=\"1.0\" "
              "byte_order=\"LittleEndian\" header_type=\"UInt64\">\n<UnstructuredGrid>\n"
           << a.xml.str() << "</UnstructuredGrid>\n<AppendedData encoding=\"raw\">\n_";
        os.write(a.data.data(), std::streamsize(a.data.size()));
        os << "\n</AppendedData>\n</VTKFile>\n";
    }
    if (par::master()) {
        std::ofstream os(dir_ + "/VTK/" + frame.str() + ".pvtu");
        os << "<?xml version=\"1.0\"?>\n<VTKFile type=\"PUnstructuredGrid\" version=\"1.0\" "
              "byte_order=\"LittleEndian\" header_type=\"UInt64\">\n<PUnstructuredGrid GhostLevel=\"0\">\n"
              "<PPoints><PDataArray type=\"Float64\" NumberOfComponents=\"3\"/></PPoints>\n<PCellData>\n";
        for (const auto& it : items_)
            os << "<PDataArray type=\"Float64\" Name=\"" << it.name << "\" NumberOfComponents=\"" << vtkType(it.nc)
               << "\"/>\n";
        os << "<PDataArray type=\"Float64\" Name=\"procId\"/>\n</PCellData>\n";
        for (int r = 0; r < par::size(); ++r) os << "<Piece Source=\"" << frame.str() << '/' << r << ".vtu\"/>\n";
        os << "</PUnstructuredGrid>\n</VTKFile>\n";
    }
    frames_.push_back({time, "VTK/" + frame.str() + ".pvtu"});
    if (par::master()) {
        std::ofstream os(dir_ + "/" + base_ + ".pvd");
        os << "<?xml version=\"1.0\"?>\n<VTKFile type=\"Collection\" version=\"0.1\">\n<Collection>\n";
        os << std::setprecision(12);
        for (auto& [t, f] : frames_) os << "<DataSet timestep=\"" << t << "\" file=\"" << f << "\"/>\n";
        os << "</Collection>\n</VTKFile>\n";
    }
    par::barrier();
}

std::string fieldChecksum(const Mesh& m, const double* data, int nc) {
    std::vector<double> local(data, data + std::size_t(m.nCells()) * nc);
    std::vector<double> g;
    if (nc == 1)
        g = m.cellOrdering().gather(local.data());
    else if (nc == 3) {
        auto v = m.cellOrdering().gather(reinterpret_cast<const Vec3*>(local.data()));
        g.assign(reinterpret_cast<const double*>(v.data()), reinterpret_cast<const double*>(v.data()) + v.size() * 3);
    }
    std::uint64_t h = 1469598103934665603ull;
    const unsigned char* p = reinterpret_cast<const unsigned char*>(g.data());
    for (std::size_t i = 0; i < g.size() * sizeof(double); ++i) {
        h ^= p[i];
        h *= 1099511628211ull;
    }
    char buf[32];
    std::snprintf(buf, sizeof(buf), "%016llx", static_cast<unsigned long long>(h));
    return par::master() ? std::string(buf) : std::string();
}

} // namespace cfd
