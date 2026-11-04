#include "fvm/io/SurfaceWriter.h"

#include "fvm/io/VtkAppender.h"

#include <filesystem>
#include <fstream>
#include <iomanip>
#include <sstream>
#include <unordered_map>

namespace cfd {

SurfaceWriter::SurfaceWriter(MeshPtr mesh, std::string dir, std::string base)
    : mesh_(std::move(mesh)), dir_(std::move(dir)), base_(std::move(base)) {
    if (par::master()) std::filesystem::create_directories(dir_ + "/VTK");
    par::barrier();
}

void SurfaceWriter::write(scalar time) {
    const Mesh& m = *mesh_;
    const VtkPiece& vp = m.vtk();
    const label nI = m.nInternalFaces();
    for (label p : patches_) {
        const Patch& pt = m.patches()[p];
        auto& frames = frames_[p];
        std::ostringstream frame;
        frame << base_ << '_' << pt.name << '_' << std::setw(5) << std::setfill('0') << frames.size();
        const std::string fdir = dir_ + "/VTK/" + frame.str();
        if (par::master()) std::filesystem::create_directories(fdir);
        par::barrier();

        // 本进程该 patch 的面：点重新紧致编号
        std::vector<Vec3> pts;
        std::vector<std::int64_t> conn, off;
        std::unordered_map<std::int64_t, std::int64_t> remap;
        for (label f = pt.start; f < pt.end(); ++f) {
            const label bf = f - nI;
            for (std::int64_t k = vp.bFaceOff[bf]; k < vp.bFaceOff[bf + 1]; ++k) {
                const std::int64_t gp = vp.bFaceConn[k];
                auto it = remap.find(gp);
                if (it == remap.end()) {
                    it = remap.emplace(gp, std::int64_t(pts.size())).first;
                    pts.push_back(vp.points[gp]);
                }
                conn.push_back(it->second);
            }
            off.push_back(std::int64_t(conn.size()));
        }
        vtkxml::Appender a;
        a.xml << "<Piece NumberOfPoints=\"" << pts.size() << "\" NumberOfVerts=\"0\" NumberOfLines=\"0\" "
                 "NumberOfStrips=\"0\" NumberOfPolys=\""
              << pt.size << "\">\n<Points>\n";
        a.array("Float64", "", 3, reinterpret_cast<const double*>(pts.data()), pts.size() * 3);
        a.xml << "</Points>\n<Polys>\n";
        a.array("Int64", "connectivity", 1, conn.data(), conn.size());
        a.array("Int64", "offsets", 1, off.data(), off.size());
        a.xml << "</Polys>\n<CellData>\n";
        std::vector<double> buf;
        for (const auto& it : items_) {
            const double* d = it.data();
            buf.assign(std::size_t(pt.size) * it.nc, 0.0);
            for (label i = 0; i < pt.size; ++i)
                for (int c = 0; c < it.nc; ++c) buf[std::size_t(i) * it.nc + c] = d[std::size_t(pt.start - nI + i) * it.nc + c];
            a.array("Float64", it.name, it.nc, buf.data(), buf.size());
        }
        std::vector<double> rk(pt.size, double(par::rank()));
        a.array("Float64", "procId", 1, rk.data(), rk.size());
        a.xml << "</CellData>\n</Piece>\n";
        {
            std::ofstream os(fdir + "/" + std::to_string(par::rank()) + ".vtp", std::ios::binary);
            os << "<?xml version=\"1.0\"?>\n<VTKFile type=\"PolyData\" version=\"1.0\" byte_order=\"LittleEndian\" "
                  "header_type=\"UInt64\">\n<PolyData>\n"
               << a.xml.str() << "</PolyData>\n<AppendedData encoding=\"raw\">\n_";
            os.write(a.data.data(), std::streamsize(a.data.size()));
            os << "\n</AppendedData>\n</VTKFile>\n";
        }
        if (par::master()) {
            std::ofstream os(dir_ + "/VTK/" + frame.str() + ".pvtp");
            os << "<?xml version=\"1.0\"?>\n<VTKFile type=\"PPolyData\" version=\"1.0\" byte_order=\"LittleEndian\" "
                  "header_type=\"UInt64\">\n<PPolyData GhostLevel=\"0\">\n"
                  "<PPoints><PDataArray type=\"Float64\" NumberOfComponents=\"3\"/></PPoints>\n<PCellData>\n";
            for (const auto& it : items_)
                os << "<PDataArray type=\"Float64\" Name=\"" << it.name << "\" NumberOfComponents=\""
                   << vtkxml::componentString(it.nc) << "\"/>\n";
            os << "<PDataArray type=\"Float64\" Name=\"procId\"/>\n</PCellData>\n";
            for (int r = 0; r < par::size(); ++r) os << "<Piece Source=\"" << frame.str() << '/' << r << ".vtp\"/>\n";
            os << "</PPolyData>\n</VTKFile>\n";
        }
        frames.push_back({time, "VTK/" + frame.str() + ".pvtp"});
        if (par::master()) {
            std::ofstream os(dir_ + "/" + base_ + "_" + pt.name + ".pvd");
            os << "<?xml version=\"1.0\"?>\n<VTKFile type=\"Collection\" version=\"0.1\">\n<Collection>\n";
            os << std::setprecision(12);
            for (auto& [t, f] : frames) os << "<DataSet timestep=\"" << t << "\" file=\"" << f << "\"/>\n";
            os << "</Collection>\n</VTKFile>\n";
        }
    }
    par::barrier();
}

} // namespace cfd
