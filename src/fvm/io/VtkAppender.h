#pragma once
// VTK XML 文件的 appended raw 二进制数据块组装（VtkWriter / SurfaceWriter 共用）：
// 每个数组在 XML 中登记 offset，数据按 [UInt64 字节数][原始字节] 依次追加。

#include <cstdint>
#include <sstream>
#include <string>
#include <vector>

namespace cfd::vtkxml {

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

inline const char* componentString(int nc) { return nc == 1 ? "1" : (nc == 3 ? "3" : "9"); }

} // namespace cfd::vtkxml
