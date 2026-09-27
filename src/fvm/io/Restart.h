#pragma once
// 续算文件：按全局编号存放的二进制数组（与进程数无关，可换进程数续算）。
// 每个数组一个文件：8 字节魔数 + nGlobal + nComponents + double 数据。

#include "fvm/solvers/IncompressibleFlow.h"

#include <map>
#include <string>

namespace cfd::restart {

void writeCellArray(const std::string& file, const Mesh& m, const double* data, int nc);
// 读入自有单元并交换幽灵单元；文件不存在返回 false
bool readCellArray(const std::string& file, const Mesh& m, double* data, int nc);
void writeFaceArray(const std::string& file, const Mesh& m, const double* data, int nc);
bool readFaceArray(const std::string& file, const Mesh& m, double* data, int nc);

// 写/读不可压求解器的完整状态（U、p、φ 及旧时间层、时间信息、体力状态）。
// extra 可附带算例自己的标量（如统计累计时间）
void write(const std::string& dir, IncompressibleFlow& flow, const std::map<std::string, double>& extra = {});
bool read(const std::string& dir, IncompressibleFlow& flow, std::map<std::string, double>* extra = nullptr);

} // namespace cfd::restart
