#pragma once
// MPI 通信薄封装。求解器其余部分只通过这里访问 MPI；未启用 MPI 时退化为单进程实现。
//
// 可复现归约：
//   ExactSum 用定点“超级累加器”精确求和（与求和顺序无关），因此任意进程数下
//   全局和（点积、范数、体积平均等）逐位相同，配合按全局顺序排列的局部面循环，
//   保证 np=1..N 的结果与串行逐位一致。
//
// 可复现模式默认关闭（普通 MPI_Allreduce，更快）。环境变量 CFD_REPRODUCIBLE=1 或
// setReproducible(true) 打开后：全局和改为精确求和，线性求解器把随分区变化的
// 预条件/光顺（DIC、DILU、GAMG、Gauss-Seidel）换成与分区无关的对角/Jacobi，
// 从而任意进程数结果逐位相同，用于核对并行正确性。

#include "fvm/core/Types.h"

#include <cstring>
#include <string>
#include <vector>

namespace cfd::par {

class Environment {
public:
    Environment(int& argc, char**& argv);
    ~Environment();
    Environment(const Environment&) = delete;
    Environment& operator=(const Environment&) = delete;
};

int rank();
int size();
inline bool master() { return rank() == 0; }
inline bool parallel() { return size() > 1; }
[[noreturn]] void abort(int code = 1);
void barrier();
double wallTime();

void setReproducible(bool on);
bool reproducible();

// ---------------------------------------------------------- 精确求和
class ExactSum {
public:
    static constexpr int NBINS = 72;
    ExactSum() { clear(); }
    void clear();
    inline void add(double v);
    // 与另一个累加器合并（精确）
    void merge(const ExactSum& o);
    // 本进程结果（不做全局归约）
    double value() const;
    // 全局精确归约后的结果（所有进程相同）
    double allReduce() const;
    // 多个累加器一次性全局归约
    static void allReduce(ExactSum* sums, int n, double* out);

private:
    void normalize() const;
    mutable std::int64_t bins_[NBINS];
    mutable std::int64_t pending_ = 0;
    bool nonFinite_ = false;
    double nonFiniteValue_ = 0;
};

inline void ExactSum::add(double v) {
    std::uint64_t bits;
    std::memcpy(&bits, &v, sizeof(bits));
    const int ex = int((bits >> 52) & 0x7ff);
    if (ex == 0x7ff) {  // Inf / NaN
        nonFinite_ = true;
        nonFiniteValue_ += v;
        return;
    }
    std::int64_t M = std::int64_t(bits & ((std::uint64_t(1) << 52) - 1));
    if (ex != 0) M |= std::int64_t(1) << 52;
    if (M == 0) return;
    // v = ±M·2^(max(ex,1)−1075)，定点位置 p = max(ex,1) − 1
    const int p = (ex == 0 ? 1 : ex) - 1;
    const int i = p >> 5, s = p & 31;
    // 无分支拆分（s=0 时 32−s=32 的移位对 int64 合法）
    const std::int64_t sgn = -std::int64_t(bits >> 63);  // 0 或 −1
    std::int64_t lo = (M & ((std::int64_t(1) << (32 - s)) - 1)) << s;
    const std::int64_t r = M >> (32 - s);
    std::int64_t mid = r & 0xFFFFFFFF;
    std::int64_t hi = r >> 32;
    lo = (lo ^ sgn) - sgn;
    mid = (mid ^ sgn) - sgn;
    hi = (hi ^ sgn) - sgn;
    bins_[i] += lo;
    bins_[i + 1] += mid;
    bins_[i + 2] += hi;
    if (++pending_ >= (std::int64_t(1) << 29)) normalize();
}

// 求和累加器：可复现模式下用 ExactSum，否则普通 double 累加 + MPI_Allreduce
class SumAcc {
public:
    SumAcc() : exact_(reproducible()) {}
    void add(double v) {
        if (exact_)
            ex_.add(v);
        else
            plain_ += v;
    }
    double allReduce() {
        double r;
        allReduce(this, 1, &r);
        return r;
    }
    static void allReduce(SumAcc* s, int n, double* out);

private:
    bool exact_;
    ExactSum ex_;
    double plain_ = 0;
};

// ---------------------------------------------------------- 归约
double allMax(double v);
double allMin(double v);
glabel allSum(glabel v);
// 普通/可复现求和：terms 为本进程各条目的贡献（顺序无关时可复现）
double allSumTerms(const double* terms, std::size_t n);
// 对若干个独立的和一次归约。sums[k] 是本进程第 k 个和的局部累加器
void allSum(ExactSum* sums, int n, double* out);
// 把本进程的 double 直接求和（非可复现模式下常用）
void allSumInPlace(double* values, int n);
void allMinInPlace(double* values, int n);

// ---------------------------------------------------------- 点对点
void sendBytes(const std::vector<char>& buf, int dest, int tag);
std::vector<char> recvBytes(int source, int tag);
void broadcastBytes(std::vector<char>& buf, int root = 0);
void broadcast(double* v, int n, int root = 0);
void broadcast(glabel* v, int n, int root = 0);

// ---------------------------------------------------------- 打包工具
template <class T> constexpr int nDoubles() {
    static_assert(sizeof(T) % sizeof(double) == 0, "only double-based types");
    return int(sizeof(T) / sizeof(double));
}

class Packer {
public:
    std::vector<char> buf;
    template <class T> void put(const T& v) {
        const char* p = reinterpret_cast<const char*>(&v);
        buf.insert(buf.end(), p, p + sizeof(T));
    }
    template <class T> void putVec(const std::vector<T>& v) {
        put<std::uint64_t>(v.size());
        const char* p = reinterpret_cast<const char*>(v.data());
        buf.insert(buf.end(), p, p + sizeof(T) * v.size());
    }
    void putString(const std::string& s) {
        put<std::uint64_t>(s.size());
        buf.insert(buf.end(), s.begin(), s.end());
    }
};

class Unpacker {
public:
    explicit Unpacker(const std::vector<char>& b) : buf_(b) {}
    template <class T> T get() {
        T v;
        std::memcpy(&v, buf_.data() + pos_, sizeof(T));
        pos_ += sizeof(T);
        return v;
    }
    template <class T> std::vector<T> getVec() {
        const auto n = get<std::uint64_t>();
        std::vector<T> v(n);
        std::memcpy(static_cast<void*>(v.data()), buf_.data() + pos_, sizeof(T) * n);
        pos_ += sizeof(T) * n;
        return v;
    }
    std::string getString() {
        const auto n = get<std::uint64_t>();
        std::string s(buf_.data() + pos_, buf_.data() + pos_ + n);
        pos_ += n;
        return s;
    }

private:
    const std::vector<char>& buf_;
    std::size_t pos_ = 0;
};

// ---------------------------------------------------------- Halo 交换
// 每个邻居进程对应 send（自有单元局部号）与 recv（幽灵单元局部号）两张表，
// 两端按全局单元编号升序一一对应。
class Halo {
public:
    void setup(std::vector<int> nbrs, std::vector<std::vector<label>> send, std::vector<std::vector<label>> recv);
    bool empty() const { return nbrs_.empty(); }
    const std::vector<int>& neighbours() const { return nbrs_; }
    const std::vector<std::vector<label>>& sendLists() const { return send_; }
    const std::vector<std::vector<label>>& recvLists() const { return recv_; }

    template <class T> void exchange(std::vector<T>& data) const { exchange(data.data()); }
    template <class T> void exchange(T* data) const {
        if (nbrs_.empty()) return;
        constexpr int nc = nDoubles<T>();
        prepare(nc);
        for (std::size_t n = 0; n < nbrs_.size(); ++n) {
            double* b = sendBuf_[n].data();
            for (std::size_t k = 0; k < send_[n].size(); ++k) std::memcpy(b + k * nc, &data[send_[n][k]], sizeof(T));
        }
        communicate(nc);
        for (std::size_t n = 0; n < nbrs_.size(); ++n) {
            const double* b = recvBuf_[n].data();
            for (std::size_t k = 0; k < recv_[n].size(); ++k)
                std::memcpy(static_cast<void*>(&data[recv_[n][k]]), b + k * nc, sizeof(T));
        }
    }

private:
    void prepare(int nc) const;
    void communicate(int nc) const;

    std::vector<int> nbrs_;
    std::vector<std::vector<label>> send_, recv_;
    mutable std::vector<std::vector<double>> sendBuf_, recvBuf_;
};

// ---------------------------------------------------------- 全局排序（收集/分发）
// 本进程 nLocal 个条目在全局序列中的位置（仅 0 号保存全部进程的位置）。
class GlobalOrdering {
public:
    void setup(label nLocal, glabel nGlobal, std::vector<glabel> positionsOnMaster);
    label nLocal() const { return nLocal_; }
    glabel nGlobal() const { return nGlobal_; }

    // 收集到 0 号进程并按全局顺序排列；其他进程返回空
    template <class T> std::vector<T> gather(const T* local) const {
        constexpr int nc = nDoubles<T>();
        std::vector<double> packed(std::size_t(nLocal_) * nc);
        std::memcpy(packed.data(), local, sizeof(T) * nLocal_);
        auto g = gatherDoubles(packed, nc);
        std::vector<T> out;
        if (master()) {
            out.resize(nGlobal_);
            std::memcpy(static_cast<void*>(out.data()), g.data(), sizeof(T) * nGlobal_);
        }
        return out;
    }
    // 从 0 号进程的全局数组分发到各进程（local 长度 >= nLocal）
    template <class T> void scatter(const std::vector<T>& global, T* local) const {
        constexpr int nc = nDoubles<T>();
        std::vector<double> g;
        if (master()) {
            g.resize(std::size_t(nGlobal_) * nc);
            std::memcpy(g.data(), global.data(), sizeof(T) * nGlobal_);
        }
        auto l = scatterDoubles(g, nc);
        std::memcpy(static_cast<void*>(local), l.data(), sizeof(T) * nLocal_);
    }

private:
    std::vector<double> gatherDoubles(const std::vector<double>& local, int nc) const;
    std::vector<double> scatterDoubles(const std::vector<double>& global, int nc) const;
    label nLocal_ = 0;
    glabel nGlobal_ = 0;
    std::vector<int> counts_;
    std::vector<glabel> positions_;
};

} // namespace cfd::par
