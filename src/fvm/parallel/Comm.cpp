#include "fvm/parallel/Comm.h"

#include <chrono>
#include <cstdlib>
#include <exception>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <streambuf>

#ifdef FVM_USE_MPI
#include <mpi.h>
#endif

namespace cfd::par {

namespace {
class NullBuffer : public std::streambuf {
protected:
    int overflow(int c) override { return traits_type::not_eof(c); }
    std::streamsize xsputn(const char*, std::streamsize n) override { return n; }
};
NullBuffer nullBuffer;
std::streambuf* savedCout = nullptr;
bool ownsMpi = false;
int gRank = 0;
int gSize = 1;
bool gReproducible = false;
const auto t0 = std::chrono::steady_clock::now();
} // namespace

Environment::Environment(int& argc, char**& argv) {
#ifdef FVM_USE_MPI
    int init = 0;
    MPI_Initialized(&init);
    if (!init) {
        MPI_Init(&argc, &argv);
        ownsMpi = true;
    }
    MPI_Comm_rank(MPI_COMM_WORLD, &gRank);
    MPI_Comm_size(MPI_COMM_WORLD, &gSize);
#else
    (void)argc;
    (void)argv;
#endif
    if (const char* e = std::getenv("CFD_REPRODUCIBLE")) gReproducible = std::string(e) != "0" && std::string(e) != "";
    if (gRank != 0) savedCout = std::cout.rdbuf(&nullBuffer);
    if (gSize > 1) {
        std::set_terminate([]() {
            std::cerr << "[rank " << gRank << "] terminate called" << std::endl;
            if (auto e = std::current_exception()) {
                try {
                    std::rethrow_exception(e);
                } catch (const std::exception& ex) {
                    std::cerr << "  what(): " << ex.what() << std::endl;
                } catch (...) {
                }
            }
            abort(1);
        });
    }
}

Environment::~Environment() {
    std::cout.flush();
    if (savedCout) {
        std::cout.rdbuf(savedCout);
        savedCout = nullptr;
    }
#ifdef FVM_USE_MPI
    int fin = 0;
    MPI_Finalized(&fin);
    if (ownsMpi && !fin) MPI_Finalize();
#endif
}

int rank() { return gRank; }
int size() { return gSize; }
void setReproducible(bool on) { gReproducible = on; }
bool reproducible() { return gReproducible; }

void abort(int code) {
#ifdef FVM_USE_MPI
    int init = 0;
    MPI_Initialized(&init);
    if (init) MPI_Abort(MPI_COMM_WORLD, code);
#endif
    std::exit(code);
}

void barrier() {
#ifdef FVM_USE_MPI
    if (gSize > 1) MPI_Barrier(MPI_COMM_WORLD);
#endif
}

double wallTime() {
    return std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
}

// ================================================================ ExactSum
// 数值 x = M·2^E（M 为 ≤53 位整数，E ≥ −1074）。以 32 位为一格定点存放在 int64 数组中，
// 每次加法把 M 左移到格边界后拆成三段分别加到相邻三格；int64 余量允许 2^30 次加法后再进位。

void ExactSum::clear() {
    for (auto& b : bins_) b = 0;
    pending_ = 0;
    nonFinite_ = false;
    nonFiniteValue_ = 0;
}

void ExactSum::normalize() const {
    for (int k = 0; k < NBINS - 1; ++k) {
        const std::int64_t carry = bins_[k] >> 32;
        bins_[k] -= carry * (std::int64_t(1) << 32);
        bins_[k + 1] += carry;
    }
    pending_ = 0;
}

void ExactSum::merge(const ExactSum& o) {
    normalize();
    o.normalize();
    for (int k = 0; k < NBINS; ++k) bins_[k] += o.bins_[k];
    normalize();
    if (o.nonFinite_) {
        nonFinite_ = true;
        nonFiniteValue_ += o.nonFiniteValue_;
    }
}

namespace {
double binsToDouble(std::int64_t* b, int n) {
    // 规范化后只有最高格可能为负（代表符号）
    for (int k = 0; k < n - 1; ++k) {
        const std::int64_t c = b[k] >> 32;
        b[k] -= c * (std::int64_t(1) << 32);
        b[k + 1] += c;
    }
    double sign = 1.0;
    if (b[n - 1] < 0) {
        sign = -1.0;
        for (int k = 0; k < n; ++k) b[k] = -b[k];
        for (int k = 0; k < n - 1; ++k) {
            const std::int64_t c = b[k] >> 32;
            b[k] -= c * (std::int64_t(1) << 32);
            b[k + 1] += c;
        }
    }
    double r = 0.0;
    for (int k = n - 1; k >= 0; --k)
        if (b[k] != 0) r += std::ldexp(static_cast<double>(b[k]), 32 * k - 1074);
    return sign * r;
}
} // namespace

double ExactSum::value() const {
    if (nonFinite_) return nonFiniteValue_;
    std::int64_t b[NBINS];
    for (int k = 0; k < NBINS; ++k) b[k] = bins_[k];
    return binsToDouble(b, NBINS);
}

double ExactSum::allReduce() const {
    double out;
    allReduce(const_cast<ExactSum*>(this), 1, &out);
    return out;
}

void ExactSum::allReduce(ExactSum* sums, int n, double* out) {
    if (gSize == 1) {
        for (int i = 0; i < n; ++i) out[i] = sums[i].value();
        return;
    }
#ifdef FVM_USE_MPI
    const int stride = NBINS + 1;
    std::vector<std::int64_t> buf(std::size_t(n) * stride), res(buf.size());
    std::vector<double> nf(n), nfRes(n);
    for (int i = 0; i < n; ++i) {
        sums[i].normalize();
        for (int k = 0; k < NBINS; ++k) buf[i * stride + k] = sums[i].bins_[k];
        buf[i * stride + NBINS] = sums[i].nonFinite_ ? 1 : 0;
        nf[i] = sums[i].nonFiniteValue_;
    }
    MPI_Allreduce(buf.data(), res.data(), int(buf.size()), MPI_INT64_T, MPI_SUM, MPI_COMM_WORLD);
    MPI_Allreduce(nf.data(), nfRes.data(), n, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    for (int i = 0; i < n; ++i) {
        if (res[i * stride + NBINS] > 0)
            out[i] = nfRes[i];
        else
            out[i] = binsToDouble(&res[i * stride], NBINS);
    }
#endif
}

// ================================================================ 归约
double allMax(double v) {
#ifdef FVM_USE_MPI
    if (gSize > 1) {
        double r;
        MPI_Allreduce(&v, &r, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
        return r;
    }
#endif
    return v;
}

double allMin(double v) {
#ifdef FVM_USE_MPI
    if (gSize > 1) {
        double r;
        MPI_Allreduce(&v, &r, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);
        return r;
    }
#endif
    return v;
}

glabel allSum(glabel v) {
#ifdef FVM_USE_MPI
    if (gSize > 1) {
        std::int64_t r, x = v;
        MPI_Allreduce(&x, &r, 1, MPI_INT64_T, MPI_SUM, MPI_COMM_WORLD);
        return r;
    }
#endif
    return v;
}

void allSumInPlace(double* values, int n) {
#ifdef FVM_USE_MPI
    if (gSize > 1) MPI_Allreduce(MPI_IN_PLACE, values, n, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
#else
    (void)values;
    (void)n;
#endif
}

void allSum(ExactSum* sums, int n, double* out) { ExactSum::allReduce(sums, n, out); }

void SumAcc::allReduce(SumAcc* s, int n, double* out) {
    if (n == 0) return;
    if (s[0].exact_) {
        std::vector<ExactSum> ex;
        ex.reserve(n);
        for (int i = 0; i < n; ++i) ex.push_back(s[i].ex_);
        ExactSum::allReduce(ex.data(), n, out);
    } else {
        for (int i = 0; i < n; ++i) out[i] = s[i].plain_;
        allSumInPlace(out, n);
    }
}

double allSumTerms(const double* terms, std::size_t n) {
    SumAcc s;
    for (std::size_t i = 0; i < n; ++i) s.add(terms[i]);
    return s.allReduce();
}

// ================================================================ 点对点
#ifdef FVM_USE_MPI
namespace {
constexpr std::size_t CHUNK = std::size_t(1) << 30;
}
#endif

void sendBytes(const std::vector<char>& buf, int dest, int tag) {
#ifdef FVM_USE_MPI
    std::uint64_t n = buf.size();
    MPI_Send(&n, 1, MPI_UINT64_T, dest, tag, MPI_COMM_WORLD);
    for (std::size_t o = 0; o < buf.size(); o += CHUNK)
        MPI_Send(buf.data() + o, int(std::min(CHUNK, buf.size() - o)), MPI_BYTE, dest, tag, MPI_COMM_WORLD);
#else
    (void)buf;
    (void)dest;
    (void)tag;
    throw std::runtime_error("MPI not enabled");
#endif
}

std::vector<char> recvBytes(int source, int tag) {
#ifdef FVM_USE_MPI
    std::uint64_t n = 0;
    MPI_Recv(&n, 1, MPI_UINT64_T, source, tag, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    std::vector<char> buf(n);
    for (std::size_t o = 0; o < buf.size(); o += CHUNK)
        MPI_Recv(buf.data() + o, int(std::min(CHUNK, buf.size() - o)), MPI_BYTE, source, tag, MPI_COMM_WORLD,
                 MPI_STATUS_IGNORE);
    return buf;
#else
    (void)source;
    (void)tag;
    throw std::runtime_error("MPI not enabled");
#endif
}

void broadcastBytes(std::vector<char>& buf, int root) {
#ifdef FVM_USE_MPI
    if (gSize == 1) return;
    std::uint64_t n = buf.size();
    MPI_Bcast(&n, 1, MPI_UINT64_T, root, MPI_COMM_WORLD);
    buf.resize(n);
    for (std::size_t o = 0; o < buf.size(); o += CHUNK)
        MPI_Bcast(buf.data() + o, int(std::min(CHUNK, buf.size() - o)), MPI_BYTE, root, MPI_COMM_WORLD);
#else
    (void)buf;
    (void)root;
#endif
}

void broadcast(double* v, int n, int root) {
#ifdef FVM_USE_MPI
    if (gSize > 1) MPI_Bcast(v, n, MPI_DOUBLE, root, MPI_COMM_WORLD);
#else
    (void)v;
    (void)n;
    (void)root;
#endif
}

void broadcast(glabel* v, int n, int root) {
#ifdef FVM_USE_MPI
    if (gSize > 1) MPI_Bcast(v, n, MPI_INT64_T, root, MPI_COMM_WORLD);
#else
    (void)v;
    (void)n;
    (void)root;
#endif
}

// ================================================================ Halo
void Halo::setup(std::vector<int> nbrs, std::vector<std::vector<label>> send, std::vector<std::vector<label>> recv) {
    nbrs_ = std::move(nbrs);
    send_ = std::move(send);
    recv_ = std::move(recv);
    sendBuf_.assign(nbrs_.size(), {});
    recvBuf_.assign(nbrs_.size(), {});
}

void Halo::prepare(int nc) const {
    for (std::size_t n = 0; n < nbrs_.size(); ++n) {
        sendBuf_[n].resize(send_[n].size() * nc);
        recvBuf_[n].resize(recv_[n].size() * nc);
    }
}

void Halo::communicate(int nc) const {
    (void)nc;
#ifdef FVM_USE_MPI
    constexpr int TAG = 7001;
    const std::size_t nn = nbrs_.size();
    std::vector<MPI_Request> req(2 * nn);
    for (std::size_t n = 0; n < nn; ++n)
        MPI_Irecv(recvBuf_[n].data(), int(recvBuf_[n].size()), MPI_DOUBLE, nbrs_[n], TAG, MPI_COMM_WORLD, &req[n]);
    for (std::size_t n = 0; n < nn; ++n)
        MPI_Isend(sendBuf_[n].data(), int(sendBuf_[n].size()), MPI_DOUBLE, nbrs_[n], TAG, MPI_COMM_WORLD,
                  &req[nn + n]);
    MPI_Waitall(int(req.size()), req.data(), MPI_STATUSES_IGNORE);
#endif
}

// ================================================================ GlobalOrdering
void GlobalOrdering::setup(label nLocal, glabel nGlobal, std::vector<glabel> positionsOnMaster) {
    nLocal_ = nLocal;
    nGlobal_ = nGlobal;
    positions_ = std::move(positionsOnMaster);
    counts_.clear();
#ifdef FVM_USE_MPI
    if (gSize > 1) {
        int c = nLocal_;
        if (gRank == 0) counts_.resize(gSize);
        MPI_Gather(&c, 1, MPI_INT, counts_.data(), 1, MPI_INT, 0, MPI_COMM_WORLD);
    }
#endif
}

std::vector<double> GlobalOrdering::gatherDoubles(const std::vector<double>& local, int nc) const {
    if (gSize == 1) return local;
#ifdef FVM_USE_MPI
    std::vector<int> cnt, disp;
    std::vector<double> all;
    if (gRank == 0) {
        cnt.resize(gSize);
        disp.resize(gSize);
        int off = 0;
        for (int p = 0; p < gSize; ++p) {
            cnt[p] = counts_[p] * nc;
            disp[p] = off;
            off += cnt[p];
        }
        all.resize(off);
    }
    MPI_Gatherv(local.data(), int(local.size()), MPI_DOUBLE, all.data(), cnt.data(), disp.data(), MPI_DOUBLE, 0,
                MPI_COMM_WORLD);
    std::vector<double> ordered;
    if (gRank == 0) {
        ordered.resize(std::size_t(nGlobal_) * nc);
        for (std::size_t k = 0; k < positions_.size(); ++k)
            for (int c = 0; c < nc; ++c) ordered[positions_[k] * nc + c] = all[k * nc + c];
    }
    return ordered;
#else
    return local;
#endif
}

std::vector<double> GlobalOrdering::scatterDoubles(const std::vector<double>& global, int nc) const {
    if (gSize == 1) return global;
#ifdef FVM_USE_MPI
    std::vector<int> cnt, disp;
    std::vector<double> all;
    if (gRank == 0) {
        cnt.resize(gSize);
        disp.resize(gSize);
        int off = 0;
        for (int p = 0; p < gSize; ++p) {
            cnt[p] = counts_[p] * nc;
            disp[p] = off;
            off += cnt[p];
        }
        all.resize(off);
        for (std::size_t k = 0; k < positions_.size(); ++k)
            for (int c = 0; c < nc; ++c) all[k * nc + c] = global[positions_[k] * nc + c];
    }
    std::vector<double> local(std::size_t(nLocal_) * nc);
    MPI_Scatterv(all.data(), cnt.data(), disp.data(), MPI_DOUBLE, local.data(), int(local.size()), MPI_DOUBLE, 0,
                 MPI_COMM_WORLD);
    return local;
#else
    return global;
#endif
}

} // namespace cfd::par
