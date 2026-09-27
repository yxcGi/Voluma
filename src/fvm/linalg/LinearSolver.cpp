#include "fvm/linalg/LinearSolver.h"

#include "fvm/linalg/GAMG.h"

#include <algorithm>
#include <cmath>
#include <iostream>
#include <map>
#include <memory>
#include <stdexcept>

namespace cfd {

namespace {

scalar gSumProd(const std::vector<scalar>& a, const std::vector<scalar>& b, label n) {
    par::SumAcc s;
    for (label c = 0; c < n; ++c) s.add(a[c] * b[c]);
    return s.allReduce();
}

scalar gSumMag(const std::vector<scalar>& a, label n) {
    par::SumAcc s;
    for (label c = 0; c < n; ++c) s.add(std::abs(a[c]));
    return s.allReduce();
}

// 预条件器使用的面顺序：两侧均为自有单元的内部面，(l,u)=(min,max)，按 (u,l) 升序
struct PrecondFaces {
    std::vector<label> face, l, u;
    std::vector<bool> ownerIsL;
};

PrecondFaces buildPrecondFaces(const Mesh& m) {
    PrecondFaces pf;
    std::vector<std::array<label, 3>> tmp;
    for (label f = 0; f < m.nInternalFaces(); ++f) {
        const label o = m.owner()[f], n = m.neighbour()[f];
        if (o >= m.nCells() || n >= m.nCells() || o == n) continue;
        tmp.push_back({std::max(o, n), std::min(o, n), f});
    }
    std::sort(tmp.begin(), tmp.end());
    for (auto& t : tmp) {
        pf.face.push_back(t[2]);
        pf.u.push_back(t[0]);
        pf.l.push_back(t[1]);
        pf.ownerIsL.push_back(m.owner()[t[2]] == t[1]);
    }
    return pf;
}

const PrecondFaces& precondFaces(const Mesh& m) {
    return m.cached<PrecondFaces>("precondFaces", [&] { return buildPrecondFaces(m); });
}

class Preconditioner {
public:
    Preconditioner(const LduSystem& s, const std::string& type) : s_(s), type_(type) {
        const label n = s.mesh.nCells();
        rD_.assign(n, 0.0);
        if (type_ == "none") return;
        if (type_ == "GAMG") {
            gamg_ = std::make_unique<GamgPreconditioner>(s);
            return;
        }
        for (label c = 0; c < n; ++c) rD_[c] = s.diag[c];
        if (type_ == "DIC" || type_ == "DILU") {
            const auto& pf = precondFaces(s.mesh);
            for (std::size_t k = 0; k < pf.face.size(); ++k) {
                const label f = pf.face[k];
                const scalar Alu = pf.ownerIsL[k] ? s.upper[f] : s.lower[f];
                const scalar Aul = pf.ownerIsL[k] ? s.lower[f] : s.upper[f];
                rD_[pf.u[k]] -= Alu * Aul / rD_[pf.l[k]];
            }
        } else if (type_ != "diagonal") {
            throw std::runtime_error("unknown preconditioner " + type_);
        }
        for (label c = 0; c < n; ++c) rD_[c] = 1.0 / rD_[c];
    }

    void apply(const std::vector<scalar>& r, std::vector<scalar>& w) const {
        const label n = s_.mesh.nCells();
        if (type_ == "none") {
            for (label c = 0; c < n; ++c) w[c] = r[c];
            return;
        }
        if (gamg_) {
            gamg_->apply(r, w);
            return;
        }
        for (label c = 0; c < n; ++c) w[c] = rD_[c] * r[c];
        if (type_ == "diagonal") return;
        const auto& pf = precondFaces(s_.mesh);
        const std::size_t nf = pf.face.size();
        for (std::size_t k = 0; k < nf; ++k) {
            const label f = pf.face[k];
            const scalar Aul = pf.ownerIsL[k] ? s_.lower[f] : s_.upper[f];
            w[pf.u[k]] -= rD_[pf.u[k]] * Aul * w[pf.l[k]];
        }
        for (std::size_t k = nf; k-- > 0;) {
            const label f = pf.face[k];
            const scalar Alu = pf.ownerIsL[k] ? s_.upper[f] : s_.lower[f];
            w[pf.l[k]] -= rD_[pf.l[k]] * Alu * w[pf.u[k]];
        }
    }

private:
    const LduSystem& s_;
    std::string type_;
    std::vector<scalar> rD_;
    std::unique_ptr<GamgPreconditioner> gamg_;
};

bool converged(const SolverControls& c, const SolverPerformance& p) {
    if (p.iterations < c.minIter) return false;
    return p.finalResidual < c.tolerance || (c.relTol > 0 && p.finalResidual < c.relTol * p.initialResidual);
}

} // namespace

void amul(const LduSystem& s, std::vector<scalar>& x, std::vector<scalar>& y) {
    const Mesh& m = s.mesh;
    m.halo().exchange(x);
    const label nt = m.nTotalCells();
    for (label c = 0; c < nt; ++c) y[c] = s.diag[c] * x[c];
    const auto& own = m.owner();
    const auto& nei = m.neighbour();
    const label nI = m.nInternalFaces();
    for (label f = 0; f < nI; ++f) {
        y[own[f]] += s.upper[f] * x[nei[f]];
        y[nei[f]] += s.lower[f] * x[own[f]];
    }
}

SolverPerformance solveLdu(const LduSystem& s, std::vector<scalar>& x, const SolverControls& ctrlIn) {
    const Mesh& m = s.mesh;
    SolverControls ctrl = ctrlIn;
    if (par::reproducible()) {
        // 可复现模式：只用与分区无关的算法
        if (ctrl.preconditioner != "none") ctrl.preconditioner = "diagonal";
        if (ctrl.solver == "GaussSeidel" || ctrl.solver == "symGaussSeidel") ctrl.solver = "Jacobi";
        static bool noted = false;
        if (!noted && (ctrl.preconditioner != ctrlIn.preconditioner || ctrl.solver != ctrlIn.solver)) {
            noted = true;
            std::cout << "reproducible mode: " << ctrlIn.solver << '/' << ctrlIn.preconditioner << " -> "
                      << ctrl.solver << '/' << ctrl.preconditioner << " (partition-independent)\n";
        }
    }
    const label n = m.nCells(), nt = m.nTotalCells();
    SolverPerformance perf;
    perf.solver = ctrl.solver;

    std::vector<scalar> wA(nt), rA(nt);
    amul(s, x, wA);
    for (label c = 0; c < n; ++c) rA[c] = s.b[c] - wA[c];

    // OpenFOAM 式归一化因子
    par::SumAcc xs;
    for (label c = 0; c < n; ++c) xs.add(x[c]);
    const scalar xRef = xs.allReduce() / scalar(m.nGlobalCells());
    std::vector<scalar> xr(nt, xRef), pA(nt);
    amul(s, xr, pA);
    par::SumAcc nf;
    for (label c = 0; c < n; ++c) nf.add(std::abs(wA[c] - pA[c]) + std::abs(s.b[c] - pA[c]));
    const scalar normFactor = nf.allReduce() + 1e-20;

    perf.initialResidual = gSumMag(rA, n) / normFactor;
    perf.finalResidual = perf.initialResidual;
    if (perf.initialResidual < ctrl.tolerance && ctrl.minIter <= 0) {
        perf.converged = true;
        return perf;
    }

    if (ctrl.solver == "PCG") {
        if (!s.symmetric) throw std::runtime_error("PCG requires a symmetric matrix; use PBiCGStab");
        Preconditioner P(s, ctrl.preconditioner);
        std::vector<scalar> w(nt), p(nt, 0.0);
        scalar wArA = GREAT, wArAold;
        do {
            wArAold = wArA;
            P.apply(rA, w);
            wArA = gSumProd(w, rA, n);
            if (perf.iterations == 0) {
                for (label c = 0; c < n; ++c) p[c] = w[c];
            } else {
                const scalar beta = wArA / wArAold;
                for (label c = 0; c < n; ++c) p[c] = w[c] + beta * p[c];
            }
            amul(s, p, wA);
            const scalar wApA = gSumProd(wA, p, n);
            if (std::abs(wApA) < VSMALL) break;
            const scalar alpha = wArA / wApA;
            for (label c = 0; c < n; ++c) {
                x[c] += alpha * p[c];
                rA[c] -= alpha * wA[c];
            }
            perf.finalResidual = gSumMag(rA, n) / normFactor;
            ++perf.iterations;
        } while (perf.iterations < ctrl.maxIter && !converged(ctrl, perf));
    } else if (ctrl.solver == "PBiCGStab") {
        Preconditioner P(s, ctrl.preconditioner == "DIC" ? "DILU" : ctrl.preconditioner);
        std::vector<scalar> rA0(rA), pA2(nt, 0.0), yA(nt), AyA(nt, 0.0), sA(nt), zA(nt), tA(nt);
        scalar rA0rA = 0, alpha = 0, omega = 0;
        do {
            const scalar rA0rAold = rA0rA;
            rA0rA = gSumProd(rA0, rA, n);
            if (perf.iterations == 0) {
                for (label c = 0; c < n; ++c) pA2[c] = rA[c];
            } else {
                if (std::abs(rA0rAold) < VSMALL || std::abs(omega) < VSMALL) break;
                const scalar beta = (rA0rA / rA0rAold) * (alpha / omega);
                for (label c = 0; c < n; ++c) pA2[c] = rA[c] + beta * (pA2[c] - omega * AyA[c]);
            }
            P.apply(pA2, yA);
            amul(s, yA, AyA);
            const scalar rA0AyA = gSumProd(rA0, AyA, n);
            if (std::abs(rA0AyA) < VSMALL) break;
            alpha = rA0rA / rA0AyA;
            for (label c = 0; c < n; ++c) sA[c] = rA[c] - alpha * AyA[c];
            ++perf.iterations;
            perf.finalResidual = gSumMag(sA, n) / normFactor;
            if (converged(ctrl, perf)) {
                for (label c = 0; c < n; ++c) x[c] += alpha * yA[c];
                break;
            }
            P.apply(sA, zA);
            amul(s, zA, tA);
            const scalar tAtA = gSumProd(tA, tA, n);
            omega = tAtA > VSMALL ? gSumProd(tA, sA, n) / tAtA : 0.0;
            for (label c = 0; c < n; ++c) {
                x[c] += alpha * yA[c] + omega * zA[c];
                rA[c] = sA[c] - omega * tA[c];
            }
            perf.finalResidual = gSumMag(rA, n) / normFactor;
        } while (perf.iterations < ctrl.maxIter && !converged(ctrl, perf));
    } else if (ctrl.solver == "GaussSeidel" || ctrl.solver == "symGaussSeidel" || ctrl.solver == "Jacobi") {
        const auto& off = m.cellFaceOffsets();
        const auto& cf = m.cellFaces();
        const auto& own = m.owner();
        const auto& nei = m.neighbour();
        const label nI = m.nInternalFaces();
        const bool jac = ctrl.solver == "Jacobi";
        const bool sym = ctrl.solver == "symGaussSeidel";
        std::vector<scalar> xOld;
        auto sweepCell = [&](label c, const std::vector<scalar>& xs) {
            scalar r = s.b[c];
            for (label k = off[c]; k < off[c + 1]; ++k) {
                const label f = cf[k];
                if (f >= nI) continue;
                if (own[f] == c && nei[f] == c) {
                    r -= (s.upper[f] + s.lower[f]) * xs[c];
                } else if (own[f] == c) {
                    r -= s.upper[f] * xs[nei[f]];
                } else {
                    r -= s.lower[f] * xs[own[f]];
                }
            }
            return r / s.diag[c];
        };
        do {
            for (int sw = 0; sw < ctrl.nSweeps; ++sw) {
                m.halo().exchange(x);
                if (jac) {
                    xOld = x;
                    for (label c = 0; c < n; ++c) x[c] = sweepCell(c, xOld);
                } else {
                    for (label c = 0; c < n; ++c) x[c] = sweepCell(c, x);
                    if (sym)
                        for (label c = n; c-- > 0;) x[c] = sweepCell(c, x);
                }
            }
            perf.iterations += ctrl.nSweeps;
            amul(s, x, wA);
            for (label c = 0; c < n; ++c) rA[c] = s.b[c] - wA[c];
            perf.finalResidual = gSumMag(rA, n) / normFactor;
        } while (perf.iterations < ctrl.maxIter && !converged(ctrl, perf));
    } else {
        throw std::runtime_error("unknown linear solver " + ctrl.solver);
    }
    m.halo().exchange(x);
    perf.converged = converged(ctrl, perf);
    return perf;
}

} // namespace cfd
