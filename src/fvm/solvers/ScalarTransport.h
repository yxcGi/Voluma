#pragma once
// 标量输运（温度、浓度、被动标量）：
//
//   ∂T/∂t + ∇·(φT) − ∇·(Γ∇T) = Q + Sp·T
//
// φ 为面通量（可来自 IncompressibleFlow，也可由给定速度场生成，或为空即纯扩散）；
// Γ 可为常数或面场（湍流时 Γ = ν/Pr + ν_t/Pr_t）。稳态时 ddt 为 steadyState，
// 一次 step 做 nNonOrthCorr+1 次组装求解（非正交修正显式迭代）。

#include "fvm/discretization/Fvm.h"

#include <string>
#include <vector>

namespace cfd {

class ScalarTransport {
public:
    ScalarTransport(MeshPtr mesh, const std::string& name, scalar diffusivity)
        : mesh_(std::move(mesh)), T_(mesh_, name), gamma_(diffusivity) {}

    ScalarField& T() { return T_; }
    const ScalarField& T() const { return T_; }
    const Mesh& mesh() const { return *mesh_; }

    scalar diffusivity() const { return gamma_; }
    void setDiffusivity(scalar g) { gamma_ = g; }
    // 面上扩散系数（非空时覆盖常数 Γ）
    std::vector<scalar>& gammaFaces() { return gammaFaces_; }

    // 单位体积源项 Q（显式）与隐式系数 Sp（Sp<0 有利于稳定），长度 nCells；空表示无
    std::vector<scalar>& Q() { return Q_; }
    std::vector<scalar>& Sp() { return Sp_; }
    void setUniformSource(scalar q) { Q_.assign(mesh_->nTotalCells(), q); }

    ConvectionScheme divScheme = ConvectionScheme::parse("linearUpwind");
    DdtScheme ddtScheme = DdtScheme::Backward;
    bool steady = false;
    int nNonOrthCorr = 0;
    scalar relax = 1.0;  // 稳态欠松弛
    bool laplacianNonOrthCorr = true;
    SolverControls controls{"PBiCGStab", "DILU", 1e-10, 0.0, 1000, 0, 1};
    bool verbose = false;

    void initialize() {
        T_.updateBCs(ts_.time);
        T_.correctBoundaryConditions();
    }

    // 推进一步（稳态为一次迭代）。phi 为空指针表示无对流。
    // 当与流动求解器耦合时，传入流动的 TimeState 使时间层一致。
    SolverPerformance step(scalar dt, const std::vector<scalar>* phi, const TimeState* flowTime = nullptr) {
        const DdtScheme ddtS = steady ? DdtScheme::Steady : ddtScheme;
        if (flowTime) {
            ts_ = *flowTime;
        } else if (!steady) {
            ts_.dt0 = ts_.timeIndex > 0 ? ts_.dt : -1.0;
            ts_.dt = dt;
            ts_.time += dt;
            ++ts_.timeIndex;
        } else {
            ts_.time += 1.0;
            ++ts_.timeIndex;
        }
        if (!steady) T_.storeOld();
        T_.updateBCs(ts_.time);
        T_.correctBoundaryConditions();
        SolverPerformance first;
        for (int no = 0; no <= nNonOrthCorr; ++no) {
            FvMatrix<scalar> M = fvm::ddt(T_, ddtS, ts_);
            if (phi) M += fvm::div(*phi, T_, divScheme);
            if (gammaFaces_.empty())
                M -= fvm::laplacian(gamma_, T_, laplacianNonOrthCorr);
            else
                M -= fvm::laplacian(gammaFaces_, T_, laplacianNonOrthCorr);
            if (!Q_.empty()) M -= Q_;
            if (!Sp_.empty()) {
                std::vector<scalar> negSp(Sp_.size());
                for (std::size_t c = 0; c < Sp_.size(); ++c) negSp[c] = -Sp_[c];
                M += fvm::Sp(negSp, T_);
            }
            if (steady) M.relax(relax);
            auto perf = M.solve(controls, verbose);
            if (no == 0 && !perf.empty()) first = perf[0];
        }
        return first;
    }

    const TimeState& time() const { return ts_; }
    TimeState& time() { return ts_; }

private:
    MeshPtr mesh_;
    ScalarField T_;
    scalar gamma_;
    std::vector<scalar> gammaFaces_, Q_, Sp_;
    TimeState ts_;
};

} // namespace cfd
