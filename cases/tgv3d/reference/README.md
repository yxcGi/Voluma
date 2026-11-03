# Taylor-Green 涡 Re = 1600 参考数据

`TGV_Re1600.dat`：Dairay, Lamballais, Laizet & Vassilicos (2017), *Numerical dissipation vs. subgrid-scale modelling
for large eddy simulation*, J. Comput. Phys. 337, 252–274 的 DNS（Incompact3d，六阶紧致格式，512³ 均匀网格，(2π)³ 周期域）。
数据来源：Zenodo，DOI 10.5281/zenodo.2577239，许可 CC BY 4.0。

列（空格分隔，`#` 为注释）：t，E_k = ½⟨u²+v²+w²⟩，ε_t = −dE_k/dt，ε（由速度梯度计算的粘性耗散），
拟涡能 ζ（Dζ = 2νε），⟨u²⟩ ⟨v²⟩ ⟨w²⟩，以及 9 个速度梯度分量的均方值。

与 Voluma 的 `monitor.csv`（t,k,epsilon,epsilon_sgs,enstrophy）对比：`python3 tools/compare_tgv.py output/monitor.csv`。
