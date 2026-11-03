# 周期山 Re ≈ 10600 参考数据（ERCOFTAC 知识库 UFR 3-30）

来源：https://kbwiki.ercoftac.org/w/index.php/UFR_3-30_Test_Case（文件原载于 kbwiki-images.s3.amazonaws.com）。

- `exp/UFR3-30_X_10600_data_CR-0NN.dat`：Rapp (2009) PIV 实验，Re = 10600，NN = 001..010 对应站位
  x/h = 0.05, 0.5, 1, 2, 3, 4, 5, 6, 7, 8；列 `y/h, u/u_b, v/u_b, u'u'/u_b², v'v'/u_b², u'v'/u_b²`（逗号分隔，山体内为 0）。
  引用：Ch. Rapp, Experimentelle Studie der turbulenten Strömung über periodische Hügel, Dissertation, TU München, 2009。
- `les/UFR3-30_C_10595_data_MB-0NN.dat`：Breuer 等 (2009) LES（LESOCC，case 9，Re = 10595），同样十个站位；
  空格分隔，列 `y/h, u/u_b, v/u_b, u'u', v'v', u'v', k`（见文件头）。
  引用：Breuer, Peller, Rapp, Manhart, Flow over periodic hills — numerical and experimental study in a wide range of
  Reynolds numbers, Computers & Fluids 38(2), 433–457, 2009。

对比：`python3 tools/compare_hill.py output/hill --ref cases/periodicHill/reference/exp --out hill.png`。
