#!/usr/bin/env python3
"""三维单位方腔四面体网格（Gmsh OCC）：顶面 y = 1 为 movingWall，其余五面 fixedWalls。
    python3 mesh.py uniform [h=1/24] [out=cavity_tet.msh]          均匀尺寸 h（h = 1/24 约 5.8 万单元，1/32 约 13 万）
    python3 mesh.py wall [hWall=1/32] [hCore=1/14] [out=cavity_tet_wall.msh]
        壁面附近加密：距壁 0.02 以内尺寸 hWall，0.2 以外 hCore，中间线性过渡（Distance + Threshold 尺寸场）；
        默认约 6 万单元，(1/64, 1/22) 约 36 万单元。
"""
import sys

import gmsh

mode = sys.argv[1] if len(sys.argv) > 1 else 'uniform'
if mode == 'uniform':
    h = float(sys.argv[2]) if len(sys.argv) > 2 else 1 / 24
    out = sys.argv[3] if len(sys.argv) > 3 else 'cavity_tet.msh'
else:
    hWall = float(sys.argv[2]) if len(sys.argv) > 2 else 1 / 32
    hCore = float(sys.argv[3]) if len(sys.argv) > 3 else 1 / 14
    out = sys.argv[4] if len(sys.argv) > 4 else 'cavity_tet_wall.msh'
gmsh.initialize()
gmsh.option.setNumber('General.Terminal', 0)
gmsh.model.add('cavity')
occ = gmsh.model.occ
occ.addBox(0, 0, 0, 1, 1, 1)
occ.synchronize()
lid, walls = [], []
for d, t in gmsh.model.getEntities(2):
    com = occ.getCenterOfMass(d, t)
    (lid if abs(com[1] - 1) < 1e-6 else walls).append(t)
gmsh.model.addPhysicalGroup(2, lid, name='movingWall')
gmsh.model.addPhysicalGroup(2, walls, name='fixedWalls')
gmsh.model.addPhysicalGroup(3, [1], name='fluid')
if mode == 'uniform':
    gmsh.option.setNumber('Mesh.CharacteristicLengthMax', h)
    gmsh.option.setNumber('Mesh.CharacteristicLengthMin', 0.5 * h)
else:
    f = gmsh.model.mesh.field
    dist = f.add('Distance')
    f.setNumbers(dist, 'SurfacesList', lid + walls)
    f.setNumber(dist, 'Sampling', 200)
    thr = f.add('Threshold')
    f.setNumber(thr, 'InField', dist)
    f.setNumber(thr, 'SizeMin', hWall)
    f.setNumber(thr, 'SizeMax', hCore)
    f.setNumber(thr, 'DistMin', 0.02)
    f.setNumber(thr, 'DistMax', 0.2)
    f.setAsBackgroundMesh(thr)
    gmsh.option.setNumber('Mesh.MeshSizeExtendFromBoundary', 0)
    gmsh.option.setNumber('Mesh.MeshSizeFromPoints', 0)
    gmsh.option.setNumber('Mesh.MeshSizeFromCurvature', 0)
gmsh.option.setNumber('Mesh.Optimize', 1)
gmsh.option.setNumber('Mesh.OptimizeNetgen', 1)
gmsh.model.mesh.generate(3)
gmsh.write(out)
gmsh.finalize()
print(out, 'written')
