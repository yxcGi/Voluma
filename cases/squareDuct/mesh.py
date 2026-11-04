#!/usr/bin/env python3
"""方形直管四面体网格（Gmsh OCC 长方体，非结构四面体）：截面 1 × 1，长 L，入口 z = 0（inlet）、出口 z = L（outlet）、四侧 wall。
    python3 mesh.py [h=1/12] [L=3] [out=duct_tet.msh]
    默认约 2.3 万四面体（最大非正交约 54°）；h = 1/18 约 7.3 万（网格收敛性检验用）。
"""
import sys

import gmsh

h = float(sys.argv[1]) if len(sys.argv) > 1 else 1 / 12
L = float(sys.argv[2]) if len(sys.argv) > 2 else 3.0
out = sys.argv[3] if len(sys.argv) > 3 else 'duct_tet.msh'
gmsh.initialize()
gmsh.option.setNumber('General.Terminal', 0)
gmsh.model.add('duct')
occ = gmsh.model.occ
occ.addBox(0, 0, 0, 1, 1, L)
occ.synchronize()
inlet, outlet, wall = [], [], []
for d, t in gmsh.model.getEntities(2):
    com = occ.getCenterOfMass(d, t)
    if abs(com[2]) < 1e-6:
        inlet.append(t)
    elif abs(com[2] - L) < 1e-6:
        outlet.append(t)
    else:
        wall.append(t)
gmsh.model.addPhysicalGroup(2, inlet, name='inlet')
gmsh.model.addPhysicalGroup(2, outlet, name='outlet')
gmsh.model.addPhysicalGroup(2, wall, name='wall')
gmsh.model.addPhysicalGroup(3, [1], name='fluid')
gmsh.option.setNumber('Mesh.CharacteristicLengthMax', h)
gmsh.option.setNumber('Mesh.CharacteristicLengthMin', 0.5 * h)
gmsh.option.setNumber('Mesh.Optimize', 1)
gmsh.option.setNumber('Mesh.OptimizeNetgen', 1)
gmsh.model.mesh.generate(3)
gmsh.write(out)
gmsh.finalize()
print(out, 'written')
