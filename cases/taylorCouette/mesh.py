#!/usr/bin/env python3
"""Taylor-Couette 二维环形区域网格（Gmsh，结构化四边形）：内半径 R1 = 0.5，外半径 R2 = 1。
边界 innerWall / outerWall（名字含 wall，导入时自动为 wall 类型）。
    python3 mesh.py [nR=20] [nTheta=96]  → annulus.msh（Gmsh 4.1 ASCII）
"""
import sys

import gmsh

nR = int(sys.argv[1]) if len(sys.argv) > 1 else 20
nT = int(sys.argv[2]) if len(sys.argv) > 2 else 96
R1, R2 = 0.5, 1.0
gmsh.initialize()
gmsh.option.setNumber('General.Terminal', 0)
gmsh.model.add('annulus')
geo = gmsh.model.geo
c = geo.addPoint(0, 0, 0)
inner = [geo.addPoint(R1 * x, R1 * y, 0) for x, y in ((1, 0), (0, 1), (-1, 0), (0, -1))]
outer = [geo.addPoint(R2 * x, R2 * y, 0) for x, y in ((1, 0), (0, 1), (-1, 0), (0, -1))]
arcsIn = [geo.addCircleArc(inner[i], c, inner[(i + 1) % 4]) for i in range(4)]
arcsOut = [geo.addCircleArc(outer[i], c, outer[(i + 1) % 4]) for i in range(4)]
rad = [geo.addLine(inner[i], outer[i]) for i in range(4)]
surfs = []
for i in range(4):
    j = (i + 1) % 4
    loop = geo.addCurveLoop([arcsIn[i], rad[j], -arcsOut[i], -rad[i]])
    s = geo.addPlaneSurface([loop])
    surfs.append(s)
geo.synchronize()
for a in arcsIn + arcsOut:
    gmsh.model.mesh.setTransfiniteCurve(a, nT // 4 + 1)
for r in rad:
    gmsh.model.mesh.setTransfiniteCurve(r, nR + 1)
for s in surfs:
    gmsh.model.mesh.setTransfiniteSurface(s)
    gmsh.model.mesh.setRecombine(2, s)
gmsh.model.addPhysicalGroup(1, arcsIn, name='innerWall')
gmsh.model.addPhysicalGroup(1, arcsOut, name='outerWall')
gmsh.model.addPhysicalGroup(2, surfs, name='fluid')
gmsh.model.mesh.generate(2)
gmsh.write('annulus.msh')
gmsh.finalize()
print('annulus.msh written:', nR, 'x', nT, 'quads')
