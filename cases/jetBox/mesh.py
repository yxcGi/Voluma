#!/usr/bin/env python3
"""封闭箱体射流二维网格（Gmsh 三角形）：箱 2 × 1，左壁中部 0.4 < y < 0.6 为 inlet，右壁为 outlet，其余为 wall。
    python3 mesh.py [h=0.025] → jetBox.msh
"""
import sys

import gmsh

h = float(sys.argv[1]) if len(sys.argv) > 1 else 0.025
gmsh.initialize()
gmsh.option.setNumber('General.Terminal', 0)
gmsh.model.add('jetBox')
geo = gmsh.model.geo
P = [geo.addPoint(x, y, 0, h) for x, y in ((0, 0), (2, 0), (2, 1), (0, 1), (0, 0.6), (0, 0.4))]
bottom = geo.addLine(P[0], P[1])
outlet = geo.addLine(P[1], P[2])
top = geo.addLine(P[2], P[3])
leftTop = geo.addLine(P[3], P[4])
inlet = geo.addLine(P[4], P[5])
leftBot = geo.addLine(P[5], P[0])
loop = geo.addCurveLoop([bottom, outlet, top, leftTop, inlet, leftBot])
s = geo.addPlaneSurface([loop])
geo.synchronize()
gmsh.model.mesh.setTransfiniteCurve(inlet, 9)
gmsh.model.addPhysicalGroup(1, [inlet], name='inlet')
gmsh.model.addPhysicalGroup(1, [outlet], name='outlet')
gmsh.model.addPhysicalGroup(1, [bottom, top, leftTop, leftBot], name='wall')
gmsh.model.addPhysicalGroup(2, [s], name='fluid')
gmsh.model.mesh.generate(2)
gmsh.write('jetBox.msh')
gmsh.finalize()
print('jetBox.msh written')
