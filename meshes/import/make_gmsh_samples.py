# 生成 meshes/import 下的 Gmsh 测试网格（需要 pip install gmsh）：
#   mixed2D_v22.msh   单位方形，三角形+四边形混合，Gmsh 2.2
#   tet10_v41.msh     单位立方体，二阶四面体，Gmsh 4.1
#   hybrid3D_v41.msh  单位立方体，左半六面体、右半四面体，界面处金字塔，Gmsh 4.1
# Fluent 与二进制 polyMesh 样例由 tools/polymesh_to_fluent.py 从 meshes/cavity/polyMesh 转换：
#   python3 -c "exec(open('tools/polymesh_to_fluent.py').read().split('if __name__')[0]); d='meshes/cavity/polyMesh';
#               fluent(d,'cavity_fluent.msh'); fluent(d,'cavity_fluent_bin.msh',True); binfoam(d,'cavity_binary_polyMesh')"
import gmsh, sys
def cube3d(name, algo, ver, recombine=False, order=1, lc=0.12):
    gmsh.initialize(); gmsh.option.setNumber("General.Terminal", 0)
    gmsh.model.add(name)
    gmsh.model.occ.addBox(0,0,0,1,1,1); gmsh.model.occ.synchronize()
    gmsh.option.setNumber("Mesh.MeshSizeMax", lc); gmsh.option.setNumber("Mesh.MeshSizeMin", lc)
    surfs = gmsh.model.getEntities(2)
    top=[s[1] for s in surfs if abs(gmsh.model.occ.getCenterOfMass(2,s[1])[2]-1)<1e-9]
    others=[s[1] for s in surfs if s[1] not in top]
    gmsh.model.addPhysicalGroup(2, top, name="movingWall")
    gmsh.model.addPhysicalGroup(2, others, name="fixedWalls")
    gmsh.model.addPhysicalGroup(3, [1], name="fluid")
    if recombine:
        gmsh.option.setNumber("Mesh.Recombine3DAll", 1); gmsh.option.setNumber("Mesh.RecombineAll", 1)
        gmsh.option.setNumber("Mesh.Recombine3DLevel", 2)
    gmsh.option.setNumber("Mesh.Algorithm3D", algo)
    gmsh.model.mesh.generate(3)
    if order>1: gmsh.model.mesh.setOrder(order)
    gmsh.option.setNumber("Mesh.MshFileVersion", ver)
    gmsh.write(name+".msh"); 
    types, tags, _ = gmsh.model.mesh.getElements(3)
    print(name, {gmsh.model.mesh.getElementProperties(t)[0]: len(tg) for t,tg in zip(types,tags)})
    gmsh.finalize()
def sq2d(name, ver):
    gmsh.initialize(); gmsh.option.setNumber("General.Terminal", 0)
    gmsh.model.add(name)
    gmsh.model.occ.addRectangle(0,0,0,0.5,1); gmsh.model.occ.addRectangle(0.5,0,0,0.5,1)
    gmsh.model.occ.fragment([(2,1)],[(2,2)]); gmsh.model.occ.synchronize()
    gmsh.option.setNumber("Mesh.MeshSizeMax", 0.1)
    gmsh.model.mesh.setRecombine(2, 2)
    curves = gmsh.model.getEntities(1)
    top=[]; walls=[]
    for d,t in curves:
        c=gmsh.model.occ.getCenterOfMass(1,t)
        if abs(c[1]-1)<1e-9: top.append(t)
        elif abs(c[0]-0.5)<1e-9 and 0<c[1]<1: pass
        else: walls.append(t)
    gmsh.model.addPhysicalGroup(1, top, name="movingWall")
    gmsh.model.addPhysicalGroup(1, walls, name="fixedWalls")
    gmsh.model.addPhysicalGroup(2, [1,2], name="fluid")
    gmsh.model.mesh.generate(2)
    gmsh.option.setNumber("Mesh.MshFileVersion", ver)
    gmsh.write(name+".msh")
    types, tags, _ = gmsh.model.mesh.getElements(2)
    print(name, {gmsh.model.mesh.getElementProperties(t)[0]: len(tg) for t,tg in zip(types,tags)})
    gmsh.finalize()

sq2d("mixed2D_v22", 2.2)
cube3d("tet10_v41", 1, 4.1, order=2, lc=0.5)
# 混合网格：左半 transfinite 六面体，右半四面体（界面四边形 → 金字塔），外加拉伸的三棱柱/六面体层
gmsh.initialize(); gmsh.option.setNumber("General.Terminal", 0)
gmsh.model.add("hybrid")
b1=gmsh.model.occ.addBox(0,0,0,0.5,1,1); b2=gmsh.model.occ.addBox(0.5,0,0,0.5,1,1)
gmsh.model.occ.fragment([(3,b1)],[(3,b2)]); gmsh.model.occ.synchronize()
vols=gmsh.model.getEntities(3)
left=[v for d,v in vols if gmsh.model.occ.getCenterOfMass(3,v)[0]<0.5][0]
right=[v for d,v in vols if v!=left][0]
for d,c in gmsh.model.getBoundary([(3,left)],combined=False,oriented=False,recursive=True):
    if d==1: gmsh.model.mesh.setTransfiniteCurve(c, 4)
for d,s in gmsh.model.getBoundary([(3,left)],combined=False,oriented=False):
    gmsh.model.mesh.setTransfiniteSurface(s); gmsh.model.mesh.setRecombine(2,s)
gmsh.model.mesh.setTransfiniteVolume(left)
gmsh.option.setNumber("Mesh.MeshSizeMax", 0.3)
surfs=gmsh.model.getEntities(2)
top=[s for d,s in surfs if abs(gmsh.model.occ.getCenterOfMass(2,s)[2]-1)<1e-9]
iface=[s for d,s in surfs if abs(gmsh.model.occ.getCenterOfMass(2,s)[0]-0.5)<1e-9]
others=[s for d,s in surfs if s not in top and s not in iface]
gmsh.model.addPhysicalGroup(2, top, name="movingWall")
gmsh.model.addPhysicalGroup(2, others, name="fixedWalls")
gmsh.model.addPhysicalGroup(3, [left,right], name="fluid")
gmsh.model.mesh.generate(3)
gmsh.write("hybrid3D_v41.msh")
types, tags, _ = gmsh.model.mesh.getElements(3)
print({gmsh.model.mesh.getElementProperties(t)[0]: len(tg) for t,tg in zip(types,tags)})
gmsh.finalize()
