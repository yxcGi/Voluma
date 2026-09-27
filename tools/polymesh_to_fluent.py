# polyMesh(ASCII) → Fluent .msh（三维，ASCII 或二进制段）以及二进制 polyMesh；用于生成网格导入测试样例
import re, sys, struct, os
def toks(path):
    t=open(path).read()
    t=re.sub(r'/\*.*?\*/','',t,flags=re.S); t=re.sub(r'//[^\n]*','',t)
    t=t[t.index('}',t.index('FoamFile'))+1:]
    return t
def readlist(path):
    t=toks(path); i=re.search(r'(\d+)\s*\(',t); return t[i.end():]
def points(d):
    t=readlist(d+'/points'); return [tuple(map(float,m)) for m in re.findall(r'\(\s*([^\s()]+)\s+([^\s()]+)\s+([^\s()]+)\s*\)',t)]
def faces(d):
    t=readlist(d+'/faces'); return [list(map(int,m.split())) for m in re.findall(r'\d+\s*\(([^)]*)\)',t)]
def labels(d,n):
    t=readlist(d+'/'+n); return list(map(int,t[:t.index(')')].split()))
def boundary(d):
    t=toks(d+'/boundary'); out=[]
    for m in re.finditer(r'(\w+)\s*\{([^}]*)\}',t):
        body=m.group(2); g=lambda k: re.search(k+r'\s+(\S+);',body).group(1)
        out.append((m.group(1),g('type'),int(g('startFace')),int(g('nFaces'))))
    return out
def fluent(d, out, binary=False):
    P=points(d); F=faces(d); own=labels(d,'owner'); nei=labels(d,'neighbour'); B=boundary(d)
    nc=max(own+nei)+1; nI=len(nei)
    w=open(out,'wb')
    W=lambda s: w.write(s.encode())
    W('(0 "test mesh")\n(2 3)\n')
    W('(10 (0 1 %x 0 3))\n(12 (0 1 %x 0))\n(13 (0 1 %x 0))\n'%(len(P),nc,len(F)))
    if binary:
        W('(3010 (1 1 %x 1 3)('%len(P)); w.write(b''.join(struct.pack('<3d',*p) for p in P)); W(')End of Binary Section 3010)\n')
    else:
        W('(10 (1 1 %x 1 3)(\n'%len(P)); [W('%.17g %.17g %.17g\n'%p) for p in P]; W('))\n')
    W('(12 (2 1 %x 1 0))\n'%nc)
    zones=[(3,'interior','interior-1',0,nI,2)]
    bcmap={'wall':3,'patch':5,'symmetryPlane':7,'symmetry':7,'empty':3}
    for k,(n,t,s,c) in enumerate(B): zones.append((4+k, 'wall' if t in('wall','empty') else ('symmetry' if 'symm' in t else 'pressure-outlet'), n, s, c, bcmap.get(t,5)))
    first=1
    for zid,zt,zn,s,c,bc in zones:
        if c==0: continue
        # Fluent 约定：右手法向指向 c0；polyMesh 法向指向 neighbour，这里反转节点顺序并把 owner 作 c1 以模拟
        if binary:
            W('(2013 (%x %x %x %x 0)('%(zid,first,first+c-1,bc))
            for f in range(s,s+c):
                nodes=F[f]; c0=own[f]+1; c1=nei[f]+1 if f<nI else 0
                w.write(struct.pack('<%di'%(len(nodes)+3), len(nodes), *[x+1 for x in nodes], c0, c1))
            W(')End of Binary Section 2013)\n')
        else:
            W('(13 (%x %x %x %x 0)(\n'%(zid,first,first+c-1,bc))
            for f in range(s,s+c):
                nodes=F[f][::-1]; c0=own[f]+1; c1=nei[f]+1 if f<nI else 0
                W(' '.join('%x'%x for x in [len(nodes)]+[x+1 for x in nodes]+[c1 if f<nI else c0, c0 if f<nI else 0])+'\n')
            W('))\n')
        first+=c
    W('(45 (2 fluid fluid)())\n')
    for zid,zt,zn,s,c,bc in zones: W('(45 (%d %s %s)())\n'%(zid,zt,zn))
    w.close()
def binfoam(d, out):
    os.makedirs(out, exist_ok=True)
    P=points(d); F=faces(d); own=labels(d,'owner'); nei=labels(d,'neighbour')
    hdr=lambda cls,obj: ('FoamFile\n{\n    version     2.0;\n    format      binary;\n    arch        "LSB;label=32;scalar=64";\n    class       %s;\n    object      %s;\n}\n// * * * //\n\n'%(cls,obj)).encode()
    with open(out+'/points','wb') as w: w.write(hdr('vectorField','points')+b'%d\n('%len(P)+b''.join(struct.pack('<3d',*p) for p in P)+b')\n')
    off=[0]
    for f in F: off.append(off[-1]+len(f))
    flat=[x for f in F for x in f]
    with open(out+'/faces','wb') as w: w.write(hdr('faceCompactList','faces')+b'%d\n('%len(off)+struct.pack('<%di'%len(off),*off)+b')\n\n'+b'%d\n('%len(flat)+struct.pack('<%di'%len(flat),*flat)+b')\n')
    for n,L in (('owner',own),('neighbour',nei)):
        with open(out+'/'+n,'wb') as w: w.write(hdr('labelList',n)+b'%d\n('%len(L)+struct.pack('<%di'%len(L),*L)+b')\n')
    b=open(d+'/boundary').read().replace('ascii','binary')
    open(out+'/boundary','w').write(b)
if __name__=='__main__':
    d=sys.argv[1]
    fluent(d, sys.argv[2]+'.msh'); fluent(d, sys.argv[2]+'_bin.msh', True); binfoam(d, sys.argv[2]+'_binfoam')
