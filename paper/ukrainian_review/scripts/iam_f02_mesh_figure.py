"""IAM-UA8-F02: verification and vector mesh layers for the proposed new Figure 2 (structured crack-tip mesh).

Input: figures_src/iam_f02/data (regenerated P21 mesh of the accepted M0 trajectory, see DATA_INDEX.json) and paper/data/evidence.json.
1. Independent verification (numpy only): data hashes; path = accepted vertices P0..P21; local frame of the last segment; plate/hole
   geometry; positive element areas; structured core (12678 T3, 6486 nodes, first-ring radius h0, outer radius r_core); exact mirror
   symmetry of every upper core element (node map and coordinates bitwise); crack-face node pairs (distinct ids, identical coordinates,
   never in the same element) and shared ligament nodes; EDI support inside the core and spanning the annulus [r_in, r_out];
   dimensional parameters against the production record; boundary clearance. -> figures_src/iam_f02/IAM_F02_VERIFICATION.json
2. Vector mesh layers (PDF, page = axis box at the printed size, line widths in printed pt; no text) for panels (a) global, (b) structured
   region, (c) symmetry detail; CSV point sets and a TeX constants file for the PGFPlots panels (labels in STIX Two).
usage: python -I iam_f02_mesh_figure.py"""
import gzip, hashlib, json, math, os, sys, zlib
import numpy as np

here = os.path.dirname(os.path.abspath(__file__)); ua = os.path.dirname(here); paper = os.path.dirname(ua)
SRC = os.path.join(ua, 'figures_src', 'iam_f02'); DATA = os.path.join(SRC, 'data'); OUT = os.path.join(SRC, 'layers')
os.makedirs(OUT, exist_ok=True)
MM2PT = 72 / 25.4
# ---------------------------------------------------------------- layout (printed sizes, mm); text width 166 mm
TW = 166.0
LAYOUT = {
    'a': {'slot': 0.485 * TW, 'L': 16.0, 'R': 4.5, 'win': (0.0, 300.0, -100.0, 100.0)},
    'b': {'slot': 0.485 * TW, 'L': 14.5, 'R': 2.5, 'win': (-3.4, 3.4, -3.4, 3.4)},
    'c': {'slot': 0.94 * TW, 'L': 14.5, 'R': 2.5, 'win': (-0.40, 0.16, -0.12, 0.12)},
}
for k, v in LAYOUT.items():
    x0, x1, y0, y1 = v['win']; v['W'] = v['slot'] - v['L'] - v['R']; v['H'] = v['W'] * (y1 - y0) / (x1 - x0)
problems, checks, detail = [], {}, {}


def chk(name, ok, info=None):
    checks[name] = bool(ok)
    if info is not None:
        detail[name] = info
    if not ok:
        problems.append(name)


# ---------------------------------------------------------------- data
DI = json.load(open(os.path.join(DATA, 'DATA_INDEX.json'), encoding='utf-8'))
A = {}
for x in DI['arrays']:
    b = gzip.open(os.path.join(DATA, x['file'])).read()
    chk('hash_' + x['name'], hashlib.sha256(b).hexdigest() == x['rawSHA256'])
    A[x['name']] = np.frombuffer(b, dtype='<f8').reshape([int(s) for s in x['shape']], order='F')
X3, XL3, T3 = A['X3_mm'], A['XL3_mm'], A['T3'].astype(int) - 1
up, lo, tip = A['upperNodes'].ravel().astype(int) - 1, A['lowerNodes'].ravel().astype(int) - 1, int(A['tipNode'].ravel()[0]) - 1
path = A['pathGlobal_mm']; pe = A['pairedElementIDs'].ravel().astype(int) - 1; ee = A['exteriorElementIDs'].ravel().astype(int) - 1
sup = A['EDIsupportIDs'].ravel().astype(int) - 1; pn = A['pairedNodeIDs'].ravel().astype(int) - 1; mir = A['pairedMirrorLocal'].ravel().astype(int) - 1
E = json.load(open(os.path.join(paper, 'data', 'evidence.json'), encoding='utf-8'))
Q21 = [q for q in E['qualification'] if q['segment'] == 21][0]
# ---------------------------------------------------------------- 1 verification
V = 1e3 * np.array(E['vertices_m'][:22])
chk('pathEqualsAcceptedP0toP21', path.shape == (22, 2) and float(np.max(np.abs(path - V))) <= 1e-9, float(np.max(np.abs(path - V))))
seg = np.diff(path, axis=0); L = np.linalg.norm(seg, axis=1)
chk('incrementAll4mm', float(np.max(np.abs(L - 4.0))) <= 1e-9, {'segments': len(L), 'maxDev_mm': float(np.max(np.abs(L - 4.0)))})
e1 = seg[-1] / L[-1]; e2 = np.array([-e1[1], e1[0]]); R = np.column_stack([e1, e2]); tipG = path[-1]
XLr = (X3 - tipG) @ R
chk('localFrameLastSegment', float(np.max(np.abs(XLr - XL3))) <= 1e-9 and float(np.linalg.norm(XL3[tip])) == 0.0 and float(np.linalg.norm(X3[tip] - tipG)) <= 1e-9,
    {'maxLocalDiff_mm': float(np.max(np.abs(XLr - XL3))), 'thetaLast_deg': math.degrees(math.atan2(e1[1], e1[0])), 'tip_mm': tipG.tolist()})
area = lambda X, T: 0.5 * ((X[T[:, 1], 0] - X[T[:, 0], 0]) * (X[T[:, 2], 1] - X[T[:, 0], 1]) - (X[T[:, 2], 0] - X[T[:, 0], 0]) * (X[T[:, 1], 1] - X[T[:, 0], 1]))
aG = area(X3, T3)
chk('positiveAreas', bool(np.all(aG > 0)) and len(T3) == Q21['T3_elements'] and len(X3) == Q21['T3_nodes'], {'T3': len(T3), 'nodes': len(X3), 'minArea_mm2': float(aG.min())})
hc, hr = 1e3 * np.array(E['C']['hole']['center']), 1e3 * E['C']['hole']['r']
dh = np.linalg.norm(X3 - hc, axis=1)
# the hole is the 480-gon of the frozen Stage-I carrier (vertices on R); the only other boundary nodes inside R belong to the crack-mouth
# junction of the carrier (checked: within 1 mm of the mouth P0)
onR = np.abs(dh - hr) <= 1e-9; innerExtra = (dh < hr - 1e-9)
chk('plateAndHole', float(X3[:, 0].min()) >= -1e-9 and float(X3[:, 0].max()) <= 300 + 1e-9 and float(np.abs(X3[:, 1]).max()) <= 100 + 1e-9
    and int(onR.sum()) == E['C']['hole']['npoly'] and bool(np.all(np.linalg.norm(X3[innerExtra] - path[0], axis=1) <= 1.0)) and float(dh.min()) >= hr - 0.002,
    {'xRange': [float(X3[:, 0].min()), float(X3[:, 0].max())], 'yRange': [float(X3[:, 1].min()), float(X3[:, 1].max())], 'nodesOnR': int(onR.sum()),
     'nodesInsideR_nearMouth': int(innerExtra.sum()), 'minDistToHoleCentre_mm': float(dh.min())})
clear = min(300 - tipG[0], 100 - abs(tipG[1]), tipG[0], 100 + tipG[1])
chk('boundaryClearance', abs(clear - 1e3 * Q21['physical_clearance_m']) <= 1e-6, {'clearance_mm': clear, 'record_mm': 1e3 * Q21['physical_clearance_m']})
# structured core
nUp = len(mir); nEl = len(pe); inv = -np.ones(len(X3), int); inv[pn] = np.arange(len(pn))
Zc = XL3[pn]; rc = np.linalg.norm(Zc, axis=1)
h0 = 1e3 * Q21['hTip_m']; ri, ro, rp = 1e3 * Q21['rInner_m'], 1e3 * Q21['rOuter_m'], 1e3 * Q21['rCore_m']
chk('coreCounts', nEl == Q21['core_T3_elements'] == 12678 and len(pn) == 6486 and nUp == 3273 and len(sup) == Q21['EDI_elements'] == 11316 and len(ee) == Q21['exterior_T3_elements'],
    {'coreElements': nEl, 'coreNodes': len(pn), 'upperNodes': nUp, 'EDIsupport': len(sup), 'exterior': len(ee)})
rFirst = float(np.min(rc[rc > 0]))
chk('coreRadii', abs(rFirst - h0) <= 1e-12 and abs(float(rc.max()) - rp) <= 1e-12 and abs(h0 - 0.0270123254) <= 1e-12 and ri == 0.4 and abs(ro - 2.6) <= 1e-12 and rp == 3.0,
    {'firstRing_mm': rFirst, 'outerRing_mm': float(rc.max()), 'h0_record_mm': h0, 'rIn_mm': ri, 'rOut_mm': ro, 'rCore_mm': rp, 'increment_mm': 1e3 * E['increment_m']})
ups = np.arange(nUp); low = mir
chk('mirrorNodesBitwise', bool(np.all(Zc[low, 0] == Zc[ups, 0])) and bool(np.all(Zc[low, 1] == -Zc[ups, 1])) and bool(np.all(Zc[ups, 1] >= 0)))
coreTri = inv[T3[pe]]; assert np.all(coreTri >= 0)
key = {tuple(sorted(t)): i for i, t in enumerate(coreTri.tolist())}
upperEl = [i for i, t in enumerate(coreTri) if np.all(t < nUp)]
mis = 0
mirrorOf = {}
for i in upperEl:
    m = tuple(sorted(mir[coreTri[i]].tolist()))
    if m in key:
        mirrorOf[i] = key[m]
    else:
        mis += 1
lowerOnly = nEl - len(upperEl)
selfMirror = sum(1 for i, j in mirrorOf.items() if i == j)
chk('mirrorElementsComplete', mis == 0 and len(upperEl) == nEl // 2 and selfMirror == 0 and len(set(mirrorOf.values())) == len(upperEl),
    {'upperElements': len(upperEl), 'mirroredFound': len(mirrorOf), 'missing': mis, 'selfMirror': selfMirror})
axis = np.where(Zc[:nUp, 1] == 0)[0]; face = axis[Zc[axis, 0] < 0]; lig = axis[Zc[axis, 0] >= 0]
chk('crackFacePairsAndSharedLigament', bool(np.all(mir[face] != face)) and bool(np.all(mir[lig] == lig)) and bool(np.all(Zc[mir[face]] == Zc[face])),
    {'corePairedFaceNodes': int(len(face)), 'coreSharedLigamentNodes': int(len(lig))})
fu = pn[face]; fl = pn[mir[face]]
elemSets = [set(t) for t in T3.tolist()]
bad = sum(1 for u, l in zip(fu, fl) for s in elemSets if u in s and l in s) if False else None
nodeEl = {}
for e, t in enumerate(T3.tolist()):
    for n in t:
        nodeEl.setdefault(n, set()).add(e)
together = sum(1 for u, l in zip(fu.tolist(), fl.tolist()) if nodeEl.get(u, set()) & nodeEl.get(l, set()))
upSide = all(np.mean(XL3[T3[e]], axis=0)[1] > 0 for u in fu.tolist() for e in nodeEl[u])
loSide = all(np.mean(XL3[T3[e]], axis=0)[1] < 0 for l in fl.tolist() for e in nodeEl[l])
chk('crackFacesTopologicallyDistinct', together == 0 and upSide and loSide and len(up) == len(lo) and len(set(up.tolist()) & set(lo.tolist())) == 1
    and float(np.max(np.abs(XL3[up] - XL3[lo]))) == 0.0 if len(up) == len(lo) else False,
    {'globalFacePairs': int(len(up)), 'sharedFaceNodes(tip)': int(len(set(up.tolist()) & set(lo.tolist()))), 'pairsSharingAnElement': together})
supSet = set(sup.tolist())
rN = np.linalg.norm(XL3[T3], axis=2)
spans = [(float(rN[e].min()), float(rN[e].max())) for e in sup]
coreSet = set(pe.tolist())
# production EDI support = core elements where the nodal weight q (1 for r <= r_in, linear to 0 at r_out) is not identically zero,
# i.e. elements with a node at r < r_out; the elements where q varies (gradient non-zero) exclude those with every node at r <= r_in
disk = set(e for e in pe.tolist() if rN[e].min() < ro - 1e-12)
qConst1 = set(e for e in sup.tolist() if rN[e].max() <= ri + 1e-12)
qVar = sorted(supSet - qConst1)
chk('EDIsupportEqualsQnonzeroDisk', supSet <= coreSet and supSet == disk, {'support': len(sup), 'qNonzeroDisk': len(disk), 'maxOuterRadius_mm': max(b for a, b in spans)})
chk('EDIvaryingQelements', len(qVar) == 10278, {'elementsWithVaryingQ': len(qVar), 'elementsWithQ1': len(qConst1), 'codeSkipConstantFingerprint': 10278})
# ---------------------------------------------------------------- 2 layers
def edges_of(T):
    e = np.vstack([T[:, [0, 1]], T[:, [1, 2]], T[:, [2, 0]]]); e.sort(axis=1)
    return np.unique(e, axis=0)


def in_window(P, Q, win, pad):
    x0, x1, y0, y1 = win
    return ~(((P[:, 0] < x0 - pad) & (Q[:, 0] < x0 - pad)) | ((P[:, 0] > x1 + pad) & (Q[:, 0] > x1 + pad)) | ((P[:, 1] < y0 - pad) & (Q[:, 1] < y0 - pad)) | ((P[:, 1] > y1 + pad) & (Q[:, 1] > y1 + pad)))


class Layer:
    def __init__(self, name, panel):
        self.name, self.p = name, LAYOUT[panel]; self.ops = []
        x0, x1, y0, y1 = self.p['win']; self.sx = self.p['W'] * MM2PT / (x1 - x0); self.sy = self.p['H'] * MM2PT / (y1 - y0)
        assert abs(self.sx / self.sy - 1) < 1e-12

    def m(self, P):
        x0, _, y0, _ = self.p['win']; return (P[:, 0] - x0) * self.sx, (P[:, 1] - y0) * self.sy

    def lines(self, P, Q, gray, w):
        px, py = self.m(P); qx, qy = self.m(Q)
        self.ops.append(f'{gray:.3f} G {w:.3f} w 1 J 1 j')
        self.ops.extend(f'{a:.3f} {b:.3f} m {c:.3f} {d:.3f} l' for a, b, c, d in zip(px, py, qx, qy))
        self.ops.append('S')

    def fill(self, X, T, gray):
        self.ops.append(f'{gray:.3f} g')
        for t in T:
            x, y = self.m(X[t])
            self.ops.append(f'{x[0]:.3f} {y[0]:.3f} m {x[1]:.3f} {y[1]:.3f} l {x[2]:.3f} {y[2]:.3f} l h')
        self.ops.append('f')

    def write(self):
        W, H = self.p['W'] * MM2PT, self.p['H'] * MM2PT
        content = (f'q 0 0 {W:.4f} {H:.4f} re W n\n' + '\n'.join(self.ops) + '\nQ\n').encode('ascii')
        z = zlib.compress(content, 9)
        objs = [b'<< /Type /Catalog /Pages 2 0 R >>', b'<< /Type /Pages /Kids [3 0 R] /Count 1 >>',
                f'<< /Type /Page /Parent 2 0 R /MediaBox [0 0 {W:.4f} {H:.4f}] /Contents 4 0 R /Resources << >> >>'.encode(),
                b'<< /Length %d /Filter /FlateDecode >>\nstream\n' % len(z) + z + b'\nendstream']
        out = bytearray(b'%PDF-1.5\n%\xe2\xe3\xcf\xd3\n'); offs = []
        for i, o in enumerate(objs, 1):
            offs.append(len(out)); out += b'%d 0 obj\n' % i + o + b'\nendobj\n'
        xref = len(out); out += b'xref\n0 %d\n0000000000 65535 f \n' % (len(objs) + 1)
        for o in offs:
            out += b'%010d 00000 n \n' % o
        out += b'trailer\n<< /Size %d /Root 1 0 R >>\nstartxref\n%d\n%%%%EOF\n' % (len(objs) + 1, xref)
        f = os.path.join(OUT, self.name + '.pdf'); open(f, 'wb').write(out)
        return {'file': os.path.relpath(f, ua).replace('\\', '/'), 'sha256': hashlib.sha256(bytes(out)).hexdigest(), 'bytes': len(out), 'W_mm': self.p['W'], 'H_mm': self.p['H']}


E3 = edges_of(T3); coreEl = np.zeros(len(T3), bool); coreEl[pe] = True
coreEdges = edges_of(T3[coreEl]); coreEdgeSet = set(map(tuple, coreEdges.tolist()))
isCore = np.array([tuple(e) in coreEdgeSet for e in E3.tolist()])
layers = {}
# (a) global, plate coordinates
La = Layer('iam_f02_a_global_mesh', 'a')
La.lines(X3[E3[:, 0]], X3[E3[:, 1]], 0.45, 0.10)
layers['a'] = La.write()
# (b) structured region, local coordinates: EDI support fill, exterior (light) and core (dark) edges, core boundary
Lb = Layer('iam_f02_b_structured_region', 'b'); win = LAYOUT['b']['win']
Lb.fill(XL3, T3[qVar], 0.86)   # interaction-integral domain: elements where q varies (r_in <= r <= r_out)
P, Qq = XL3[E3[:, 0]], XL3[E3[:, 1]]; vis = in_window(P, Qq, win, 0.2)
Lb.lines(P[vis & ~isCore], Qq[vis & ~isCore], 0.62, 0.12)
Lb.lines(P[vis & isCore], Qq[vis & isCore], 0.30, 0.10)
ce = np.vstack([T3[pe][:, [0, 1]], T3[pe][:, [1, 2]], T3[pe][:, [2, 0]]]); ce.sort(axis=1)
u_, cnt = np.unique(ce, axis=0, return_counts=True); bnd = u_[cnt == 1]
bnd = bnd[np.abs(np.linalg.norm(XL3[bnd[:, 0]], axis=1) - rp) < 1e-9]   # outer ring only (the crack faces are also free edges)
Lb.lines(XL3[bnd[:, 0]], XL3[bnd[:, 1]], 0.0, 0.90)
layers['b'] = Lb.write()
detail['coreOuterBoundaryEdges'] = int(len(bnd))
# (c) symmetry detail: three representative mirror pairs shaded, all edges in the window
Lc = Layer('iam_f02_c_symmetry_detail', 'c'); winc = LAYOUT['c']['win']
# representative upper core elements (explicit geometric definitions), each shaded with its mirror element:
#  1 the tip-fan element touching the crack face at (-h0, 0); 2 the tip-fan element touching the ligament at (+h0, 0) (shared node);
#  3 the element with an edge on the crack face nearest to x1 = -0.15 mm
tipLocal = int(np.where(rc == 0)[0][0])
onAxis = lambda n: Zc[n, 1] == 0
fan = [i for i in upperEl if tipLocal in coreTri[i]]
f1 = [i for i in fan if any(onAxis(n) and Zc[n, 0] < 0 for n in coreTri[i] if n != tipLocal)]
f2 = [i for i in fan if any(onAxis(n) and Zc[n, 0] > 0 for n in coreTri[i] if n != tipLocal)]
edgeFace = [i for i in upperEl if sum(1 for n in coreTri[i] if onAxis(n) and Zc[n, 0] < 0) == 2]
f3 = min(edgeFace, key=lambda i: abs(Zc[coreTri[i]][:, 0].mean() + 0.15))
assert len(fan) == 3 and len(f1) == 1 and len(f2) == 1, (len(fan), len(f1), len(f2))
pick = [f1[0], f2[0], f3]
grays = [0.55, 0.72, 0.86]
for g, i in zip(grays, pick):
    Lc.fill(XL3, T3[[pe[i], pe[mirrorOf[i]]]], g)
P, Qq = XL3[E3[:, 0]], XL3[E3[:, 1]]; vis = in_window(P, Qq, winc, 0.05)
Lc.lines(P[vis], Qq[vis], 0.25, 0.30)
layers['c'] = Lc.write()
selected = []
for g, i in zip(grays, pick):
    a_ = XL3[T3[pe[i]]]; b_ = XL3[T3[pe[mirrorOf[i]]]]
    selected.append({'gray': g, 'upperElement': int(pe[i]) + 1, 'mirrorElement': int(pe[mirrorOf[i]]) + 1, 'upperNodes_mm': a_.tolist(), 'mirrorNodes_mm': b_.tolist(),
                     'mirrorExact': sorted(map(tuple, b_.tolist())) == sorted((x, -y) for x, y in a_.tolist())})
chk('selectedPairsMirrorExact', all(s['mirrorExact'] for s in selected), selected)
# ---------------------------------------------------------------- point sets and TeX constants
def csv(name, rows, head):
    f = os.path.join(OUT, name + '.csv'); open(f, 'w', encoding='utf-8', newline='\n').write(head + '\n' + '\n'.join(','.join(f'{v:.12g}' for v in r) for r in rows) + '\n')


pathL = (path - tipG) @ R
csv('crack_global', path.tolist(), 'x,y'); csv('crack_local', pathL.tolist(), 'x,y')
fx = Zc[face]; lx = Zc[lig]
inC = lambda P: (P[:, 0] >= winc[0]) & (P[:, 0] <= winc[1])
csv('face_pairs_c', fx[inC(fx)].tolist(), 'x,y'); csv('ligament_c', lx[inC(lx) & (lx[:, 0] > 0)].tolist(), 'x,y')
corners = lambda w: np.array([[w[0], w[2]], [w[1], w[2]], [w[1], w[3]], [w[0], w[3]], [w[0], w[2]]])
csv('window_b_in_global', (corners(LAYOUT['b']['win']) @ R.T + tipG).tolist(), 'x,y'); csv('window_c_in_local', corners(winc).tolist(), 'x,y')
tex = ['% generated by iam_f02_mesh_figure.py -- do not edit']
d = lambda n, v: tex.append('\\def\\%s{%s}' % (n, f'{v:.10g}' if isinstance(v, float) else v))
for k, nm in (('a', 'A'), ('b', 'B'), ('c', 'C')):
    p = LAYOUT[k]; d(f'iamF{nm}AxW', f'{p["W"]:.4f}mm'); d(f'iamF{nm}AxH', f'{p["H"]:.4f}mm'); d(f'iamF{nm}L', f'{p["L"]:.4f}mm'); d(f'iamF{nm}R', f'{p["R"]:.4f}mm')
    for j, s in enumerate(('xmin', 'xmax', 'ymin', 'ymax')):
        d(f'iamF{nm}{s}', float(p['win'][j]))
d('iamFtipX', float(tipG[0])); d('iamFtipY', float(tipG[1])); d('iamFholeX', float(hc[0])); d('iamFholeY', float(hc[1])); d('iamFholeR', float(hr))
d('iamFrIn', float(ri)); d('iamFrOut', float(ro)); d('iamFrCore', float(rp)); d('iamFhTip', float(h0))
open(os.path.join(OUT, 'iam_f02_constants.tex'), 'w', encoding='utf-8', newline='\n').write('\n'.join(tex) + '\n')
res = {'schema': 'IAM-UA8-F02.Verification.v1', 'pass': not problems, 'problems': problems, 'checks': checks, 'detail': detail, 'layout_mm': LAYOUT,
       'layers': layers, 'selectedMirrorPairs': selected, 'state': {'segment': 21, 'a_mm': 84.0, 'tip_mm': tipG.tolist(),
       'KIIoverKI': [r for r in E['stateRows'] if r['segment'] == 21][0]['KII_over_KI']}}
json.dump(res, open(os.path.join(SRC, 'IAM_F02_VERIFICATION.json'), 'w', encoding='utf-8', newline='\n'), indent=1, default=float)
print('IAM_F02', 'PASS' if not problems else 'FAIL', f'{sum(checks.values())}/{len(checks)}', problems)
for k, v in layers.items():
    print(k, v['file'], v['bytes'], f"{v['W_mm']:.2f}x{v['H_mm']:.2f} mm")
