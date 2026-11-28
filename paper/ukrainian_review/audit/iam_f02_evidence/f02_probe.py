"""Probe the exported P21 arrays (read-only): numbering of crack-face lists, mirror map semantics, core structure."""
import json, os, sys
import numpy as np
sys.stdout.reconfigure(encoding='utf-8')
RD = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'raw')
IX = {x['file'][:-4]: x for x in json.load(open(os.path.join(RD, 'RAW_INDEX.json'), encoding='utf-8'))}


def raw(n):
    x = IX[n]; a = np.fromfile(os.path.join(RD, x['file']), dtype='<f8'); return a.reshape([int(s) for s in x['shape']], order='F')


X3, T3, XL3, X6, T6, XL6 = raw('X3_mm'), raw('T3').astype(int) - 1, raw('XL3_mm'), raw('X6_mm'), raw('T6').astype(int) - 1, raw('XL6_mm')
up, lo, tip = raw('upperNodes').ravel().astype(int) - 1, raw('lowerNodes').ravel().astype(int) - 1, int(raw('tipNode').ravel()[0]) - 1
pe, ee, ed = raw('pairedElementIDs').ravel().astype(int) - 1, raw('exteriorElementIDs').ravel().astype(int) - 1, raw('EDIsupportIDs').ravel().astype(int) - 1
pn, mir = raw('pairedNodeIDs').ravel().astype(int) - 1, raw('pairedMirrorLocal').ravel()
print('X3', X3.shape, 'T3', T3.shape, 'X6', X6.shape, 'T6', T6.shape)
print('upper/lower', len(up), len(lo), 'tip', tip, 'max up idx', up.max(), 'paired elems', len(pe), 'exterior', len(ee), 'EDI', len(ed))
for name, X in (('XL3', XL3), ('XL6', XL6)):
    if up.max() < len(X):
        print(name, 'upper y range', X[up, 1].min(), X[up, 1].max(), 'x range', X[up, 0].min(), X[up, 0].max(), 'tip', X[tip])
print('paired node ids', len(pn), pn.min(), pn.max(), 'mirror', len(mir), mir.min(), mir.max())
# mirror semantics: try mir as local index into pn
for cand, X in (('XL3', XL3), ('XL6', XL6)):
    if pn.max() < len(X):
        P = X[pn]
        if mir.max() <= len(pn):
            Q = X[pn[mir.astype(int) - 1]]
            print(cand, 'mirror as index into pairedNodeIDs: max|x-xm|', np.max(np.abs(P[:, 0] - Q[:, 0])), 'max|y+ym|', np.max(np.abs(P[:, 1] + Q[:, 1])))
        if mir.max() <= len(X):
            Q = X[mir.astype(int) - 1]
            print(cand, 'mirror as global index: max|x-xm|', np.max(np.abs(P[:, 0] - Q[:, 0])), 'max|y+ym|', np.max(np.abs(P[:, 1] + Q[:, 1])))
# core structure
c = XL3[T3[pe]].mean(axis=1); r = np.linalg.norm(c, axis=1)
a, b, cc = XL3[T3[pe, 0]], XL3[T3[pe, 1]], XL3[T3[pe, 2]]
A = 0.5 * ((b[:, 0] - a[:, 0]) * (cc[:, 1] - a[:, 1]) - (cc[:, 0] - a[:, 0]) * (b[:, 1] - a[:, 1])); h = np.sqrt(2 * np.abs(A))
for lo_, hi_ in ((0, .05), (.05, .2), (.2, .4), (.4, 1), (1, 2), (2, 2.6), (2.6, 3.1)):
    m = (r >= lo_) & (r < hi_); print(f'r in [{lo_},{hi_}): n={m.sum()} size~{h[m].mean() if m.any() else 0:.4f} mm')
print('core max radius of nodes', np.max(np.linalg.norm(XL3[np.unique(T3[pe])], axis=1)))
print('min T3 area sign all', np.all(A > 0))
