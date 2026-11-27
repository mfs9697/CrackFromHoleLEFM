"""Independent scientific data verification of the Ukrainian figures (AAM-UA7, WP3).

Independent of the reconstruction: it does not import or read extract_ua_figure_data.m, the panel .tex sources or the
style file. Expected values are recomputed here from the authoritative records:
  paper/data/evidence.json (full-precision JSON serialization of the reference evidence; the extractor reads the MAT),
  paper/data/isolated_tip_resolution.json, figures_src/data/increment_records_raw.json (raw MAT dump, no transformation),
  manuscript-quoted values (only for the retained M1 panels, which have no surviving history).
Level 1 (tables): every plotted column of figures_src/data/*.csv equals the value recomputed from the records (exact for
copies and IEEE products, <= 4 ulp for derived norms).
Level 2 (PDF geometry): for each panel the axis transform is rebuilt from the panel's own tick labels snapped to its tick
marks; every series polyline vertex and marker centre is mapped back to data units and compared, in order, with the
expected points (tolerance 2e-4 of the axis range, covering PGF parsing and PDF coordinate rounding); series identity by
stroke colour, dash and marker shape; legend sample colour/order per entry; no point outside the axis limits; fonts
embedded STIX Two only; printed sizes. Retained-layer panels: layer operators equal the source operators of the pinned
span (Figs. 8, 9a) or differ only by neutralized label/decoration paints (Figs. 1, 2); the placed Form XObject carries the
layer unchanged; the placed layer reproduces the hole contour (170, -20, R 30), the M0 reference history, the manuscript-
quoted M1 maxima, or the crack line x1 in [-3, 0] mm, x2 = 0 (mesh panels).
usage: python verify_ua_figure_data.py [--data DIR] [--figures DIR] [--only g/n,...] [--out FILE]"""
import argparse, glob, hashlib, json, math, os, re, sys
import pymupdf
sys.stdout.reconfigure(encoding='utf-8')
from pypdf import PdfReader
from pypdf.generic import ContentStream

here = os.path.dirname(os.path.abspath(__file__)); ua = os.path.dirname(here); paper = os.path.dirname(ua); root = os.path.dirname(paper)
ap = argparse.ArgumentParser()
ap.add_argument('--data', default=os.path.join(ua, 'figures_src', 'data'))
ap.add_argument('--figures', default=os.path.join(ua, 'figures'))
ap.add_argument('--retained', default=os.path.join(ua, 'figures_src', 'retained'))
ap.add_argument('--only', default='')
ap.add_argument('--out', default='')
a = ap.parse_args()
only = set(x for x in a.only.split(',') if x)
sha = lambda p: hashlib.sha256(open(p, 'rb').read()).hexdigest()
rel = lambda p: os.path.relpath(p, root).replace(os.sep, '/')
E = json.load(open(os.path.join(paper, 'data', 'evidence.json'), encoding='utf-8'))
I = json.load(open(os.path.join(paper, 'data', 'isolated_tip_resolution.json'), encoding='utf-8'))
INC = json.load(open(os.path.join(a.data, 'increment_records_raw.json'), encoding='utf-8'))
RESULT = {'schema': 'AAM-UA7.FigureDataVerification.v1', 'records': {}, 'tables': {}, 'panels': {}, 'problems': []}
for f in ('evidence.json', 'evidence_exact.mat', 'isolated_tip_resolution.json', 'increment_plot_records.mat', 'accepted_states.csv', 'cod_fits.csv', 'qualification.csv'):
    RESULT['records'][f] = sha(os.path.join(paper, 'data', f))
RESULT['records']['increment_records_raw.json(record sha)'] = INC['sha256']
if INC['sha256'] != RESULT['records']['increment_plot_records.mat']:
    RESULT['problems'].append('raw MAT dump not produced from the current increment record')


def bad(where, msg):
    RESULT['problems'].append(f'{where}: {msg}')


def csv(name):
    p = os.path.join(a.data, name + '.csv'); L = open(p, encoding='utf-8').read().strip().split('\n'); h = L[0].split(',')
    rows = []
    for l in L[1:]:
        v = l.split(','); rows.append({k: (float(x) if re.fullmatch(r'[-+0-9.eEinfINFaN]+', x) else x) for k, x in zip(h, v)})
    return rows


def ulpclose(x, y, n=4):
    return x == y or abs(x - y) <= n * max(math.ulp(x), math.ulp(y))


# ------------------------------------------------------------------ expected values from the records
S = E['stateRows']; V = [[1e3 * p[0], 1e3 * p[1]] for p in E['vertices_m']]
Q = E['qualification']; C = E['codFits']; hole = E['C']['hole']
aLS = E['interpolatedZeroLength_mm']; xLS, yLS = 1e3 * E['interpolatedZeroPoint_m'][0], 1e3 * E['interpolatedZeroPoint_m'][1]
# independent re-derivation of the interpolated zero (linear between P21 and P22 of the reference)
q21 = [r for r in S if r['segment'] == 21][0]['KII_over_KI']; q22 = [r for r in S if r['segment'] == 22][0]['KII_over_KI']
f = q21 / (q21 - q22); RESULT['interpolationCheck'] = {'fraction_recomputed': f, 'fraction_record': E['interpolationFraction'],
    'xLS_recomputed': V[21][0] + f * (V[22][0] - V[21][0]), 'xLS_record': xLS,
    'aLS_recomputed': [r for r in S if r['segment'] == 21][0]['crack_length_mm'] + f * 1e3 * E['increment_m'], 'aLS_record': aLS}
if abs(f - E['interpolationFraction']) > 1e-12 or abs(RESULT['interpolationCheck']['xLS_recomputed'] - xLS) > 1e-9 \
        or abs(RESULT['interpolationCheck']['aLS_recomputed'] - aLS) > 1e-9:
    bad('constants', 'interpolated local symmetry not reproduced from P21/P22')
Iv = sorted(I['independent'], key=lambda r: r['segment']); Vm = [[1e3 * I['mouth_m'][0], 1e3 * I['mouth_m'][1]]] + [[1e3 * r['tip_x_m'], 1e3 * r['tip_y_m']] for r in Iv]
da = math.sqrt((E['vertices_m'][1][0] - E['vertices_m'][0][0]) ** 2 + (E['vertices_m'][1][1] - E['vertices_m'][0][1]) ** 2)
late = [r for r in S if r['crack_length_mm'] >= 60]
W = [(0.04, 0.20), (0.04, 0.30), (0.08, 0.30), (0.12, 0.30)]


def cod(deg, j, key, scale):
    rows = sorted([r for r in C if r['crack_length_mm'] >= 60 and r['degree'] == deg and abs(r['lower_r_over_DeltaA'] - W[j][0]) < 1e-12
                   and abs(r['upper_r_over_DeltaA'] - W[j][1]) < 1e-12], key=lambda r: r['crack_length_mm'])
    return [(r['crack_length_mm'], scale * r[key]) for r in rows]


Fp = INC['fullPrecision']; incm = lambda k: [float(x) for x in Fp[k]]
EXP_TABLES = {
    'fig03_vertices': {'x_mm': [v[0] for v in V], 'y_mm': [v[1] for v in V], 'k': list(range(25))},
    'fig04_states': {'a_mm': [r['crack_length_mm'] for r in S], 'KI': [r['KI_unit'] for r in S], 'KII': [r['KII_unit'] for r in S],
                     'KII_1e3': [1e3 * r['KII_unit'] for r in S], 'theta_deg': [r['theta_deg'] for r in S]},
    'fig05_late': {'a_mm': [r['crack_length_mm'] for r in late], 'q': [r['KII_over_KI'] for r in late], 'q_1e3': [1e3 * r['KII_over_KI'] for r in late],
                   'dtheta_next_deg': [r['delta_theta_next_deg'] for r in late]},
    'fig06_cod': {'a_mm': [p[0] for d in (1, 2) for j in range(4) for p in cod(d, j, 'turn_error_deg', 1)],
                  'turn_err_1e3deg': [p[1] for d in (1, 2) for j in range(4) for p in cod(d, j, 'turn_error_deg', 1e3)],
                  'ratio_err_1e5': [p[1] for d in (1, 2) for j in range(4) for p in cod(d, j, 'ratio_error', 1e5)]},
    'fig07_vertices': {'x_h0_mm': [v[0] for v in V[:24]], 'y_h0_mm': [v[1] for v in V[:24]], 'x_2h0_mm': [v[0] for v in Vm], 'y_2h0_mm': [v[1] for v in Vm]},
    'fig07_metrics': {'a_mm': [1e3 * da * r['segment'] for r in S], 'dy_um': [1e6 * (m['tip_y_m'] - r['tip_y_m']) for r, m in zip(S, Iv)],
                      'dtheta_mdeg': [1e3 * (m['theta_deg'] - r['theta_deg']) for r, m in zip(S, Iv)],
                      'q_h0_1e3': [1e3 * r['KII_over_KI'] for r in S], 'q_2h0_1e3': [1e3 * m['KII_over_KI'] for m in Iv]},
    'fig09_metrics': {'a_mm': incm('crack_length_mm'), 'dy_2m4_um': incm('dy_2minus4_um'), 'dy_1m2_um': incm('dy_1minus2_um'),
                      'dth_2m4_mdeg': incm('dtheta_2minus4_mdeg'), 'dth_1m2_mdeg': incm('dtheta_1minus2_mdeg')},
    'fig09_native': {'increment_mm': incm('increment_mm'), 'a_mm': [float(x) for x in INC['native']['crack_length_mm']],
                     'turn_density_deg_per_mm': incm('turn_density_deg_per_mm')},
    'fig10_mesh': {'s_mm': [1e3 * r['path_length_m'] for r in Q if r['segment'] <= 23], 'min_angle_deg': [r['min_angle_deg'] for r in Q if r['segment'] <= 23],
                   'max_neighbor_ratio': [r['max_neighbor_ratio'] for r in Q if r['segment'] <= 23]},
    'fig10_solver': {'a_mm': [r['crack_length_mm'] for r in S if r['segment'] >= 2], 'pcg_iterations': [r['PCG_iterations'] for r in S if r['segment'] >= 2],
                     'residual_over_gate': [r['true_rel_residual'] / 5e-10 for r in S if r['segment'] >= 2]},
}
# The raw dump's fullPrecision map is keyed by column name (%.17g strings) and crack_length_mm occurs in both MAT tables, so
# the native crack length is taken from the dump's native-table column (integer-valued multiples of the increment; exact).
DERIVED = {('fig07_metrics', 'a_mm'): 4}
for name, cols in EXP_TABLES.items():
    rows = csv(name); res = {'rows': len(rows), 'columns': {}}
    for col, exp in cols.items():
        if exp is None:
            continue
        got = [r[col] for r in rows]; n = DERIVED.get((name, col), 0)
        if len(got) != len(exp):
            res['columns'][col] = {'pass': False, 'reason': f'{len(got)} rows vs {len(exp)} expected'}; bad(name, f'{col}: row count'); continue
        if n == -1:
            ok = all(abs(g - e) <= 1e-12 * max(1, abs(e)) for g, e in zip(got, exp)); rule = 'rel <= 1e-12 (15+ digit JSON copy)'
        elif n:
            ok = all(ulpclose(g, e, n) for g, e in zip(got, exp)); rule = f'<= {n} ulp'
        else:
            ok = all(g == e for g, e in zip(got, exp)); rule = 'exact'
        maxd = max((abs(g - e) for g, e in zip(got, exp)), default=0.0)
        res['columns'][col] = {'pass': ok, 'rule': rule, 'maxAbsDiff': maxd}
        if not ok:
            bad(name, f'{col}: values differ from the record (max |d| {maxd:.3g}, rule {rule})')
    RESULT['tables'][name] = res

# ------------------------------------------------------------------ PDF geometry
COL = {'blue': (0.0, 0.447, 0.741), 'orange': (0.85, 0.325, 0.098), 'green': (0.466, 0.674, 0.188), 'purple': (0.494, 0.184, 0.556),
       'gray': (0.45, 0.45, 0.45), 'yellow': (0.929, 0.694, 0.125), 'black': (0.0, 0.0, 0.0)}


def cname(c):
    if c is None:
        return None
    return min(COL, key=lambda k: sum((x - y) ** 2 for x, y in zip(COL[k], c))) if min(sum((x - y) ** 2 for x, y in zip(COL[k], c)) for k in COL) < 1e-3 else 'other'


def num(t):
    t = t.strip().replace('−', '-')
    return float(t) if re.fullmatch(r'-?\d+(\.\d+)?', t) else None


def transform(pg):
    D = pg.get_drawings()
    frames = [d for d in D if d.get('color') is not None and cname(d['color']) == 'black' and d.get('fill') is None and d['rect'].width > 50 and d['rect'].height > 40
              and len(d['items']) <= 6]
    fr = max(frames, key=lambda d: d['rect'].width * d['rect'].height)['rect']
    mx, my = [], []
    for d in D:
        if d.get('color') is None or cname(d['color']) != 'black':
            continue
        for it in d['items']:
            if it[0] == 'l':
                p, q = it[1], it[2]
                if abs(p.x - q.x) < 1e-3 and abs(p.y - q.y) < 4.5 and min(abs(p.y - fr.y1), abs(q.y - fr.y1)) < 0.5:
                    mx.append(p.x)
                if abs(p.y - q.y) < 1e-3 and abs(p.x - q.x) < 4.5 and min(abs(p.x - fr.x0), abs(q.x - fr.x0)) < 0.5:
                    my.append(p.y)
    X, Y = [], []
    for b in pg.get_text('dict')['blocks']:
        for l in b.get('lines', []):
            t = ''.join(s['text'] for s in l['spans']); v = num(t)
            if v is None:
                continue
            bb = l['bbox']; cx, cy = (bb[0] + bb[2]) / 2, (bb[1] + bb[3]) / 2
            if 0 <= bb[1] - fr.y1 <= 12 and fr.x0 - 10 <= cx <= fr.x1 + 10:
                X.append((cx, v))
            elif 0 <= fr.x0 - bb[2] <= 12 and fr.y0 - 10 <= cy <= fr.y1 + 10:
                Y.append((cy, v))

    def fit(pairs, marks):
        pts = []
        for c, v in pairs:
            m = min(marks, key=lambda t: abs(t - c)) if marks else c
            pts.append((m if abs(m - c) < 3 else c, v))
        n = len(pts); sx = sum(p for p, _ in pts); sv = sum(v for _, v in pts); sxx = sum(p * p for p, _ in pts); sxv = sum(p * v for p, v in pts)
        b_ = (n * sxv - sx * sv) / (n * sxx - sx * sx); a_ = (sv - b_ * sx) / n
        return a_, b_, max(abs(a_ + b_ * p - v) for p, v in pts), n
    ax_, bx_, rx, nx = fit(X, mx); ay_, by_, ry, ny = fit(Y, my)
    lim = (ax_ + bx_ * fr.x0, ax_ + bx_ * fr.x1, ay_ + by_ * fr.y1, ay_ + by_ * fr.y0)
    return {'fx': lambda px: ax_ + bx_ * px, 'fy': lambda py: ay_ + by_ * py, 'frame': fr, 'lim': lim, 'tickResidual': (rx, ry), 'ticks': (nx, ny)}


def polylines(pg, color, dashed=None, minpts=2):
    out = []
    for d in pg.get_drawings():
        if d.get('fill') is not None or cname(d.get('color')) != color:
            continue
        its = d['items']
        if not its or any(i[0] != 'l' for i in its):
            continue
        dash = d.get('dashes') or ''; isd = bool(re.search(r'\[\s*[0-9.]', dash))
        if dashed is not None and isd != dashed:
            continue
        pts = [its[0][1]] + [i[2] for i in its]
        # a polyline is one connected path
        if all(abs(its[k][2].x - its[k + 1][1].x) < 1e-6 and abs(its[k][2].y - its[k + 1][1].y) < 1e-6 for k in range(len(its) - 1)) and len(pts) >= minpts:
            out.append(pts)
    return out


def mshape(d):
    """classify one drawing as a plot mark: circle (Bezier closed), square (re/qu or 4 axis-aligned lines), triangle (3 lines),
    diamond (4 oblique lines), x (2 crossing lines); None for anything else (line segments, dashes, frames)."""
    its = d['items']; k = [i[0] for i in its]
    if d.get('dashes') and re.search(r'\[\s*[0-9.]', d['dashes']):
        return None
    if set(k) == {'c'} and len(k) == 4:
        return 'circle'
    if k in (['re'], ['qu']):
        return 'square'
    if set(k) == {'l'}:
        obl = [abs(i[1].x - i[2].x) > 1e-3 and abs(i[1].y - i[2].y) > 1e-3 for i in its]
        if len(k) == 2 and all(obl):
            return 'x'
        if len(k) == 3:
            return 'triangle'
        if len(k) == 4:
            return 'diamond' if all(obl) else ('square' if not any(obl) else None)
    return None


def markers(pg, color, shape=None, size_max=8.0):
    out = []
    for d in pg.get_drawings():
        r = d['rect']
        if cname(d.get('color')) != color or r.width > size_max or r.height > size_max or r.width < 1.5 or r.height < 1.5:
            continue
        s = mshape(d)
        if s is None or (shape is not None and s != shape):
            continue
        if s in ('circle', 'square', 'x'):
            c = ((r.x0 + r.x1) / 2, (r.y0 + r.y1) / 2)
        else:   # polygon marks are symmetric about the data point: mean of the vertices
            vs = [i[1] for i in d['items']]; c = (sum(v.x for v in vs) / len(vs), sum(v.y for v in vs) / len(vs))
        out.append(c)
    return out


SHAPE = {('blue', None): 'circle', ('orange', None): 'square', ('green', None): 'triangle', ('yellow', None): 'triangle',
         ('purple', True): 'diamond', ('purple', False): 'circle', ('purple', None): 'circle'}


def match(where, got, exp, T, tolfrac=2e-4, ordered=True):
    xmin, xmax, ymin, ymax = T['lim']; tx = tolfrac * (xmax - xmin); ty = tolfrac * (ymax - ymin)
    data = [(T['fx'](p[0] if isinstance(p, tuple) else p.x), T['fy'](p[1] if isinstance(p, tuple) else p.y)) for p in got]
    if len(data) != len(exp):
        bad(where, f'{len(data)} plotted points vs {len(exp)} expected'); return {'pass': False, 'plotted': len(data), 'expected': len(exp)}
    if not ordered:
        data = sorted(data); exp = sorted(exp)
    dx = max(abs(d[0] - e[0]) for d, e in zip(data, exp)); dy = max(abs(d[1] - e[1]) for d, e in zip(data, exp))
    ok = dx <= tx and dy <= ty
    inside = all(xmin - tx <= d[0] <= xmax + tx and ymin - ty <= d[1] <= ymax + ty for d in data)
    if not ok:
        bad(where, f'coordinates differ (max |dx| {dx:.3g} > {tx:.3g} or |dy| {dy:.3g} > {ty:.3g})')
    if not inside:
        bad(where, 'a plotted point lies outside the axis limits (clipped data)')
    return {'pass': ok and inside, 'points': len(data), 'maxAbsDx': dx, 'maxAbsDy': dy, 'tolX': tx, 'tolY': ty}


def series(pg, T, where, color, exp, dashed=None, marker=True, shape=None):
    shape = shape or SHAPE.get((color, dashed if color == 'purple' else None))
    L = polylines(pg, color, dashed)
    best = [p for p in L if len(p) == len(exp)] or L
    r = {'color': color}
    if not best:
        bad(where, f'no {color} series polyline found'); r['pass'] = False; return r
    cand = [match(where + f':{color}:line', p, exp, T) for p in best[:1]]
    r['line'] = cand[0]
    if marker:
        M = markers(pg, color, shape)
        r['markerShape'] = shape
        # keep marker centres that sit on expected points (legend sample markers lie in the legend, away from data)
        xmin, xmax, ymin, ymax = T['lim']
        Md = [m for m in M if any(abs(T['fx'](m[0]) - e[0]) <= 1e-3 * (xmax - xmin) and abs(T['fy'](m[1]) - e[1]) <= 1e-3 * (ymax - ymin) for e in exp)]
        r['markers'] = match(where + f':{color}:markers', Md, exp, T, ordered=False)
    r['pass'] = r['line']['pass'] and (not marker or r['markers']['pass'])
    return r


def texts(pg):
    return [''.join(s['text'] for s in l['spans']) for b in pg.get_text('dict')['blocks'] for l in b.get('lines', [])]


def fonts_ok(pg, where):
    F = [f[3] for f in pg.get_fonts()]
    ok = all(re.search(r'STIXTwo(Text|Math)', f) for f in F) and F != []
    sizes = sorted({round(s['size'], 2) for b in pg.get_text('dict')['blocks'] for l in b.get('lines', []) for s in l['spans']})
    if not ok:
        bad(where, f'non-STIX font embedded: {F}')
    if 9.96 not in sizes or 10.96 not in sizes:
        bad(where, f'expected 10 pt ticks and 11 pt labels, sizes {sizes}')
    for b in pg.get_text('dict')['blocks']:
        for l in b.get('lines', []):
            for s in l['spans']:
                if '�' in s['text']:
                    bad(where, 'replacement character (missing glyph) in text')
    return {'fonts': F, 'sizes_pt': sizes, 'pass': ok}


def legend(pg, where, entries):
    """each legend entry text has a sample (stroke or marker) of the expected colour on its baseline, in the given order"""
    lines = [(l['bbox'], ''.join(s['text'] for s in l['spans'])) for b in pg.get_text('dict')['blocks'] for l in b.get('lines', [])]
    D = pg.get_drawings(); res = []; prev_y = -1
    for text, color in entries:
        hit = [bb for bb, t in lines if text in t]
        if not hit:
            bad(where, f'legend text "{text}" missing'); res.append({'text': text, 'pass': False}); continue
        bb = hit[0]; yc = (bb[1] + bb[3]) / 2
        if yc < prev_y:
            bad(where, f'legend order differs at "{text}"')
        prev_y = yc
        if color is None:
            res.append({'text': text, 'pass': True}); continue
        smp = [d for d in D if cname(d.get('color')) == color and d['rect'].x1 <= bb[0] + 1 and bb[0] - d['rect'].x1 < 14
               and abs((d['rect'].y0 + d['rect'].y1) / 2 - yc) < 3.5]
        ok = bool(smp)
        if not ok:
            bad(where, f'legend sample for "{text}" is not {color}')
        res.append({'text': text, 'color': color, 'pass': ok})
    return res


def panel(gn):
    return pymupdf.open(os.path.join(a.figures, gn + '.pdf'))[0]


P = {}
pts = lambda xs, ys: list(zip(xs, ys))
X3 = EXP_TABLES['fig04_states']['a_mm']


def run(gn, fn):
    if only and gn not in only:
        return
    try:
        pg = panel(gn); T = transform(pg); r = {'transform': {'limits': T['lim'], 'ticks': T['ticks'], 'tickResidual': T['tickResidual']},
                                                    'fonts': fonts_ok(pg, gn)}
        r.update(fn(pg, T)); P[gn] = r
    except Exception as ex:
        bad(gn, f'verification error {type(ex).__name__}: {ex}'); P[gn] = {'pass': False, 'error': str(ex)}


def lim_check(gn, T, exp):
    got = T['lim']; ok = all(abs(g - e) <= 1e-3 * abs(exp[1] - exp[0] if i < 2 else exp[3] - exp[2]) for i, (g, e) in enumerate(zip(got, exp)))
    if not ok:
        bad(gn, f'axis limits {got} differ from {exp}')
    return ok


acc = pts([v[0] for v in V[:24]], [v[1] for v in V[:24]])
run('trajectory/overview', lambda pg, T: {'limits': lim_check('trajectory/overview', T, (136, 302, -54, 14)),
    'accepted': series(pg, T, 'trajectory/overview', 'blue', acc, dashed=False),
    'unsolved': series(pg, T, 'trajectory/overview', 'orange', pts([V[23][0], V[24][0]], [V[23][1], V[24][1]]), dashed=True, shape='x'),
    'lastAccepted': match('overview:P23', markers(pg, 'blue', 'square'), [tuple(V[23])], T),
    'legend': legend(pg, 'trajectory/overview', [('Прийняті стани', 'blue'), ('без розв’язку', 'orange')])})
run('trajectory/late_detail', lambda pg, T: {'limits': lim_check('trajectory/late_detail', T, (258, 302, -29, -21.3)),
    'accepted': series(pg, T, 'trajectory/late_detail', 'blue', acc[15:], dashed=False),
    'unsolved': series(pg, T, 'trajectory/late_detail', 'orange', pts([V[23][0], V[24][0]], [V[23][1], V[24][1]]), dashed=True, shape='x'),
    'lastAccepted': match('late:P23', markers(pg, 'blue', 'square'), [tuple(V[23])], T),
    'qmax': match('late:qmax', [m for m in markers(pg, 'purple', 'circle', 12)], [tuple(V[17])], T),
    'ls': match('late:ls', [m for m in markers(pg, 'green', 'diamond', 12)], [(xLS, yLS)], T),
    'labels': all(t in ' '.join(texts(pg)) for t in ('без розв’язку',))})
run('intensities/KI', lambda pg, T: {'limits': lim_check('intensities/KI', T, (4, 92, 0.2, 1.4)), 's': series(pg, T, 'KI', 'blue', pts(X3, EXP_TABLES['fig04_states']['KI']))})
run('intensities/KII', lambda pg, T: {'limits': lim_check('intensities/KII', T, (4, 92, -8, 4)), 's': series(pg, T, 'KII', 'blue', pts(X3, EXP_TABLES['fig04_states']['KII_1e3']))})
run('intensities/theta', lambda pg, T: {'limits': lim_check('intensities/theta', T, (4, 92, -4, 0.2)), 's': series(pg, T, 'theta', 'blue', pts(X3, EXP_TABLES['fig04_states']['theta_deg']))})
AL = EXP_TABLES['fig05_late']['a_mm']


def lsline(pg, T, where, color='green'):
    L = [p for p in polylines(pg, color, dashed=True) if len(p) == 2 and abs(p[0].x - p[1].x) < 1e-3]
    if not L:
        bad(where, 'a_LS line missing'); return {'pass': False}
    x = T['fx'](L[0][0].x); ok = abs(x - aLS) <= 2e-4 * (T['lim'][1] - T['lim'][0])
    if not ok:
        bad(where, f'a_LS line at {x} differs from {aLS}')
    return {'pass': ok, 'x': x}


run('late/mode_mixity', lambda pg, T: {'limits': lim_check('late/mode_mixity', T, (60, 94, -6, 4)), 's': series(pg, T, 'late/q', 'blue', pts(AL, EXP_TABLES['fig05_late']['q_1e3'])),
                                       'aLS': lsline(pg, T, 'late/mode_mixity')})
run('late/turn', lambda pg, T: {'limits': lim_check('late/turn', T, (60, 94, -0.4, 0.8)), 's': series(pg, T, 'late/turn', 'blue', pts(AL, EXP_TABLES['fig05_late']['dtheta_next_deg'])),
                                'aLS': lsline(pg, T, 'late/turn')})
CODC = ['blue', 'orange', 'green', 'purple']; CODD = [False, True, True, True]
for name, deg, key, sc in (('turn_linear', 1, 'turn_error_deg', 1e3), ('turn_quadratic', 2, 'turn_error_deg', 1e3), ('mixity_linear', 1, 'ratio_error', 1e5), ('mixity_quadratic', 2, 'ratio_error', 1e5)):
    def fn(pg, T, deg=deg, key=key, sc=sc, name=name):
        r = {'xlimits': abs(T['lim'][0] - 60) < 0.01 and abs(T['lim'][1] - 94) < 0.01}
        for j in range(4):
            r[f'window{j + 1}'] = series(pg, T, f'cod/{name}:w{j + 1}', CODC[j], cod(deg, j, key, sc), dashed=CODD[j])
        if name == 'mixity_quadratic':
            r['legend'] = legend(pg, 'cod/' + name, [('[0.04, 0.20]', 'blue'), ('[0.04, 0.30]', 'orange'), ('[0.08, 0.30]', 'green'), ('[0.12, 0.30]', 'purple')])
        return r
    run('cod/' + name, fn)
A7 = EXP_TABLES['fig07_metrics']['a_mm']
run('tip_resolution_sensitivity/trajectory', lambda pg, T: {'limits': lim_check('tipres/traj', T, (135, 300, -55, 15)),
    'h0': series(pg, T, 'tipres/traj:h0', 'blue', acc, dashed=False), 'h2': series(pg, T, 'tipres/traj:2h0', 'orange', pts([v[0] for v in Vm], [v[1] for v in Vm]), dashed=True),
    'legend': legend(pg, 'tipres/traj', [('Закон M0', None), ('ℎ0', 'blue'), ('2ℎ0', 'orange'), ('max', 'purple'), ('Локальна симетрія', 'green'), ('останній прийнятий', 'blue')])})
run('tip_resolution_sensitivity/vertical_deviation', lambda pg, T: {'s': series(pg, T, 'tipres/dy', 'blue', pts(A7, EXP_TABLES['fig07_metrics']['dy_um']))})
run('tip_resolution_sensitivity/direction_deviation', lambda pg, T: {'s': series(pg, T, 'tipres/dth', 'purple', pts(A7, EXP_TABLES['fig07_metrics']['dtheta_mdeg']))})
run('tip_resolution_sensitivity/mode_mixity', lambda pg, T: {'h0': series(pg, T, 'tipres/q:h0', 'blue', pts(A7, EXP_TABLES['fig07_metrics']['q_h0_1e3']), dashed=False),
    '2h0': series(pg, T, 'tipres/q:2h0', 'orange', pts(A7, EXP_TABLES['fig07_metrics']['q_2h0_1e3']), dashed=True), 'aLS': lsline(pg, T, 'tipres/q'),
    'legend': legend(pg, 'tipres/q', [('ℎ', 'blue'), ('2ℎ', 'orange')])})
A9 = EXP_TABLES['fig09_metrics']['a_mm']; N9 = list(zip(EXP_TABLES['fig09_native']['increment_mm'], EXP_TABLES['fig09_native']['a_mm'], EXP_TABLES['fig09_native']['turn_density_deg_per_mm']))
run('increment_sensitivity/vertical_deviation', lambda pg, T: {'2m4': series(pg, T, 'inc/dy:2-4', 'orange', pts(A9, EXP_TABLES['fig09_metrics']['dy_2m4_um'])),
    '1m2': series(pg, T, 'inc/dy:1-2', 'yellow', pts(A9, EXP_TABLES['fig09_metrics']['dy_1m2_um']))})
run('increment_sensitivity/direction_deviation', lambda pg, T: {'2m4': series(pg, T, 'inc/dth:2-4', 'orange', pts(A9, EXP_TABLES['fig09_metrics']['dth_2m4_mdeg'])),
    '1m2': series(pg, T, 'inc/dth:1-2', 'yellow', pts(A9, EXP_TABLES['fig09_metrics']['dth_1m2_mdeg']))})
run('increment_sensitivity/turn_density', lambda pg, T: {k: series(pg, T, f'inc/td:{k}', c, [(x, y) for i, x, y in N9 if i == inc])
    for k, c, inc in (('4mm', 'blue', 4), ('2mm', 'orange', 2), ('1mm', 'yellow', 1))})
A10 = EXP_TABLES['fig10_solver']['a_mm']; S10 = EXP_TABLES['fig10_mesh']['s_mm']
run('quality/min_angle', lambda pg, T: {'s': series(pg, T, 'q/angle', 'blue', pts(S10, EXP_TABLES['fig10_mesh']['min_angle_deg']))})
run('quality/size_ratio', lambda pg, T: {'s': series(pg, T, 'q/ratio', 'blue', pts(S10, EXP_TABLES['fig10_mesh']['max_neighbor_ratio']))})
run('quality/iterations', lambda pg, T: {'s': series(pg, T, 'q/pcg', 'blue', pts(A10, EXP_TABLES['fig10_solver']['pcg_iterations']))})
run('quality/true_residual', lambda pg, T: {'s': series(pg, T, 'q/res', 'blue', pts(A10, EXP_TABLES['fig10_solver']['residual_over_gate']))})

# ------------------------------------------------------------------ retained-layer panels
RL = json.load(open(os.path.join(a.retained, 'RETAINED_LAYERS.json'), encoding='utf-8'))['assets']
PROF = json.load(open(os.path.join(paper, 'data', 'retained_figure_7_9', 'profiles.json'), encoding='utf-8'))
PAINT = {b'S', b's', b'f', b'F', b'f*', b'B', b'B*', b'b', b'b*'}


def opsof(path):
    rd = PdfReader(path); return ContentStream(rd.pages[0].get_contents(), rd).operations


def norm(ops):
    return [([float(x) if hasattr(x, '__float__') and not isinstance(x, (bytes, str)) else repr(x) for x in o], op) for o, op in ops]


def xobject_ops(panel_pdf):
    rd = PdfReader(panel_pdf); res = rd.pages[0].get('/Resources', {}); xs = res.get('/XObject', {}) if res else {}
    found = []
    def walk(xd, depth=0):
        for k in xd:
            x = xd[k].get_object()
            if x.get('/Subtype') == '/Form':
                found.append(ContentStream(x, rd).operations)
                sub = x.get('/Resources', {}).get('/XObject', {}) if x.get('/Resources') else {}
                if sub and depth < 3:
                    walk(sub, depth + 1)
    walk(xs)
    return found


for rec in RL:
    gn = rec['panel']
    if only and gn not in only:
        continue
    r = {'source': rec['source'], 'sourceSHA256': sha(os.path.join(root, rec['source'])), 'layerSHA256': sha(os.path.join(root, rec['layer']))}
    if r['sourceSHA256'] != rec['sourceSHA256']:
        bad(gn, 'retained source changed')
    src = norm(opsof(os.path.join(root, rec['source']))); lay = norm(opsof(os.path.join(root, rec['layer'])))
    if rec['mode'] == 'span':
        a0, b0 = PROF[rec['name']]['graph']
        if PROF[rec['name']]['sha256'] != r['sourceSHA256']:
            bad(gn, 'source differs from profiles.json')
        # operators of the span appear unchanged, in order, in the layer; outside the span no paint other than n and no text
        idx = [i for i in range(len(lay)) if lay[i] == src[a0]]
        ok = any(lay[i:i + (b0 - a0)] == src[a0:b0] for i in idx)
        outside = [op for i, (o, op) in enumerate(lay) if op in PAINT] ; span_paints = sum(1 for o, op in src[a0:b0] if op in PAINT)
        r['spanUnchanged'] = ok; r['layerPaints'] = len(outside); r['spanPaints'] = span_paints
        r['noText'] = not any(op in (b'BT', b'Tj', b'TJ') for o, op in lay)
        if not ok or len(outside) != span_paints or not r['noText']:
            bad(gn, 'layer operators are not exactly the pinned graph span')
    elif rec['mode'] == 'mesh':
        # every mesh/crack stroke kept with identical operands; neutralized paints are only gray-0.149 decoration strokes,
        # gray-0.149 glyph fills and white backgrounds (colour state re-derived here, independently of the builder)
        diff = [i for i, (s, l) in enumerate(zip(src, lay)) if s != l]; G = [None]; g = [None]; why = []
        stG, stg = [], []; dset = set(diff)
        for i, (o, op) in enumerate(src):
            if op == b'q':
                stG.append(G[0]); stg.append(g[0])
            elif op == b'Q':
                G[0] = stG.pop(); g[0] = stg.pop()
            elif op == b'G':
                G[0] = o[0]
            elif op == b'g':
                g[0] = o[0]
            if i in dset:
                why.append('deco' if op == b'S' and G[0] == 0.149 else 'glyph' if op == b'f' and g[0] == 0.149 else 'bg' if op == b'f' and g[0] == 1.0 else 'OTHER')
        kept_src = sum(1 for o, op in src if op == b'S') - why.count('deco'); kept_lay = sum(1 for o, op in lay if op == b'S')
        r['lengthEqual'] = len(src) == len(lay); r['changedOperators'] = len(diff); r['neutralized'] = {k: why.count(k) for k in ('deco', 'glyph', 'bg', 'OTHER')}
        r['allChangesAreNeutralizedPaints'] = all(lay[i][1] == b'n' for i in diff) and why.count('OTHER') == 0
        r['meshAndCrackStrokesKept'] = kept_lay; r['meshAndCrackStrokesSource'] = kept_src
        if not (r['lengthEqual'] and r['allChangesAreNeutralizedPaints'] and kept_src == kept_lay and r['neutralized'] == {'deco': 18, 'glyph': 23, 'bg': 3, 'OTHER': 0}):
            bad(gn, 'mesh layer differs from the source beyond the neutralized decoration, glyph and background paints')
    else:
        # geometry kept: every source paint that is not an outlined glyph is still a paint in the layer, same operands
        diff = [i for i, (s, l) in enumerate(zip(src, lay)) if s != l]
        r['changedOperators'] = len(diff); r['changedAllGlyphFillsToN'] = all(lay[i][1] == b'n' and src[i][1] in (b'f', b'f*', b'F') for i in diff)
        r['lengthEqual'] = len(src) == len(lay)
        if not (r['lengthEqual'] and r['changedAllGlyphFillsToN'] and r['changedOperators'] == rec['glyphPathsRemoved']):
            bad(gn, 'geometry layer differs from the source beyond the removed glyph fills')
    pdf = os.path.join(a.figures, gn + '.pdf')
    xo = xobject_ops(pdf); r['layerEmbeddedUnchanged'] = any(norm(x) == lay for x in xo)
    if not r['layerEmbeddedUnchanged']:
        bad(gn, 'the placed Form XObject does not carry the layer operators unchanged')
    pg = pymupdf.open(pdf)[0]; r['fonts'] = fonts_ok(pg, gn)
    if rec['mode'] == 'span':
        T = transform(pg); r['transform'] = {'limits': T['lim'], 'ticks': T['ticks']}
        if rec['name'] in ('m1_trajectory', 'increment_trajectory'):
            circ = [d for d in pg.get_drawings() if len(d['items']) > 100 and abs(d['rect'].width - d['rect'].height) < 1 and d['rect'].width > 80]
            c = circ[0]['rect']; cx, cy = T['fx']((c.x0 + c.x1) / 2), T['fy']((c.y0 + c.y1) / 2); rr = (T['fx'](c.x1) - T['fx'](c.x0)) / 2
            r['hole'] = {'centre': [cx, cy], 'radius': rr, 'pass': abs(cx - 170) < 0.05 and abs(cy + 20) < 0.05 and abs(rr - 30) < 0.05}
            if not r['hole']['pass']:
                bad(gn, f'placed hole contour {cx:.4f},{cy:.4f},R {rr:.4f} differs from (170,-20,30)')
            if rec['name'] == 'm1_trajectory':   # the M0 curve is the reference history (evidence.json): independent check
                L = [p for p in polylines(pg, 'blue', dashed=False) if len(p) == 24]
                r['M0_vs_reference'] = match(gn + ':M0', L[0], acc, T, tolfrac=5e-4) if L else {'pass': False}
                if not L:
                    bad(gn, 'M0 polyline (24 vertices) not found')
                r['legend'] = legend(pg, gn, [('M0', 'blue'), ('M1', 'orange'), ('для M0: ', 'purple'), ('Локальна симетрія', 'green'), ('останній прийнятий', 'blue')])
            else:
                r['legend'] = legend(pg, gn, [('Масштабовані', None), ('= 4', 'blue'), ('= 2', 'orange'), ('= 1', 'yellow')])
                # the saved plot record keeps the 4/2/1-mm tip coordinates at the 23 common crack lengths (native vertices of
                # each history): every recorded point must coincide with a vertex of the retained polyline of that colour
                xmin, xmax, ymin, ymax = T['lim']; tx = 5e-4 * (xmax - xmin); ty = 5e-4 * (ymax - ymin); r['recordVertices'] = {}
                for k, colr in (('4mm', 'blue'), ('2mm', 'orange'), ('1mm', 'yellow')):
                    L = max(polylines(pg, colr, dashed=False), key=len, default=[])
                    V_ = [(T['fx'](p.x), T['fy'](p.y)) for p in L]; recp = list(zip(incm(f'x_{k}_mm'), incm(f'y_{k}_mm')))
                    dmax = max(min(max(abs(vx - ex) / tx, abs(vy - ey) / ty) for vx, vy in V_) for ex, ey in recp) if V_ else float('inf')
                    r['recordVertices'][k] = {'polylineVertices': len(L), 'recordPoints': len(recp), 'maxDeviationOverTolerance': dmax, 'pass': dmax <= 1}
                    if dmax > 1:
                        bad(gn, f'{k} retained trajectory does not pass through the recorded common-length vertices ({dmax:.2f} x tol)')
        elif rec['name'] == 'm1_mode_mixity':
            L = [p for p in polylines(pg, 'blue', dashed=False) if len(p) == 23]
            r['M0_vs_reference'] = match(gn + ':M0', L[0], pts(A7, EXP_TABLES['fig07_metrics']['q_h0_1e3']), T, tolfrac=5e-4) if L else {'pass': False}
            if not L:
                bad(gn, 'M0 polyline (23 vertices) not found')
            r['aLS'] = lsline(pg, T, gn); r['legend'] = legend(pg, gn, [('M0', 'blue'), ('M1', 'orange')])
        else:   # manuscript-quoted maxima of the retained M1 histories (main_ua.tex: 0.06734 um, 0.07304 mdeg)
            colr, quoted = ('blue', 0.06734) if rec['name'] == 'm1_vertical_deviation' else ('purple', 0.07304)
            L = max(polylines(pg, colr), key=len); mxv = max(T['fy'](p.y) for p in L)
            ok = abs(mxv - quoted) <= 5e-6 + 2e-4 * (T['lim'][3] - T['lim'][2])
            r['maximumVsManuscript'] = {'plotted': mxv, 'quoted': quoted, 'pass': ok, 'vertices': len(L)}
            if not ok:
                bad(gn, f'plotted maximum {mxv:.6f} differs from the manuscript value {quoted}')
    elif rec['mode'] == 'mesh':
        # independent calibration check: the crack (black, single straight stroke inside the frame) runs from x1 = -3 mm
        # (core radius r_core = 3 mm, main text) to the tip x1 = 0 along x2 = 0
        T = transform(pg); r['transform'] = {'limits': T['lim'], 'ticks': T['ticks']}; fr = T['frame']
        cr = [d for d in pg.get_drawings() if cname(d.get('color')) == 'black' and len(d['items']) == 1 and d['items'][0][0] == 'l'
              and fr.x0 + 2 < d['rect'].x0 and d['rect'].x1 < fr.x1 - 2 and fr.y0 + 2 < d['rect'].y0 < fr.y1 - 2 and d['rect'].width > 40]
        if len(cr) != 1:
            bad(gn, f'{len(cr)} crack strokes found'); r['crack'] = {'pass': False}
        else:
            p0, p1 = cr[0]['items'][0][1], cr[0]['items'][0][2]; xs = sorted([T['fx'](p0.x), T['fx'](p1.x)]); ys = [T['fy'](p0.y), T['fy'](p1.y)]
            ok = abs(xs[0] + 3) < 2e-3 and abs(xs[1]) < 2e-3 and all(abs(y) < 2e-3 for y in ys)
            r['crack'] = {'x1_mm': xs, 'x2_mm': ys, 'pass': ok}
            if not ok:
                bad(gn, f'crack line at x1 {xs}, x2 {ys} instead of [-3, 0] x 0')
        txt = ' '.join(texts(pg)); r['labels'] = {'мм': txt.count('мм') == 2, 'noLatinMm': 'mm' not in txt}
        if not all(r['labels'].values()):
            bad(gn, 'mesh axis labels are not the Ukrainian "x, мм" labels')
    else:
        T = None; txt = ' '.join(texts(pg)); need = ['𝜎', '𝐵', '2𝐴', '𝑅', '𝑂', 'Вільний', 'край']
        r['labels'] = {t: (t in txt) for t in need}
        if not all(r['labels'].values()):
            bad(gn, 'a schematic label is missing')
    P[gn] = r
# ------------------------------------------------------------------ approved English assets untouched
EN = {}
for rec in RL:
    p = os.path.join(root, rec['source']); EN[rec['source']] = sha(p) == rec['sourceSHA256']
RESULT['retainedSourcesUnchanged'] = EN
if not all(EN.values()):
    bad('sources', 'an approved English source asset changed')
RESULT['panels'] = P
RESULT['pass'] = not RESULT['problems']
txt = json.dumps(RESULT, indent=1, ensure_ascii=False, default=lambda o: o if not callable(o) else None)
if a.out:
    open(a.out, 'w', encoding='utf-8', newline='\n').write(txt + '\n')
print('UA7_VERIFY', 'PASS' if RESULT['pass'] else 'FAIL', 'panels', len(P), 'tables', len(RESULT['tables']), 'problems', len(RESULT['problems']))
for p in RESULT['problems'][:40]:
    print('  -', p)
sys.exit(0 if RESULT['pass'] else 1)
