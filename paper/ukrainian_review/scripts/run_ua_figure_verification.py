"""AAM-UA7 WP3 driver: baseline independent verification + negative controls -> audit/FIGURE_DATA_VERIFICATION.json.

Baseline: verify_ua_figure_data.py on the committed figures and data tables (must PASS).
Negative controls (each must be REJECTED): a copy of figures_src is mutated in exactly one way, only the affected panel is
rebuilt with build_ua_figures.ps1 -SourceRoot/-OutputRoot into a temporary directory, and the verifier is run
  (L2) with the ORIGINAL data tables against the mutated PDF -> the PDF-vs-record comparison must fail, and
  (L1) for table mutations, with the mutated data tables -> the table-vs-record comparison must fail.
The committed figures, data and records are never modified. Windows: requires XeLaTeX (MiKTeX) and PowerShell.
usage: python run_ua_figure_verification.py [--work DIR] [--out FILE]"""
import argparse, json, os, shutil, subprocess, sys, tempfile
sys.stdout.reconfigure(encoding='utf-8')
here = os.path.dirname(os.path.abspath(__file__)); ua = os.path.dirname(here)
ap = argparse.ArgumentParser()
ap.add_argument('--work', default='')
ap.add_argument('--out', default=os.path.join(ua, 'audit', 'FIGURE_DATA_VERIFICATION.json'))
a = ap.parse_args()
work = a.work or tempfile.mkdtemp(prefix='ua7_negctl_')
os.makedirs(work, exist_ok=True); assert not os.listdir(work), 'work directory must be empty'
VER = os.path.join(here, 'verify_ua_figure_data.py'); BUILD = os.path.join(here, 'build_ua_figures.ps1')


def verify(tag, **kw):
    out = os.path.join(work, tag + '.json'); cmd = [sys.executable, VER, '--out', out]
    for k, v in kw.items():
        cmd += ['--' + k, v]
    p = subprocess.run(cmd, capture_output=True, text=True, encoding='utf-8')
    return p.returncode, json.load(open(out, encoding='utf-8'))


def csv_mutate(path, fn):
    L = open(path, encoding='utf-8').read().strip().split('\n'); h = L[0].split(','); rows = [l.split(',') for l in L[1:]]
    h, rows = fn(h, rows)
    open(path, 'w', encoding='utf-8', newline='\n').write('\n'.join([','.join(h)] + [','.join(r) for r in rows]) + '\n')


def col_map(name, f):
    def fn(h, rows):
        j = h.index(name)
        for i, r in enumerate(rows):
            r[j] = f(i, r[j])
        return h, rows
    return fn


def text_mutate(path, old, new):
    t = open(path, encoding='utf-8').read(); assert t.count(old) == 1, (path, old); open(path, 'w', encoding='utf-8', newline='\n').write(t.replace(old, new))


def layer_mutate(path):
    """change one real operand inside the retained graph span by +1 pt (a moved data vertex)"""
    from pypdf import PdfReader, PdfWriter
    from pypdf.generic import DecodedStreamObject, NameObject
    rd = PdfReader(path); data = rd.pages[0].get_contents().get_data()
    i = data.index(b' l\n', len(data) // 2); j = data.rindex(b' ', 0, i); old = data[j + 1:i]
    data = data[:j + 1] + repr(float(old) + 1.0).encode() + data[i:]
    w = PdfWriter(); w.add_page(rd.pages[0]); s = DecodedStreamObject(); s.set_data(data)
    w.pages[0][NameObject('/Contents')] = w._add_object(s); w.write(path)


def calib_shift(path, key, dv):
    import re
    t = open(path, encoding='utf-8').read(); m = re.search(r'\\def\\uaR%s\{([^}]*)\}' % key, t)
    v = float(m.group(1)) + dv; t = t[:m.start(1)] + repr(v) + t[m.end(1):]
    open(path, 'w', encoding='utf-8', newline='\n').write(t)


CONTROLS = [
    # id, kind, panel, mutation (applied to the copied figures_src), table mutated?
    ('C1_modified_value', 'modified numerical value: K_I of state 10 multiplied by 1.01', 'intensities/KI',
     lambda s: csv_mutate(os.path.join(s, 'data', 'fig04_states.csv'), col_map('KI', lambda i, v: repr(float(v) * 1.01) if i == 9 else v)), True),
    ('C2_reversed_sign', 'reversed sign: 1e3 K_II/K_I of the late states negated', 'late/mode_mixity',
     lambda s: csv_mutate(os.path.join(s, 'data', 'fig05_late.csv'), col_map('q_1e3', lambda i, v: repr(-float(v)))), True),
    ('C3_unit_conversion', 'incorrect unit conversion: 1e2 instead of 1e3 for the K_II axis (10^3 K_II)', 'intensities/KII',
     lambda s: csv_mutate(os.path.join(s, 'data', 'fig04_states.csv'), col_map('KII_1e3', lambda i, v: repr(float(v) / 10))), True),
    ('C4_missing_curve', 'missing curve: 2h0 series not plotted (legend entry kept)', 'tip_resolution_sensitivity/mode_mixity',
     lambda s: text_mutate(os.path.join(s, 'panels', 'tip_resolution_sensitivity', 'mode_mixity.tex'),
                           r'\addplot[ua alt] table[x=a_mm,y=q_2h0_1e3]{\T};', r'\addlegendimage{ua alt}'), False),
    ('C5_swapped_series', 'swapped series identity: y(2)-y(4) and y(1)-y(2) columns exchanged', 'increment_sensitivity/vertical_deviation',
     lambda s: csv_mutate(os.path.join(s, 'data', 'fig09_metrics.csv'),
                          lambda h, r: ([{'dy_2m4_um': 'dy_1m2_um', 'dy_1m2_um': 'dy_2m4_um'}.get(x, x) for x in h], r)), True),
    ('C6_unsolved_as_accepted', 'state status: unsolved P24 drawn as part of the accepted trajectory', 'trajectory/overview',
     lambda s: text_mutate(os.path.join(s, 'panels', 'trajectory', 'overview.tex'), r'{\thisrow{k}}{0:23}]{\V};', r'{\thisrow{k}}{0:24}]{\V};'), False),
    ('C7_small_coordinate_shift', 'modified numerical value (small): y of P20 scaled by 1.001 (0.025 mm)', 'trajectory/late_detail',
     lambda s: csv_mutate(os.path.join(s, 'data', 'fig03_vertices.csv'), col_map('y_mm', lambda i, v: repr(float(v) * (1.0 + 1e-3)) if i == 20 else v)), True),
    ('C8_retained_layer_vertex', 'retained vector layer altered: one drawing operand moved by 1 pt', 'm1_mesh_sensitivity/vertical_deviation',
     lambda s: layer_mutate(os.path.join(s, 'retained', 'm1_vertical_deviation_layer.pdf')), False),
    ('C9_mesh_calibration', 'retained mesh layer misplaced: left page edge of the calibration shifted by 0.05 mm', 'mesh_levels/tip_mesh_H1',
     lambda s: calib_shift(os.path.join(s, 'retained', 'tip_mesh_H1_calib.tex'), 'xL', 0.05), False),
]

res = {'schema': 'AAM-UA7.FigureDataVerification.v1', 'verifier': 'paper/ukrainian_review/scripts/verify_ua_figure_data.py'}
code, base = verify('baseline'); res['baseline'] = base; res['baselinePass'] = code == 0 and base['pass']
print('BASELINE', 'PASS' if res['baselinePass'] else 'FAIL', len(base['problems']))
src0 = os.path.join(ua, 'figures_src'); out = []
for cid, desc, panel, mut, table in CONTROLS:
    s = os.path.join(work, cid, 'figures_src'); o = os.path.join(work, cid, 'figures')
    shutil.copytree(src0, s, ignore=shutil.ignore_patterns('.build')); mut(s)
    b = subprocess.run(['powershell', '-NoProfile', '-ExecutionPolicy', 'Bypass', '-File', BUILD, '-Only', panel, '-SourceRoot', s, '-OutputRoot', o],
                       capture_output=True, text=True)
    r = {'id': cid, 'mutation': desc, 'panel': panel, 'built': b.returncode == 0}
    if b.returncode != 0:
        r['buildOutput'] = (b.stdout + b.stderr)[-800:]
    else:
        c2, v2 = verify(cid + '_L2', figures=o, only=panel)
        r['L2'] = {'rejected': c2 != 0 and not v2['pass'], 'problems': v2['problems']}
    if table:
        c1, v1 = verify(cid + '_L1', data=os.path.join(s, 'data'), only=panel)
        r['L1'] = {'rejected': c1 != 0 and not v1['pass'], 'problems': [p for p in v1['problems'] if not p.startswith(panel)][:6]}
    r['detected'] = r['built'] and r.get('L2', {}).get('rejected', False) and (not table or r['L1']['rejected'])
    out.append(r); print(cid, 'DETECTED' if r['detected'] else 'NOT DETECTED', r.get('L2', {}).get('problems', [''])[:1])
res['negativeControls'] = out; res['negativeControlsDetected'] = sum(r['detected'] for r in out); res['negativeControlsTotal'] = len(out)
res['pass'] = res['baselinePass'] and res['negativeControlsDetected'] == len(out)
open(a.out, 'w', encoding='utf-8', newline='\n').write(json.dumps(res, indent=1, ensure_ascii=False) + '\n')
print('UA7_FIGURE_VERIFICATION', 'PASS' if res['pass'] else 'FAIL', f"controls {res['negativeControlsDetected']}/{len(out)}", 'work', work)
sys.exit(0 if res['pass'] else 1)
