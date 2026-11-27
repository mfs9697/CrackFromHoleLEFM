"""Retained vector layers for the panels WITHOUT preserved numerical histories (AAM-UA7).

Fig. 8(a-d) (M1 exterior sensitivity) and Fig. 9(a) (native 4/2/1-mm trajectories) have no surviving numerical history;
the repository keeps only the accepted vector PDFs (paper/data/retained_figure_7_9, pinned in profiles.json).
Digitization is prohibited, so these curves are NOT reconstructed: the original drawing operators of the pinned
"graph" span are copied unchanged into a layer PDF; every operator outside the span keeps the graphics state but its
paint operators are replaced by 'n' (end path, no paint) and every text block (BT..ET) is removed. All language-
dependent text (axes, ticks, labels, legends, annotations) is typeset anew in STIX Two by the Ukrainian panel.

Fig. 1 (approved symbolic specimen schematic, paper/figures/geometry_loading/specimen_geometry.pdf): every geometry
path is kept unchanged; only the outlined label glyphs (filled black paths with >= 12 construction operators) are
replaced by 'n'; the labels are typeset anew in STIX Two at the measured glyph positions.

Fig. 2 (structured crack-tip meshes 4h0/2h0/h0, paper/figures/mesh_levels/tip_mesh_H2/H1/H0.pdf): mesh coordinates and
connectivity were not saved and remeshing is prohibited, so every mesh and crack stroke is kept unchanged; only the MATLAB
axes frame/ticks (stroke gray 0.149), the outlined Computer Modern tick/axis-label glyphs (fill gray 0.149, English
"[mm]") and the white backgrounds are replaced by 'n'; axes, ticks and "x, мм" labels are typeset anew in STIX Two.

Text-state operators outside BT..ET (Tf etc.) and the then unused /Font resources are removed, so no source font enters
the Ukrainian panel. Content streams are written with exact (shortest round-trip, positional) real operands.

Calibration (pt -> data): lower panels from the source tick marks/labels (linear fit, residual reported); spatial
panels from the retained hole contour (centre (170,-20) mm, R = 30 mm, uniform scale) and cross-checked with ticks;
mesh panels from the axes frame = generator limits, per axis, cross-checked by tick marks at integer values.
Outputs: figures_src/retained/<name>_layer.pdf, <name>_calib.tex (TeX macros), RETAINED_LAYERS.json (record).
usage: python build_retained_layers.py [outDir]"""
import hashlib, io, json, os, re, sys
from decimal import Decimal
from pypdf import PdfReader, PdfWriter
from pypdf.generic import ArrayObject, ContentStream, DecodedStreamObject, FloatObject, NameObject, NumberObject
import pymupdf

here = os.path.dirname(os.path.abspath(__file__)); ua = os.path.dirname(here); paper = os.path.dirname(ua); root = os.path.dirname(paper)
out = sys.argv[1] if len(sys.argv) > 1 else os.path.join(ua, 'figures_src', 'retained')
os.makedirs(out, exist_ok=True)
RET = os.path.join(paper, 'data', 'retained_figure_7_9'); prof = json.load(open(os.path.join(RET, 'profiles.json'), encoding='utf-8'))
PAINT = {b'S', b's', b'f', b'F', b'f*', b'B', b'B*', b'b', b'b*'}
CONSTR = {b'm', b'l', b'c', b'v', b'y', b'h', b're'}
sha = lambda b: hashlib.sha256(b).hexdigest()
rel = lambda p: os.path.relpath(p, root).replace(os.sep, '/')
ASSETS = [
    # name, source, mode, target panel, axis limits (x,y) of the new axis in the units of the new axis label, y scale
    ('m1_trajectory', os.path.join(RET, 'm1_trajectory.pdf'), 'span', 'm1_mesh_sensitivity/trajectory', (135, 300, -55, 15), 'circle', 1.0),
    ('m1_vertical_deviation', os.path.join(RET, 'm1_vertical_deviation.pdf'), 'span', 'm1_mesh_sensitivity/vertical_deviation', None, 'ticks', 1.0),
    ('m1_direction_deviation', os.path.join(RET, 'm1_direction_deviation.pdf'), 'span', 'm1_mesh_sensitivity/direction_deviation', None, 'ticks', 1.0),
    ('m1_mode_mixity', os.path.join(RET, 'm1_mode_mixity.pdf'), 'span', 'm1_mesh_sensitivity/mode_mixity', None, 'ticks', 1.0),
    ('increment_trajectory', os.path.join(RET, 'increment_trajectory.pdf'), 'span', 'increment_sensitivity/trajectory', (135, 300, -55, 15), 'circle', 1.0),
    ('specimen_geometry', os.path.join(paper, 'figures', 'geometry_loading', 'specimen_geometry.pdf'), 'glyphs', 'geometry_loading/specimen_geometry', None, 'page', 1.0),
] + [('tip_mesh_' + h, os.path.join(paper, 'figures', 'mesh_levels', 'tip_mesh_%s.pdf' % h), 'mesh', 'mesh_levels/tip_mesh_' + h,
      (-3.15, 0.85, -1.65, 1.65), 'frame', 1.0) for h in ('H2', 'H1', 'H0')]


def num(s):
    s = s.strip().replace('−', '-').replace('−', '-')
    return float(s) if re.fullmatch(r'-?\d+(\.\d+)?', s) else None


def segments(ops):
    """index ranges [i0, ipaint] of path segments (construction ops followed by one paint/n op)"""
    seg = []; start = None
    for i, (_, op) in enumerate(ops):
        if op in CONSTR:
            if start is None:
                start = i
        elif op in PAINT or op == b'n':
            if start is not None:
                seg.append((start, i)); start = None
    return seg


def filter_ops(ops, keep_span=None, glyph_rule=False):
    new = []; dropped_text = 0; neutralized = 0; glyphs = 0; intext = False; tf = 0
    fill_black = True  # track nonstroking colour for the glyph rule
    seg_end = {e: s for s, e in segments(ops)}
    for i, (opd, op) in enumerate(ops):
        if op == b'g':
            fill_black = float(opd[0]) == 0.0
        elif op == b'rg':
            fill_black = all(float(v) == 0.0 for v in opd)
        inside = keep_span is not None and keep_span[0] <= i < keep_span[1]
        if op == b'BT':
            intext = True
        if intext:
            dropped_text += 1
            if op == b'ET':
                intext = False
            continue
        if op in (b'Tf', b'Tc', b'Tw', b'Tz', b'TL', b'Tr', b'Ts'):   # text-state operators outside BT..ET (no painting effect)
            assert not inside, ('text-state operator inside the pinned graph span', i)
            tf += 1; continue
        if keep_span is not None and not inside and op in PAINT:
            new.append(([], b'n')); neutralized += 1; continue
        if glyph_rule and op in (b'f', b'f*', b'F') and i in seg_end and fill_black and (i - seg_end[i]) >= 12:
            new.append(([], b'n')); glyphs += 1; continue
        new.append((opd, op))
    return new, dict(textOpsRemoved=dropped_text, textStateOpsRemoved=tf, paintsNeutralizedOutsideSpan=neutralized, glyphPathsRemoved=glyphs)


def filter_mesh(ops):
    """Fig. 2 mesh panels: keep every mesh (gray) and crack (black) stroke unchanged; replace by 'n' only the paints of the
    MATLAB axes decoration (frame and ticks: stroke gray 0.149), of the outlined tick/axis-label glyphs (fill gray 0.149) and of
    the white page/axes backgrounds (fill 1). Colour state is tracked through q/Q."""
    new = []; st = [{'G': None, 'g': None}]; n = {'decorationStrokes': 0, 'glyphFills': 0, 'backgroundFills': 0}
    for opd, op in ops:
        s = st[-1]
        if op == b'q':
            st.append(dict(s))
        elif op == b'Q':
            st.pop()
        elif op == b'G':
            s['G'] = float(opd[0])
        elif op == b'g':
            s['g'] = float(opd[0])
        elif op in (b'RG', b'rg', b'K', b'k', b'CS', b'cs', b'SC', b'sc', b'SCN', b'scn'):
            raise AssertionError('unexpected colour operator in mesh panel')
        if op in (b'S', b's') and s['G'] is not None and abs(s['G'] - 0.149) < 1e-6:
            new.append(([], b'n')); n['decorationStrokes'] += 1; continue
        if op in (b'f', b'f*', b'F'):
            if abs(s['g'] - 0.149) < 1e-6:
                new.append(([], b'n')); n['glyphFills'] += 1; continue
            if s['g'] == 1.0:
                new.append(([], b'n')); n['backgroundFills'] += 1; continue
            raise AssertionError('unexpected fill in mesh panel')
        assert op not in (b'B', b'B*', b'b', b'b*', b'BT'), 'unexpected paint/text operator in mesh panel'
        new.append((opd, op))
    return new, n


def operand(x):
    if isinstance(x, FloatObject):
        r = format(Decimal(repr(float(x))), 'f')   # shortest round-trip digits, positional (PDF has no exponent form)
        if '.' in r:
            r = r.rstrip('0').rstrip('.')
        assert float(r) == float(x) and 'e' not in r.lower()
        return r.encode()
    if isinstance(x, NumberObject):
        return str(int(x)).encode()
    if isinstance(x, ArrayObject):
        return b'[' + b' '.join(operand(v) for v in x) + b']'
    bio = io.BytesIO(); x.write_to_stream(bio); return bio.getvalue()


def serialize(ops):
    out = []
    for opd, op in ops:
        assert op not in (b'INLINE IMAGE', b'BI', b'ID', b'EI'), 'inline images are not expected'
        out.append(b' '.join([operand(v) for v in opd] + [op]))
    return b'\n'.join(out) + b'\n'


def ticks_fit(pg, axes_td):
    """linear pt->value maps from tick labels (centre of each numeric label) for x (below axes) and y (left of axes)"""
    x0, y0, x1, y1 = axes_td; X = []; Y = []
    for b in pg.get_text('dict')['blocks']:
        for l in b.get('lines', []):
            for s in l['spans']:
                v = num(s['text']); bx = s['bbox']
                if v is None or abs(s['size'] - 8.5) > 0.3:   # tick labels only (8.5 pt); axis-label digits differ in size
                    continue
                cx = (bx[0] + bx[2]) / 2; cy = (bx[1] + bx[3]) / 2
                if 0 <= bx[1] - y1 <= 10 and x0 - 8 <= cx <= x1 + 8:          # directly below the axes
                    X.append((cx, v))
                elif 0 <= x0 - bx[2] <= 10 and y0 - 8 <= cy <= y1 + 8:        # right-aligned against the left axis
                    Y.append((cy, v))
    # tick marks: short strokes on the axes boundary give the exact tick positions; snap each label to the nearest mark
    marks_x = []; marks_y = []
    for d in pg.get_drawings():
        for it in d['items']:
            if it[0] == 'l':
                p, q = it[1], it[2]
                if abs(p.x - q.x) < 1e-6 and abs(p.y - q.y) < 6 and (abs(p.y - y1) < 0.6 or abs(q.y - y1) < 0.6):
                    marks_x.append(p.x)
                if abs(p.y - q.y) < 1e-6 and abs(p.x - q.x) < 6 and (abs(p.x - x0) < 0.6 or abs(q.x - x0) < 0.6):
                    marks_y.append(p.y)

    def fit(pairs, marks):
        snapped = []
        for c, v in pairs:
            m = min(marks, key=lambda t: abs(t - c)) if marks else c
            snapped.append((m if abs(m - c) < 2.0 else c, v))
        n = len(snapped); sx = sum(p for p, _ in snapped); sv = sum(v for _, v in snapped)
        sxx = sum(p * p for p, _ in snapped); sxv = sum(p * v for p, v in snapped)
        b = (n * sxv - sx * sv) / (n * sxx - sx * sx); a = (sv - b * sx) / n
        res = max(abs(a + b * p - v) for p, v in snapped)
        return a, b, res, snapped
    return fit(X, marks_x), fit(Y, marks_y)


def circle_fit(pg, bbox_td):
    """hole contour: largest closed polyline whose bounding box matches the pinned circle box"""
    best = None
    for d in pg.get_drawings():
        r = d['rect']
        if len(d['items']) > 100 and abs(r.width - r.height) < 0.5 and r.width > 50:
            if best is None or abs(r.x0 - bbox_td[0]) + abs(r.y0 - bbox_td[1]) < abs(best.x0 - bbox_td[0]) + abs(best.y0 - bbox_td[1]):
                best = r
    return best


records = []
for name, src, mode, panel, lims, calib, _ in ASSETS:
    raw = open(src, 'rb').read(); rd = PdfReader(src); page = rd.pages[0]
    ops = ContentStream(page.get_contents(), rd).operations
    rec = {'name': name, 'source': rel(src), 'sourceSHA256': sha(raw), 'panel': panel, 'mode': mode, 'operators': len(ops)}
    if mode == 'span':
        p = prof[name]; assert sha(raw) == p['sha256'], (name, 'source differs from profiles.json')
        a, b = p['graph']; new, info = filter_ops(ops, keep_span=(a, b)); rec.update(span=[a, b], **info)
        assert info['textOpsRemoved'] > 0 and info['paintsNeutralizedOutsideSpan'] > 0
    elif mode == 'mesh':
        new, info = filter_mesh(ops); rec.update(**info)
        assert info['glyphFills'] == 23 and info['decorationStrokes'] == 18 and info['backgroundFills'] == 3, info
    else:
        new, info = filter_ops(ops, glyph_rule=True); rec.update(**info)
    w = PdfWriter(); w.add_page(page); np_ = w.pages[0]
    # exact serialization: pypdf's own writer shortens real operands (e.g. -60.41794205 -> -60.4179421), which would alter
    # the retained drawing; each real is written with the shortest round-trip representation (identical IEEE double).
    sdat = DecodedStreamObject(); sdat.set_data(serialize(new))
    np_[NameObject('/Contents')] = w._add_object(sdat.flate_encode())   # lossless Flate of the exact bytes
    if '/Annots' in np_:
        del np_['/Annots']
    # with every text block removed the layer must not carry the source fonts (Times/CM subsets) into the panel PDF
    assert not any(op in (b'BT', b'Tf', b'Tj', b'TJ', b"'", b'"') for _, op in new)
    res = np_['/Resources'].get_object()
    if '/Font' in res:
        del res['/Font']; rec['fontResourcesRemoved'] = True
    lf = os.path.join(out, name + '_layer.pdf')
    with open(lf, 'wb') as fh:
        w.write(fh)
    rec['layer'] = rel(lf); rec['layerSHA256'] = sha(open(lf, 'rb').read())
    # kept operators equal the source operators (values), checked again independently by the verifier
    kept = [(o, op) for (o, op) in new if op != b'n']
    rec['keptPaintOperators'] = sum(1 for _, op in kept if op in PAINT)
    pg = pymupdf.open(src)[0]; W, H = pg.rect.width, pg.rect.height; rec['page_pt'] = [W, H]
    if calib == 'circle':
        r = circle_fit(pg, p['circle_bbox_topdown_pt'])
        s = r.width / 60.0; cx = (r.x0 + r.x1) / 2; cy = (r.y0 + r.y1) / 2
        xL = 170 - cx / s; xR = xL + W / s; yT = -20 + cy / s; yB = yT - H / s
        rec['calibration'] = {'method': 'hole contour (centre 170,-20 mm; R 30 mm)', 'circle_bbox_pt': [r.x0, r.y0, r.x1, r.y1],
                              'pt_per_mm': s, 'aspect_check_pt': r.height / 60.0}
        rec['edges'] = {'xL': xL, 'xR': xR, 'yB': yB, 'yT': yT}; rec['axis'] = {'xmin': lims[0], 'xmax': lims[1], 'ymin': lims[2], 'ymax': lims[3]}
    elif calib == 'ticks':
        p = prof[name]; ax_td = p['axes_bbox_topdown_pt']
        (axf, bxf, rx, sx_), (ayf, byf, ry, sy_) = ticks_fit(pg, ax_td)
        scale = 1e3 if name == 'm1_mode_mixity' else 1.0   # source ticks show K_II/K_I in units of 1e-3 (x10^-3 annotation)
        xL = axf; xR = axf + bxf * W; yT = ayf; yB = ayf + byf * H
        xmin, xmax = axf + bxf * ax_td[0], axf + bxf * ax_td[2]; ymax, ymin = ayf + byf * ax_td[1], ayf + byf * ax_td[3]
        rec['calibration'] = {'method': 'tick labels snapped to tick marks, linear least squares', 'xTicks': sx_, 'yTicks': sy_,
                              'xResidual': rx, 'yResidual': ry, 'axes_bbox_pt': ax_td, 'yTickUnit': '1e-3' if scale == 1e3 else '1',
                              'profileYlim': p['ylim'], 'ylimFromTicks': [ymin / scale, ymax / scale]}
        rec['edges'] = {'xL': xL, 'xR': xR, 'yB': yB, 'yT': yT}; rec['axis'] = {'xmin': xmin, 'xmax': xmax, 'ymin': ymin, 'ymax': ymax}
    elif calib == 'frame':
        # axes frame (stroke gray 0.149, long lines) = generator limits XLim_mm [-3.15 0.85], YLim_mm [-1.65 1.65] with axis equal
        # (verification/crack_path/main_tip_core_mesh_level_comparison.m defaults); cross-checked with the integer tick marks.
        xmin, xmax, ymin, ymax = lims; fr = None; tx = []; ty = []
        for d in pg.get_drawings():
            if d.get('color') is None or abs(d['color'][0] - 0.149) > 2e-3 or d.get('fill') is not None:
                continue
            r = d['rect']
            if r.width > 50 or r.height > 50:   # frame sides are zero-area rects: union by coordinates (Rect | ignores empties)
                b4 = (r.x0, r.y0, r.x1, r.y1) if fr is None else (min(fr.x0, r.x0), min(fr.y0, r.y0), max(fr.x1, r.x1), max(fr.y1, r.y1))
                fr = pymupdf.Rect(*b4)
        for d in pg.get_drawings():
            if d.get('color') is None or abs(d['color'][0] - 0.149) > 2e-3:
                continue
            for it in d['items']:
                if it[0] == 'l' and abs(it[1] - it[2]) < 6:
                    p_, q_ = it[1], it[2]
                    if abs(p_.x - q_.x) < 1e-6 and min(abs(p_.y - fr.y1), abs(q_.y - fr.y1)) < 0.6:
                        tx.append(p_.x)
                    if abs(p_.y - q_.y) < 1e-6 and min(abs(p_.x - fr.x0), abs(q_.x - fr.x0)) < 0.6:
                        ty.append(p_.y)
        sx = fr.width / (xmax - xmin); sy = fr.height / (ymax - ymin)
        xt = sorted({round(xmin + (t - fr.x0) / sx, 6) + 0.0 for t in tx}); yt = sorted({round(ymax - (t - fr.y0) / sy, 6) + 0.0 for t in ty})
        res_x = max(abs(v - round(v)) for v in xt); res_y = max(abs(v - round(v)) for v in yt)
        # The MATLAB export is 0.25 % anisotropic (content cm 4.8017254 vs 4.81354911); each axis is calibrated separately,
        # so every vertex lands on its own coordinate in the exactly equal-scale Ukrainian axis (recorded, not hidden).
        assert abs(sx / sy - 1) < 5e-3 and res_x < 2e-3 and res_y < 2e-3, (sx, sy, xt, yt)
        assert xt == [-3.0, -2.0, -1.0, 0.0] and yt == [-1.0, 0.0, 1.0], (xt, yt)
        xL = xmin - fr.x0 / sx; xR = xL + W / sx; yT = ymax + fr.y0 / sy; yB = yT - H / sy
        rec['calibration'] = {'method': 'axes frame = generator limits (XLim/YLim defaults), tick marks at integers as cross-check',
                              'frame_pt': [fr.x0, fr.y0, fr.x1, fr.y1], 'pt_per_mm_x': sx, 'pt_per_mm_y': sy,
                              'sourceAnisotropy': sx / sy - 1,
                              'tickValuesFromMarks_x': xt, 'tickValuesFromMarks_y': yt, 'tickResidual': [res_x, res_y]}
        rec['edges'] = {'xL': xL, 'xR': xR, 'yB': yB, 'yT': yT}; rec['axis'] = {'xmin': xmin, 'xmax': xmax, 'ymin': ymin, 'ymax': ymax}
    else:
        rec['edges'] = {'W_pt': W, 'H_pt': H}
    # TeX macros
    L = ['% generated by build_retained_layers.py -- do not edit', '\\def\\uaRetainedLayer{../../retained/%s_layer.pdf}' % name]
    for k, v in rec.get('edges', {}).items():
        L.append('\\def\\uaR%s{%r}' % (k.replace('_', ''), float(v)))
    for k, v in rec.get('axis', {}).items():
        L.append('\\def\\uaR%s{%r}' % (k, float(v)))
    open(os.path.join(out, name + '_calib.tex'), 'w', encoding='utf-8', newline='\n').write('\n'.join(L) + '\n')
    records.append(rec)
    print(name, {k: rec[k] for k in ('textOpsRemoved', 'paintsNeutralizedOutsideSpan', 'glyphPathsRemoved', 'keptPaintOperators') if k in rec},
          rec.get('axis'), rec.get('calibration', {}).get('xResidual'), rec.get('calibration', {}).get('yResidual'))
json.dump({'schema': 'AAM-UA7.RetainedLayers.v1', 'assets': records}, open(os.path.join(out, 'RETAINED_LAYERS.json'), 'w', encoding='utf-8', newline='\n'), indent=1)
