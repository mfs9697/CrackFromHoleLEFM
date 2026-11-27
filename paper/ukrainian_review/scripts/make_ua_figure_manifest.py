"""Write audit/FIGURE_SOURCE_MANIFEST.json: one record per graphic of the Ukrainian edition (AAM-UA7).

For each of the 31 PDFs included by main_ua.tex: figure/panel, method (reconstructed from records or retained vector
layer with text typeset anew), authoritative sources with SHA-256, plotted data tables, panel source, approved English
asset (pinned, unchanged), Ukrainian output SHA-256, printed slot and measured width, embedded fonts, printed text sizes,
and the result of the independent verification (audit/FIGURE_DATA_VERIFICATION.json).
usage: python make_ua_figure_manifest.py"""
import hashlib, json, os, re, subprocess
import pymupdf

here = os.path.dirname(os.path.abspath(__file__)); ua = os.path.dirname(here); paper = os.path.dirname(ua); root = os.path.dirname(paper)
sha = lambda p: hashlib.sha256(open(p, 'rb').read()).hexdigest()
rel = lambda p: os.path.relpath(p, root).replace(os.sep, '/')
mm = lambda pt: pt / 72 * 25.4
BASE = '7597cba1a378e6f436a3bfab304f85702902ee4e'
DM = json.load(open(os.path.join(ua, 'figures_src', 'data', 'DATA_MANIFEST.json'), encoding='utf-8'))
RL = {r['panel']: r for r in json.load(open(os.path.join(ua, 'figures_src', 'retained', 'RETAINED_LAYERS.json'), encoding='utf-8'))['assets']}
VER = json.load(open(os.path.join(ua, 'audit', 'FIGURE_DATA_VERIFICATION.json'), encoding='utf-8'))
TEX = open(os.path.join(ua, 'main_ua.tex'), encoding='utf-8').read()
FIG = {  # panel -> (figure, panel letter, printed slot mm)
    'geometry_loading/specimen_geometry': (1, '', 129.48),
    'mesh_levels/tip_mesh_H2': (2, 'a', 52.29), 'mesh_levels/tip_mesh_H1': (2, 'b', 52.29), 'mesh_levels/tip_mesh_H0': (2, 'c', 52.29),
    'trajectory/overview': (3, 'a', 124.832), 'trajectory/late_detail': (3, 'b', 124.832),
    'intensities/KI': (4, 'a', 52.29), 'intensities/KII': (4, 'b', 52.29), 'intensities/theta': (4, 'c', 52.29),
    'late/mode_mixity': (5, 'a', 64.408), 'late/turn': (5, 'b', 64.408),
    'cod/turn_linear': (6, 'a', 64.408), 'cod/turn_quadratic': (6, 'b', 64.408), 'cod/mixity_linear': (6, 'c', 64.408), 'cod/mixity_quadratic': (6, 'd', 64.408),
    'tip_resolution_sensitivity/trajectory': (7, 'a', 156.04), 'tip_resolution_sensitivity/vertical_deviation': (7, 'b', 52.29),
    'tip_resolution_sensitivity/direction_deviation': (7, 'c', 52.29), 'tip_resolution_sensitivity/mode_mixity': (7, 'd', 52.29),
    'm1_mesh_sensitivity/trajectory': (8, 'a', 156.04), 'm1_mesh_sensitivity/vertical_deviation': (8, 'b', 52.29),
    'm1_mesh_sensitivity/direction_deviation': (8, 'c', 52.29), 'm1_mesh_sensitivity/mode_mixity': (8, 'd', 52.29),
    'increment_sensitivity/trajectory': (9, 'a', 156.04), 'increment_sensitivity/vertical_deviation': (9, 'b', 52.29),
    'increment_sensitivity/direction_deviation': (9, 'c', 52.29), 'increment_sensitivity/turn_density': (9, 'd', 52.29),
    'quality/min_angle': (10, 'a', 64.408), 'quality/size_ratio': (10, 'b', 64.408), 'quality/iterations': (10, 'c', 64.408), 'quality/true_residual': (10, 'd', 64.408),
}
JUST = {
    'span': 'No numerical history survives (paper/data/retained_figure_7_9, pinned in profiles.json) and digitization is prohibited: the '
            'drawing operators of the pinned graph span are kept unchanged; axes, ticks, labels and legend are typeset anew in STIX Two.',
    'glyphs': 'Symbolic schematic without numerical parameters (illustrative crack): every geometry path of the approved composition is '
              'kept; the outlined label glyphs were removed and the labels typeset anew in STIX Two.',
    'mesh': 'Mesh coordinates/connectivity were not saved and remeshing is prohibited: every mesh and crack stroke is kept unchanged; '
            'the MATLAB frame/ticks, outlined Computer Modern labels ("[mm]") and backgrounds were neutralized; axes typeset anew in STIX Two.',
}
TABSRC = {'fig03': ['evidence'], 'fig04': ['evidence'], 'fig05': ['evidence'], 'fig06': ['evidence'], 'fig07': ['evidence', 'isolated'],
          'fig09': ['increment'], 'fig10': ['evidence']}
graphics = re.findall(r'\\includegraphics\[[^\]]*\]\{figures/([^}]+)\.pdf\}', TEX)
assert len(graphics) == 31 and set(graphics) == set(FIG), set(graphics) ^ set(FIG)
assets = []
for gn in graphics:
    f, p, slot = FIG[gn]; out = os.path.join(ua, 'figures', gn + '.pdf'); en = os.path.join(paper, 'figures', gn + '.pdf')
    en_pinned = subprocess.run(['git', 'rev-parse', f'{BASE}:paper/figures/{gn}.pdf'], cwd=root, capture_output=True, text=True).stdout.strip()
    en_blob = subprocess.run(['git', 'hash-object', en], cwd=root, capture_output=True, text=True).stdout.strip()
    src = os.path.join(ua, 'figures_src', 'panels', gn + '.tex'); t = open(src, encoding='utf-8').read()
    pg = pymupdf.open(out)[0]
    sizes = sorted({round(s['size'], 2) for b in pg.get_text('dict')['blocks'] for l in b.get('lines', []) for s in l['spans']})
    rec = {'path': f'paper/ukrainian_review/figures/{gn}.pdf', 'figure': f, 'panel': p, 'panelId': gn,
           'englishAsset': {'path': f'paper/figures/{gn}.pdf', 'sha256': sha(en), 'gitBlobEqualsPinnedBase': en_pinned == en_blob},
           'ukrainianSHA256': sha(out), 'panelSource': rel(src), 'panelSourceSHA256': sha(src),
           'printedSlot_mm': slot, 'measuredWidth_mm': round(mm(pg.rect.width), 3), 'height_mm': round(mm(pg.rect.height), 3),
           'includeScale': round(mm(pg.rect.width) / slot, 5), 'fonts': sorted({x[3].split('+')[-1] for x in pg.get_fonts()}),
           'textSizes_pt': sizes}
    tabs = re.findall(r'data/(\w+)\.csv', t)
    if gn in RL:
        r = RL[gn]
        rec['method'] = 'retained vector layer, text typeset anew' if r['mode'] != 'mesh' else 'retained mesh layer, axes typeset anew'
        rec['retainedLayer'] = {'source': r['source'], 'sourceSHA256': r['sourceSHA256'], 'layer': r['layer'], 'layerSHA256': r['layerSHA256'],
                                'mode': r['mode'], 'calibration': r.get('calibration', {}).get('method')}
        rec['justification'] = JUST[r['mode']]
    else:
        rec['method'] = 'reconstructed from authoritative records (PGFPlots, XeLaTeX)'
    rec['dataTables'] = [{'table': x, 'csv': f'paper/ukrainian_review/figures_src/data/{x}.csv',
                          'sha256': sha(os.path.join(ua, 'figures_src', 'data', x + '.csv')),
                          'sources': [DM['sources'][k]['path'] for k in TABSRC[x[:5]]],
                          'transformation': DM['tables'][x].get('transformation')} for x in dict.fromkeys(tabs)]
    v = VER['baseline']['panels'].get(gn, {}); rec['independentVerification'] = 'PASS' if VER['baseline']['pass'] and gn in VER['baseline']['panels'] else 'MISSING'
    assets.append(rec)
res = {'schema': 'AAM-UA7.FigureSourceManifest.v1', 'pinnedEnglishBase': BASE,
       'backend': 'PGFPlots 1.18 in standalone XeLaTeX documents; STIX Two Text and STIX Two Math (stix2-otf), the manuscript fonts',
       'extraction': 'scripts/extract_ua_figure_data.m (MATLAB, read-only on records) -> figures_src/data/*.csv (%.17g) + DATA_MANIFEST.json',
       'records': DM.get('sources', DM.get('records')),
       'counts': {'graphics': len(assets), 'reconstructedFromRecords': sum(a['method'].startswith('reconstructed') for a in assets),
                  'retainedVectorLayers': sum(a['method'].startswith('retained') for a in assets),
                  'unchangedEnglishFiles': 0},
       'assets': assets}
open(os.path.join(ua, 'audit', 'FIGURE_SOURCE_MANIFEST.json'), 'w', encoding='utf-8', newline='\n').write(json.dumps(res, indent=1, ensure_ascii=False) + '\n')
print('FIGURE_SOURCE_MANIFEST', res['counts'], 'all verified', all(a['independentVerification'] == 'PASS' for a in assets),
      'english pinned', all(a['englishAsset']['gitBlobEqualsPinnedBase'] for a in assets), 'scale', sorted({a['includeScale'] for a in assets}))
