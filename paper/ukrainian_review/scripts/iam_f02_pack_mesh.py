"""IAM-UA8-F02: pack the regenerated P21 mesh arrays (raw little-endian float64 exports written by
iam_f02_regenerate_mesh.m) into figures_src/iam_f02/data/*.f64.gz with a SHA-256 index. Deterministic gzip
(mtime 0). Only the arrays needed for the figure and its verification are packed.
usage: python -I iam_f02_pack_mesh.py <rawExportDir> [<outDir>]"""
import gzip, hashlib, json, os, sys
here = os.path.dirname(os.path.abspath(__file__)); ua = os.path.dirname(here)
src = sys.argv[1]; out = sys.argv[2] if len(sys.argv) > 2 else os.path.join(ua, 'figures_src', 'iam_f02', 'data')
os.makedirs(out, exist_ok=True)
IX = {x['file'][:-4]: x for x in json.load(open(os.path.join(src, 'RAW_INDEX.json'), encoding='utf-8'))}
KEEP = ['X3_mm', 'XL3_mm', 'T3', 'upperNodes', 'lowerNodes', 'tipNode', 'pathGlobal_mm', 'pairedElementIDs', 'exteriorElementIDs',
        'EDIsupportIDs', 'pairedNodeIDs', 'pairedMirrorLocal']
index = []
for k in KEEP:
    b = open(os.path.join(src, IX[k]['file']), 'rb').read()
    shape = IX[k]['shape'] if isinstance(IX[k]['shape'], list) else [IX[k]['shape']]
    f = os.path.join(out, k + '.f64.gz')
    assert not os.path.exists(f), f
    with open(f, 'wb') as fh:
        with gzip.GzipFile(filename='', mode='wb', fileobj=fh, mtime=0, compresslevel=9) as gz:
            gz.write(b)
    index.append({'name': k, 'file': k + '.f64.gz', 'shape': shape, 'order': 'column-major (MATLAB)', 'rawSHA256': hashlib.sha256(b).hexdigest(),
                  'gzSHA256': hashlib.sha256(open(f, 'rb').read()).hexdigest()})
json.dump({'schema': 'IAM-UA8-F02.MeshData.v1', 'state': 'P21 (21 segments, a = 84 mm) of the accepted M0 trajectory',
           'origin': 'deterministic regeneration with the pinned historical generator cfdf0f110010d688a4e3c48f6d88a00fd17dc698 '
                     '(qualify_incremental_crack_candidate, RunSynthetic=false; no FEM, no SIF/EDI evaluation); see audit/IAM_F02_MESH_FIGURE_REPORT.md',
           'units': 'coordinates in mm; indices 1-based (MATLAB)', 'arrays': index},
          open(os.path.join(out, 'DATA_INDEX.json'), 'w', encoding='utf-8', newline='\n'), indent=1)
print('packed', len(index), 'arrays into', out)
