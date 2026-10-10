"""Current read-only portable manuscript audit; local archive status is explicit."""
from pathlib import Path
import argparse,json,runpy
paper=Path(__file__).resolve().parent
parser=argparse.ArgumentParser()
parser.add_argument('--require-local-archive',action='store_true',help='Fail if the original full reference archive is unavailable.')
args=parser.parse_args()
e=json.loads((paper/'data/evidence.json').read_text())
assert len(e['stateRows'])==23 and len(e['codFits'])==176
assert e['regression']['theta2_pass'] and e['regression']['step2_pass']
assert e['audit']['acceptedPhysicalCount']==22 and e['audit']['qualificationCount']==23
assert e['audit']['qualifiedUnsolvedSegment']==24
archive=paper.parent/'verification/crack_path/final_clean_run'
available=(archive/'path_run_state.mat').is_file()
if args.require_local_archive and not available:
    raise SystemExit('Original full reference archive unavailable; no local field revalidation performed.')
runpy.run_path(str(paper/'verify_final_presentation.py'))
print('Local reference archive:', 'available' if available else 'unavailable; committed evidence verified, raw fields not revalidated')
print('No reference manifests, numerical evidence, or solver sources were overwritten.')
