"""Audit paper data, source identity, and manuscript guardrails. No FE solve."""
from pathlib import Path
import csv
import hashlib
import json
import re
import subprocess

paper=Path(__file__).resolve().parent
root=paper.parent
run=root/'verification/crack_path/final_clean_run'
e=json.loads((paper/'data/evidence.json').read_text())
m=json.loads((paper/'data/derived_metrics.json').read_text())
text=(paper/'main.tex').read_text()
assert e['regression']['theta2_pass'] and e['regression']['step2_pass']
assert len(e['stateRows'])==23 and len(e['codFits'])==176
assert e['audit']['acceptedPhysicalCount']==22 and e['audit']['qualificationCount']==23
assert e['audit']['noPhysicalSolvePerformed']
assert all(x['segment']<=23 for x in e['stateRows'])
assert all(x['segment']<=23 for x in e['codFits'])
assert m['all_COD_signs_P21_positive'] and m['all_COD_signs_P22_negative']
assert m['source_evidence_sha256']==hashlib.sha256((paper/'data/evidence.json').read_bytes()).hexdigest()
assert '% Evidence JSON SHA256: '+m['source_evidence_sha256'] in text
assert text.endswith('\\end{document}\n')
assert not re.search(r'\\(?:input|include|includegraphics)\b',text)
assert not re.search(r'\\cite\w*\b',text)
abstract=text.split(r'\begin{abstract}')[1].split(r'\end{abstract}')[0]
assert 150<=len(abstract.split())<=220
assert text.count(r'\begin{figure}')==5
assert 'qualified but unsolved' in text and 'not an independently solved' in text

files=[('atomic_state',run/'path_run_state.mat')]
files += [('compact_physical',run/f'step_{k:03d}_physical_small.mat') for k in range(2,24)]
files += [('qualification',run/f'step_{k:03d}_qualification_small.mat') for k in range(2,25)]
files += [('configuration_source',run/'step_002_physical_solved.mat'),
          ('restored_stage1_source',paper/'data/accepted_stage1_source.mat')]
sources=['run_incremental_crack_path.m','qualify_incremental_crack_candidate.m',
         'solve_incremental_crack_tip.m','kink_angle_LEFM_MTS.m','SIF_LEFM_interaction_EDI.m',
         'native_COD_polyline_audit.m','verification/sif_audit/native_COD_audit.m',
         'T3toT6_fast.m','stif_assem.m','integr.m',
         'verification/crack_path/main_stage1_freeze_starting_state.m',
         'plot_final_clean_run_results.m','plot_final_clean_run_additional_results.m',
         'verification/crack_path/INCREMENTAL_PATH.md','verification/sif_audit/AUDIT_CLOSURE.md',
         'Documents/CrPathLEFM.tex']
files += [('scientific_source',root/name) for name in sources]
oldfile=root/'verification/crack_path/LOCAL_NUMERICAL_ARCHIVE_SHA256.csv'
old={}
if oldfile.exists():
    with oldfile.open(encoding='utf-8-sig',newline='') as f:
        old={x['RelativePath'].replace('\\','/'):x['SHA256'].lower() for x in csv.DictReader(f)}
rows=[];checked=0
for kind,p in files:
    assert p.is_file(),p
    rel=p.relative_to(root).as_posix()
    digest=hashlib.sha256(p.read_bytes()).hexdigest()
    if rel in old:
        assert old[rel]==digest, f'Stop: pre-existing archive fingerprint differs: {rel}'
        checked+=1
    rows.append({'role':kind,'repository_relative_path':rel,'bytes':p.stat().st_size,'sha256':digest})
with (paper/'data/source_manifest.csv').open('w',newline='',encoding='utf-8') as f:
    w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)
changed=subprocess.check_output(['git','diff','--name-only',e['sourceCommit'],'--',*sources],cwd=root,text=True)
assert not changed.strip(), f'Stop: existing scientific sources changed: {changed}'
report={'audit_date':'2026-10-07','source_commit':e['sourceCommit'],
        'source_manifest_count':len(rows),'prior_archive_hashes_matched':checked,
        'abstract_word_count':len(abstract.split()),'main_figures':4,'supplementary_figures':1,
        'accepted_states':23,'accepted_compact_physical_files':22,
        'qualification_files':23,'P24_physical_result_present':False,
        'existing_scientific_sources_changed':False,'physical_solve_performed':False,
        'bibliography_entries_added':0,'status':'PASS'}
(paper/'data/manuscript_checks.json').write_text(json.dumps(report,indent=2)+'\n')
print(json.dumps(report,indent=2))
