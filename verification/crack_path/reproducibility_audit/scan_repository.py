from pathlib import Path
import argparse,csv,hashlib,json,re,subprocess
parser=argparse.ArgumentParser();parser.add_argument('--output-dir',type=Path,required=True)
args=parser.parse_args()
root=Path(__file__).resolve().parents[3]
out=args.output_dir;out.mkdir(parents=True,exist_ok=True)
suffixes={'.m','.py','.sh','.ps1','.bat','.cmd','.js','.mjs','.tex','.bib','.md','.json','.yml','.yaml','.toml'}
patterns={
 'increment_4mm':r'(?i)4\s*[- ]?\s*mm|(?:a0|increment|da)\w*\s*[=:]\s*(?:0?\.004|4\b)',
 'core_scale':r'CoreScale|coreScale|coreMeshControls|[\'\"]Scale[\'\"]\s*,\s*1\b',
 'fixed_states_lengths':r'\b(?:23|92)\b|\bP_?\{?(?:17|21|22|23)\}?\b',
 'reference_fingerprints':r'\b(?:12678|11316|10278|49518|44130|40146|3318|2976|2700|47828|96606|24389)\b|38\s*[;,/]\s*55|expectedNative|expectedCore|expectedEDI',
 'reference_paths':r'final_clean_run|incremental_run|accepted_states|evidence_exact|stage1_starting_state|accepted_stage1_source|tip_2h0_independent_run',
 'length_grid':r'TargetLength|MaxSegments|ismembertol|round\(|floor\(|ceil\(|\bmod\(|interp1|local_rows_at_lengths',
 'parameter_flow':r'IncrementMM|A0Override|HTipOverA0|rInner|rOuter|rCore|transitionLength|farCap|ExteriorCalibration',
 'reuse_provenance':r'Resume|Checkpoint|checkpoint|exteriorMeshControls|isReferenceProductionExterior|schemaVersion',
}
compiled={k:re.compile(v) for k,v in patterns.items()}
tracked=set(subprocess.check_output(['git','ls-files'],cwd=root,text=True).splitlines())
inventory=[];hits=[]
for p in sorted(root.rglob('*')):
 if not p.is_file() or p.suffix.lower() not in suffixes:continue
 rel=p.relative_to(root).as_posix()
 if '.git' in p.relative_to(root).parts or 'build' in p.relative_to(root).parts or '__pycache__' in p.parts:continue
 raw=p.read_bytes();text=raw.decode('utf-8-sig',errors='replace');lines=text.splitlines()
 role='data_snapshot' if ('/data/' in '/'+rel and p.suffix.lower() in {'.json','.yml','.yaml','.toml'}) else 'source_or_document'
 inventory.append({'path':rel,'extension':p.suffix.lower(),'tracked':rel in tracked,'role':role,'lines':len(lines),'sha256':hashlib.sha256(raw).hexdigest()})
 for n,line in enumerate(lines,1):
  categories=[k for k,pat in compiled.items() if pat.search(line)]
  if categories:hits.append({'path':rel,'line':n,'categories':';'.join(categories),'role':role,'text':line[:600]})
for name,rows in [('source_inventory.csv',inventory),('assumption_hits.csv',hits)]:
 with (out/name).open('w',newline='',encoding='utf-8') as f:
  w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)
summary={'base_commit':subprocess.check_output(['git','rev-parse','HEAD'],cwd=root,text=True).strip(),
 'files':len(inventory),'lines':sum(x['lines'] for x in inventory),'hits':len(hits),
 'by_extension':{s:sum(x['extension']==s for x in inventory) for s in sorted(suffixes)},
 'by_category':{k:sum(k in x['categories'].split(';') for x in hits) for k in patterns}}
(out/'scan_summary.json').write_text(json.dumps(summary,indent=2)+'\n')
print(json.dumps(summary,indent=2))
