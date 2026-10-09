"""Portable publication/table checks; no archive access or numerical experiments."""
from pathlib import Path
import csv,hashlib,json,math,re

paper=Path(__file__).resolve().parent
root=paper.parent
d=json.loads((paper/'data/isolated_tip_resolution.json').read_text())
manifest=json.loads((paper/'data/isolated_tip_resolution_manifest.json').read_text())
metrics=json.loads((paper/'figures/tip_resolution_sensitivity/isolated_tip_resolution_metrics.json').read_text())
reference=json.loads((paper/'data/evidence.json').read_text())['stateRows']
text=(paper/'main.tex').read_text();count=0
assert d['datasetId']=='isolated_core_tip_resolution' and d['exteriorScale']==1
assert d['authoritativeStudy']=='isolated_tip_resolution_study_20261009T200424169'
for item in manifest['portableFiles']+manifest['referenceInputs']:
    path=root/item['path']
    assert hashlib.sha256(path.read_bytes()).hexdigest()==item['sha256'],path
    count+=1
t=d['independent'];f=d['fixed']
assert len(t)==23 and len(f)==8 and [x['segment'] for x in t]==list(range(1,24))
for x in t+f:
    assert x['pass'] and x['qualificationPassed'] and x['syntheticPassed'] and all(x['physicalGates'].values())
    assert x['exterior_scale']==1 and x['pcg_flag']==0
    assert abs(x['KII_over_KI']-x['KII_unit']/x['KI_unit'])<1e-14
    count+=1
for x in t: assert x['core_scale']==2
for scale in (2,.5): assert [x['segment'] for x in f if x['core_scale']==scale]==[17,21,22,23]
dx=[1e6*(a['tip_x_m']-b['tip_x_m']) for a,b in zip(t,reference)]
dy=[1e6*(a['tip_y_m']-b['tip_y_m']) for a,b in zip(t,reference)]
dt=[1e3*(a['theta_deg']-b['theta_deg']) for a,b in zip(t,reference)]
dq=[a['KII_over_KI']-b['KII_over_KI'] for a,b in zip(t,reference)]
values={'maxAbsDy_um':max(map(abs,dy)),'maxTipSeparation_um':max(map(math.hypot,dx,dy)),
        'maxAbsDtheta_mdeg':max(map(abs,dt)),'maxAbsDq':max(map(abs,dq))}
for name,value in values.items(): assert abs(value-metrics[name])<1e-12;count+=1
def crossing(rows):
    for a,b in zip(rows,rows[1:]):
        if a['KII_over_KI']*b['KII_over_KI']<=0:
            return a['crack_length_mm']-(b['crack_length_mm']-a['crack_length_mm'])*a['KII_over_KI']/(b['KII_over_KI']-a['KII_over_KI'])
    raise AssertionError('No interpolated crossing')
assert abs(crossing(t)-metrics['isolatedLocalSymmetry_mm'])<1e-11
assert abs(crossing(reference)-metrics['referenceLocalSymmetry_mm'])<1e-11
count+=2
def sci(x,decimals=5,sign=False):
    mantissa,exponent=format(x,('+' if sign else '')+f'.{decimals}e').split('e')
    return mantissa+rf'\times10^{{{int(exponent)}}}'
table=text.split(r'\label{tab:tipmesh}')[1].split(r'\end{table}')[0]
for x in metrics['fixed']:
    required=[sci(x['q21'],sign=True),format(x['turn21_deg'],'.8f'),sci(x['q22']),
              format(x['turn22_deg'],'+.8f'),format(x['crossing_mm'],'.6f')]
    for value in required: assert value in table,value;count+=1
summary=text.split(r'\label{tab:independent-sensitivity}')[1].split(r'\end{table}')[0]
expected=(f"Isolated-core $2h_0$ & {values['maxAbsDy_um']:.5f} & {values['maxTipSeparation_um']:.5f} & "
          f"{values['maxAbsDtheta_mdeg']:.5f} & ${sci(values['maxAbsDq'],2)}$ & ${metrics['localSymmetryDifference_um']:.3f}$")
assert expected in summary;count+=1
assert r'Coarser exterior M1 & 0.06734 & 0.06742 & 0.07304 & $2.89\times10^{-7}$ & $-0.322$' in summary
assert '0.123~' not in text and 'tip-scaled' not in text and 'decouple tip and exterior' not in text
assert r'\texttt{ExteriorScale}=1' in text and 'connectivity need not be identical' in text
assert format(metrics['isolatedLocalSymmetry_mm'],'.9f') in text
assert format(metrics['fixed'][-1]['P17KIChange_percent'],'.5f') in text
assert sci(metrics['fixed'][-1]['P17qDifference'],2) in text
assert sci(metrics['fixed'][-1]['P17turnDifference_deg'],2) in text
assert 'not an independently solved' in text and 'qualified but unsolved' in text
plot=(root/'verification/crack_path/plot_tip2h0_vs_reference_publication.m').read_text()
assert 'tip2h0_states.csv' not in plot and 'tip_2h0_independent_run' not in plot
assert 'HistoricalInputRejected' in plot and 'load_isolated_tip_publication_data' in plot
inputs=re.findall(r'\\input\{(figures/[^}]+)\}',text)
assert inputs==['figures/mesh_levels/figure','figures/trajectory/figure','figures/intensities/figure',
                'figures/late/figure','figures/cod/figure','figures/tip_resolution_sensitivity/figure',
                'figures/m1_mesh_sensitivity/figure','figures/increment_sensitivity/figure','figures/quality/figure']
panel_count=0
for item in inputs:
    wrapper=(paper/(item+'.tex')).read_text()
    for asset in re.findall(r'\\includegraphics(?:\[[^]]*\])?\{([^}]+)\}',wrapper):
        assert (paper/asset).is_file(),asset;panel_count+=1
assert panel_count==30
print(f'PASS: {count} independent data/table checks; 31 accepted records; 9 figures / 30 panels; historical fallback absent.')
