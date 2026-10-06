"""Refresh only audited inline data in main.tex; no mechanics calculation."""
from pathlib import Path
import hashlib
import json
import re
import csv

HERE = Path(__file__).resolve().parent
e = json.loads((HERE / 'data/evidence.json').read_text(encoding='utf-8'))
t = e['stateRows']
c = e['codFits']
q = e['qualification']
s = e.get('stage1Summary')
if not s:
    raise SystemExit('Exact Stage-I summary missing: re-export using the accepted FrozenStateFile.')
assert len(t) == 23 and len(c) == 176 and len(q) == 23
assert e['audit']['qualifiedUnsolvedSegment'] == 24
assert all(x['segment'] <= 23 for x in t)

def number(x):
    return format(float(x), '.17g')

def sci(x, digits=5):
    a, b = format(float(x), f'.{digits-1}e').split('e')
    return rf'{a}\times10^{{{int(b)}}}'

def macro(name, value):
    return rf'\newcommand{{\{name}}}{{{value}}}'

last = t[-1]
v = e['vertices_m']
peak = t[16]
late = [x for x in c if x['segment'] == 23 and x['degree'] == 2]
near = next(x for x in c if x['segment'] == 21 and x['degree'] == 1
            and x['lower_r_over_DeltaA'] == .12 and x['upper_r_over_DeltaA'] == .3)
metrics = {
    'last_horizontal_ligament_mm': 1000*(e['C']['A']-last['tip_x_m']),
    'minimum_theta_deg': min(x['theta_deg'] for x in t),
    'minimum_theta_segment': min(t, key=lambda x:x['theta_deg'])['segment'],
    'COD_max_absolute_turn_difference_deg': max(abs(x['turn_error_deg']) for x in c),
    'COD_P23_quadratic_turn_min_deg': min(x['delta_theta_next_MTS_deg'] for x in late),
    'COD_P23_quadratic_turn_max_deg': max(x['delta_theta_next_MTS_deg'] for x in late),
    'COD_P21_rear_linear_ratio_relative_difference_percent': 100*near['ratio_error']/near['EDI_ratio'],
    'COD_P21_rear_linear_absolute_turn_difference_deg': abs(near['turn_error_deg']),
    'all_COD_signs_P21_positive': all(x['ratio_COD']>0 for x in c if x['segment']==21),
    'all_COD_signs_P22_negative': all(x['ratio_COD']<0 for x in c if x['segment']==22),
    'positive_KI_monotone': all(b['KI_unit']>a['KI_unit'] for a,b in zip(t,t[1:])),
    'max_tiny_mixed_relative_error': max(x[3] for x in e['audit']['syntheticErrors']),
    'source_evidence_sha256': hashlib.sha256((HERE/'data/evidence.json').read_bytes()).hexdigest(),
}
assert metrics['all_COD_signs_P21_positive'] and metrics['all_COD_signs_P22_negative']
assert metrics['positive_KI_monotone']

lines = [f'% Snapshot: {e["sourceCommit"]}; no P24 physical data.',
         f'% Evidence JSON SHA256: {metrics["source_evidence_sha256"]}']
values = {
    'PlateAMM': format(e['C']['A']*1000,'.8g'),
    'PlateBMM': format(e['C']['B']*1000,'.8g'),
    'HoleCxMM': format(s['hole_x_m']*1000,'.8g'),
    'HoleCyMM': format(s['hole_y_m']*1000,'.8g'),
    'HoleRMM': format(s['hole_R_m']*1000,'.8g'),
    'PhiStar': format(s['phi_star_deg'],'.10f'),
    'MouthX': format(s['x_star_m'],'.12g'), 'MouthY': format(s['y_star_m'],'.12g'),
    'NormalX': format(s['nmat_x'],'.12g'), 'NormalY': format(s['nmat_y'],'.12g'),
    'PeakStress': format(s['sigma_tt_peak_unit'],'.6f'),
    'LoadFactor': format(s['lambda_ini'],'.6f'),
    'LastXMM': format(last['tip_x_m']*1000,'.6f'),
    'LastYMM': format(last['tip_y_m']*1000,'.6f'),
    'NextXMM': number(v[24][0]*1000),'NextYMM': number(v[24][1]*1000),
    'PeakXMM': number(peak['tip_x_m']*1000),'PeakYMM': number(peak['tip_y_m']*1000),
    'LSXMM': number(e['interpolatedZeroPoint_m'][0]*1000),
    'LSYMM': number(e['interpolatedZeroPoint_m'][1]*1000),
    'LSLengthMM': format(e['interpolatedZeroLength_mm'],'.2f'),
    'LSPlotMM': number(e['interpolatedZeroLength_mm']),
    'LastLigamentMM': format(metrics['last_horizontal_ligament_mm'],'.4f'),
    'MinimumTheta': format(metrics['minimum_theta_deg'],'.6f'),
    'MaxCODError': format(metrics['COD_max_absolute_turn_difference_deg'],'.6f'),
    'CODQuadMin': format(metrics['COD_P23_quadratic_turn_min_deg'],'.6f'),
    'CODQuadMax': format(metrics['COD_P23_quadratic_turn_max_deg'],'.6f'),
    'CODNearZeroPercent': format(metrics['COD_P21_rear_linear_ratio_relative_difference_percent'],'.1f'),
    'CODNearZeroTurnDiff': format(metrics['COD_P21_rear_linear_absolute_turn_difference_deg'],'.6f'),
    'SyntheticRelativeError': sci(metrics['max_tiny_mixed_relative_error'],3),
    'MaxTrueResidual': sci(e['audit']['maxTrueResidual'],5),
}
lines += [macro(k,x) for k,x in values.items()]

def coords(points):
    return ' '.join(f'({number(x*1000)},{number(y*1000)})' for x,y in points)
lines += [macro('AcceptedPath',coords(v[:24])), macro('LateAcceptedPath',coords(v[15:24])),
          macro('UnsolvedPath',coords(v[23:25]))]

def table(name, header, rows):
    return '\\pgfplotstableread[col sep=space]{\n'+header+'\n'+'\n'.join(
        ' '.join(number(x) for x in row) for row in rows)+'\n}\\'+name

lines.append(table('historydata','a KI KII q theta turn',[
    [x['crack_length_mm'],x['KI_unit'],x['KII_unit'],x['KII_over_KI'],x['theta_deg'],x['delta_theta_next_deg']]
    for x in t]))
windows=[(.04,.20),(.04,.30),(.08,.30),(.12,.30)]
lookup={(x['segment'],x['degree'],x['lower_r_over_DeltaA'],x['upper_r_over_DeltaA']):x for x in c}
header=['a','qEDI','tEDI']+[f'{f}d{d}w{w}' for f in ['q','t'] for d in [1,2] for w in range(1,5)]
rows=[]
for x in t[1:]:
    row=[x['crack_length_mm'],x['KII_over_KI'],x['delta_theta_next_deg']]
    for field in ['ratio_COD','delta_theta_next_MTS_deg']:
        for d in [1,2]:
            for lo,hi in windows:
                row.append(lookup[(x['segment'],d,lo,hi)][field])
    rows.append(row)
lines.append(table('coddata',' '.join(header),rows))
qm={x['segment']:x for x in q}
lines.append(table('qualitydata','a angle ratio iter res',[
    [x['crack_length_mm'],qm[x['segment']]['min_angle_deg'],qm[x['segment']]['max_neighbor_ratio'],
     x['PCG_iterations'],x['true_rel_residual']] for x in t[1:]]))

colors=['pathblue','diagorange','fitgreen','fitpurple']
for deg, word in [(1,'One'),(2,'Two')]:
    turns=[];ratios=[];legends=[]
    for w,color in enumerate(colors,1):
        turns.append(rf'\addplot[{color},mark=*] table[x=a,y expr={{\thisrow{{td{deg}w{w}}}-\thisrow{{tEDI}}}}] {{\coddata}};')
        p=rf'\addplot[{color},mark=*] table[x=a,y expr={{100000*(\thisrow{{qd{deg}w{w}}}-\thisrow{{qEDI}})}}] {{\coddata}};'
        ratios.append(p)
        lo,hi=windows[w-1]
        legends += [p,rf'\addlegendentry{{$[{lo:.2f},{hi:.2f}]\Delta a$}}']
    zero=r'\addplot[black,dashed,thin,forget plot] coordinates {(60,0) (92,0)};'
    lines.append(macro('CODTurnPlots'+word,'\n'.join(turns+[zero])))
    lines.append(macro('CODRatioPlots'+word,'\n'.join(ratios+[zero])))
    if deg==2:lines.append(macro('CODRatioPlotsTwoLegend','\n'.join(legends+[zero])))

rows=[]
for k in [1,17,21,22,23]:
    x=t[k-1]
    rows.append(f'$P_{{{k}}}$ & {k*4} & {x["KI_unit"]:.6f} & ${sci(x["KII_unit"],5)}$ & '
                f'${sci(x["KII_over_KI"],5)}$ & {x["theta_deg"]:.6f} & {x["delta_theta_next_deg"]:+.6f} '+r'\\')
lines.append(macro('CharacteristicRows','\n'.join(rows)))

main=HERE/'main.tex'
text=main.read_text(encoding='utf-8')
start='% BEGIN AUDITED DATA: populated by sync_manuscript_data.py'
end='% END AUDITED DATA'
assert text.count(start)==text.count(end)==1
text=text.split(start)[0]+start+'\n'+'\n'.join(lines)+'\n'+end+text.split(end)[1]
main.write_text(text,encoding='utf-8')
(HERE/'data/derived_metrics.json').write_text(json.dumps(metrics,indent=2)+'\n',encoding='utf-8')
print('Updated standalone manuscript: 23 accepted states; P24 geometry only; no solve.')
