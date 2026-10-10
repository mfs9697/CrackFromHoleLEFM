"""Round manuscript displays only; accepted JSON/MAT/CSV evidence is untouched."""
from pathlib import Path
import json,re
paper=Path(__file__).resolve().parent
e=json.loads((paper/'data/evidence.json').read_text())
p=paper/'main.tex';text=p.read_text()
display={'PhiStar':'-1.5606','MouthX':'0.199989','MouthY':'-0.0208170',
         'NormalX':'0.999629','NormalY':'-0.0272345','PeakStress':'3.7141','LoadFactor':'80.774',
         'LastXMM':'291.840','LastYMM':'-25.579','LastLigamentMM':'8.16','MinimumTheta':'-3.646',
         'MaxCODError':'0.01024','CODQuadMin':'0.6204','CODQuadMax':'0.6225','CODNearZeroTurnDiff':'0.00245',
         'MaxTrueResidual':r'1.00\times10^{-10}'}
for name,value in display.items():
    pattern=r'(\\newcommand\{\\'+name+r'\}\{)[^}]+(\})'
    # Braced scientific notation needs an explicit complete command match.
    if name=='MaxTrueResidual':
        text=re.sub(r'\\newcommand\{\\MaxTrueResidual\}\{.*\}',lambda _:rf'\newcommand{{\MaxTrueResidual}}{{{value}}}',text,count=1)
    else:
        text,n=re.subn(pattern,lambda m:m[1]+value+m[2],text,count=1);assert n==1,name
def sci(x):
    a,b=f'{x:.3e}'.split('e');return a+rf'\times10^{{{int(b)}}}'
rows=[]
for k in (1,17,21,22,23):
    x=e['stateRows'][k-1]
    rows.append(f'$P_{{{k}}}$ & {4*k} & {x["KI_unit"]:.5g} & ${sci(x["KII_unit"])}$ & '
                f'${sci(x["KII_over_KI"])}$ & {x["theta_deg"]:.4f} & {x["delta_theta_next_deg"]:+.5f} '+r'\\')
start=text.index(r'\newcommand{\CharacteristicRows}')
end=text.index('% END AUDITED DATA',start)
text=text[:start]+r'\newcommand{\CharacteristicRows}{'+'\n'.join(rows)+'}\n'+text[end:]
head,body=text.split('% END AUDITED DATA',1)
for old,new in {'0.366480':'0.36648','1.244360':'1.2444','84.140272':'84.14027',
                '84.146370':'84.14637','84.149410':'84.14941',
                '84.146369912':'84.146370','84.146276379':'84.146276',
                '0.425410776939':'0.42541','-0.1141524':'-0.11415',
                '-0.1141583':'-0.11416','-0.1141509':'-0.11415'}.items():
    body=body.replace(old,new)
# Long data-derived angle/SIF literals lose spurious digits. Reproducibility
# constants defining the frozen mesh family retain their exact configuration.
keep={'0.00675308135','0.0135061627'}
def shorten(m):
    token=m[0]
    if token.lstrip('+-') in keep:return token
    result=format(float(token),'+.5g' if token.startswith('+') else '.5g')
    if 'e' in result:
        mantissa,exponent=result.split('e');return mantissa+rf'\times10^{{{int(exponent)}}}'
    return result
body=re.sub(r'(?<![\d.])[+-]?\d+\.\d{7,}(?!\d)',shorten,body)
body=re.sub(r'([+-]?\d+\.\d+)[eE]([+-]?\d+)',lambda m:m[1]+rf'\times10^{{{int(m[2])}}}',body)
body=body.replace('$-0.000000713$',r'$-7.13\times10^{-7}$')
body=body.replace('The final\ncolumn is a linear interpolation between solved mode-mixity values, not an\nindependently solved crack state.',
'''The final column is an interpolation between solved values, not a solved
crack state; its extra digits resolve numerical mesh-comparison shifts only.''')
p.write_text(head+'% END AUDITED DATA'+body,newline='\n')
print('Updated displayed precision only; full-precision evidence remains unchanged.')
