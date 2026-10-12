"""Before/after contact sheets and TeX proofs, using existing vector assets."""
import argparse
import json
import os
import re
import subprocess
from pathlib import Path
from PIL import Image,ImageDraw

PAPER=Path(__file__).resolve().parent


def prepare(root):
    root=Path(root);manifest=json.loads((root/'baseline_manifest.json').read_text())
    paths=[r['asset'] for r in manifest['assets']]
    for first in range(0,len(paths),8):
        sheet=Image.new('RGB',(1400,1500),'white');dr=ImageDraw.Draw(sheet)
        for j,asset in enumerate(paths[first:first+8]):
            x=(j%2)*700;y=(j//2)*375
            dr.text((x+5,y+4),asset.replace('figures/',''),fill='black')
            for k,state in enumerate(['before','after']):
                image=Image.open((root/state/asset).with_suffix('.png')).convert('RGB');image.thumbnail((340,325))
                sheet.paste(image,(x+k*350,y+45));dr.text((x+k*350+5,y+25),state,fill='black')
        sheet.save(root/f'comparison-{first//8+1}.png')
    before=subprocess.check_output(['git','show',manifest['baseCommit']+':paper/main.tex'],cwd=PAPER.parent,text=True)
    after=(PAPER/'main.tex').read_text()
    for state,text in [('before',before),('after',after)]:
        for asset in paths:text=text.replace('{'+asset+'}','{'+state+'/'+asset+'}')
        try:bib=os.path.relpath(PAPER/'references',root).replace('\\','/')
        except ValueError:bib=(PAPER/'references').as_posix()
        text=text.replace(r'\bibliography{references}',r'\bibliography{'+bib+'}')
        (root/f'manuscript_{state}.tex').write_text(text)
        preamble=text.split(r'\begin{document}')[0]
        figures=re.findall(r'\\begin\{figure\}.*?\\end\{figure\}',text,re.S)
        assert len(figures)==10
        preview=preamble+r'\begin{document}\setcounter{figure}{6}'+'\n'
        preview+='\n\\clearpage\n'.join(figures[6:9])+r'\end{document}'+'\n'
        (root/f'figures_7_9_{state}.tex').write_text(preview)
    print('Prepared four before/after panel contact sheets and full/native-size sensitivity TeX proofs')


if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('root');prepare(ap.parse_args().root)
