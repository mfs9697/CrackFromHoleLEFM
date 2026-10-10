"""Current portable manuscript verification; read-only and no numerical solve."""
from pathlib import Path
from pypdf import PdfReader
import csv,hashlib,json,re,runpy,math
paper=Path(__file__).resolve().parent;root=paper.parent
runpy.run_path(str(paper/'verify_isolated_tip_correction.py'))
text=(paper/'main.tex').read_text();abstract=text.split(r'\begin{abstract}')[1].split(r'\end{abstract}')[0]
assert 150<=len(abstract.split())<=200
assert not re.search(r'\d|\\cite|\\ref|CoreScale|ExteriorScale',abstract)
assert 'hidelinks' not in text and 'colorlinks=true' in text and 'referenceblue' in text
assert r'\input{figures/' not in text and len(re.findall(r'\\begin\{figure\}',text))==10
labels=re.findall(r'\\label\{([^}]+)\}',text);assert len(labels)==len(set(labels))
for target in re.findall(r'\\(?:eqref|ref)\{([^}]+)\}',text):assert target in labels,target
for symbol in (r'\Theta_k',r'\eta_c',r'\eta_e',r'\delta\theta',r'\dd\Omega'):assert symbol in text
assert 'CoreScale' not in text and 'ExteriorScale' not in text
assert r'\Omega_{\mathrm{rect}}=[0,2A]\times[-B/2,B/2]' in text
assert r'y_c-B/2' in text and r'2A-x_{23}' in text
assert r'E=210$ GPa' in text
assert not any(ord(x)<32 and x not in '\n\r\t' for x in text)
register=list(csv.DictReader((paper/'notation_register.csv').open()))
assert len(register)>=40 and all(x['first_source_line'] for x in register)
pdf=PdfReader(paper/'main.pdf');pages=len(pdf.pages);all_text=''.join(x.extract_text() for x in pdf.pages)
assert '??' not in all_text
links=[]
for page in pdf.pages:
    for item in page.get('/Annots',[]):
        annotation=item.get_object()
        if annotation.get('/Subtype')=='/Link':
            links.append(annotation)
            assert annotation.get('/Border',[0,0,0])[-1]==0
            dest=annotation.get('/Dest')
            if isinstance(dest,str):assert dest in pdf.named_destinations,dest
assert len(links)>40
schematic=PdfReader(paper/'figures/geometry_loading/specimen_geometry.pdf').pages[0]
assert not any(x.get_object().get('/Subtype')=='/Image' for x in schematic['/Resources'].get('/XObject',{}).values())
# The effective scale is determined from physical PDF pages, not FontSize alone.
for folder,target in [('geometry_loading',.78*162),('tip_resolution_sensitivity',.315*162)]:
    name='specimen_geometry' if folder=='geometry_loading' else 'vertical_deviation'
    page=PdfReader(paper/'figures'/folder/(name+'.pdf')).pages[0]
    width=float(page.mediabox.width)*25.4/72
    assert abs(width-target)<.15,(folder,width,target)
print(f'PASS: {len(abstract.split())}-word abstract; {len(register)} notation families; {len(labels)} unique labels; '
      f'{len(links)} working unboxed links; {pages} PDF pages; vector schematic and final physical export sizes.')
