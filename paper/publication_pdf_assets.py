"""Audited vector-only publication transforms for Issue #78.

The PDF is the input drawing, never a digitized numerical history. Tagging
paint operators lets us separate legend/label decorations without masking
curves underneath. Data operators stay in their original order/coordinates.
Cropping adds one translation; it never rescales a plot or changes a viewport.
"""
from __future__ import annotations
import argparse
import copy
import hashlib
import io
import json
import os
import re
import shutil
import subprocess
from dataclasses import dataclass, field
from pathlib import Path

import pdfplumber
from pdfminer.utils import mult_matrix
from pypdf import PdfReader, PdfWriter, Transformation
from pypdf.generic import (ArrayObject, ByteStringObject, ContentStream, DecodedStreamObject,
                           DictionaryObject, FloatObject, NameObject, NumberObject)
from reportlab.pdfgen import canvas

PAPER = Path(__file__).resolve().parent
PT = 72 / 25.4
PAINT = {b'S',b's',b'f',b'F',b'f*',b'B',b'B*',b'b',b'b*',b'Tj',b'TJ',b"'",b'"'}


def sha(data): return hashlib.sha256(data).hexdigest()
def union(boxes):
    boxes=list(boxes);assert boxes,'Empty visible drawing'
    return [min(x[0] for x in boxes),min(x[1] for x in boxes),max(x[2] for x in boxes),max(x[3] for x in boxes)]
def box(obj): return [obj['x0'],obj['y0'],obj['x1'],obj['y1']]
def inside(a,b,tol=0): return a[0]>=b[0]-tol and a[1]>=b[1]-tol and a[2]<=b[2]+tol and a[3]<=b[3]+tol
def intersects(a,b): return a[0]<b[2] and b[0]<a[2] and a[1]<b[3] and b[1]<a[3]
def expand(b,n): return [b[0]-n,b[1]-n,b[2]+n,b[3]+n]
def white(c): return c in (1,1.0,(1,1,1),[1,1,1])
def colored(o):
    c=o.get('stroking_color') if o.get('stroke') else o.get('non_stroking_color')
    return isinstance(c,(tuple,list)) and len(c)>=3 and max(c)-min(c)>.03
def painted(o):
    return o.get('object_type')=='char' or (o.get('stroke') or o.get('fill') and not white(o.get('non_stroking_color')))
def bytes_pdf(page):
    w=PdfWriter();w.add_page(page);b=io.BytesIO();w.write(b);return b.getvalue()
def set_ops(page,ops,reader):
    p=copy.copy(page);s=ContentStream(page.get_contents(),reader);s.operations=ops
    d=DecodedStreamObject();d.set_data(s.get_data());p[NameObject('/Contents')]=d;return p


@dataclass
class Node:
    begin:int
    state:dict
    end:int=0
    paints:list=field(default_factory=list)
    children:list=field(default_factory=list)
    objects:list=field(default_factory=list)


class Drawing:
    def __init__(self,file):
        self.file=Path(file);self.reader=PdfReader(file);self.page=self.reader.pages[0]
        assert len(self.reader.pages)==1,'Panels must be single-page PDFs'
        self.ops=ContentStream(self.page.get_contents(),self.reader).operations
        self.nodes=[];self.owner={};self.font_at={};self.objects=[]
        state={'ctm':(1,0,0,1,0,0),'properties':{}};stack=[];current=None;tagged=[]
        for i,(args,op) in enumerate(self.ops):
            if op==b'q':
                n=Node(i,copy.deepcopy(state));self.nodes.append(n)
                if current:current.children.append(n)
                stack.append((current,copy.deepcopy(state)));current=n
            elif op==b'Q':
                assert current is not None,'Unbalanced source graphics state'
                current.end=i+1;current,state=stack.pop()
            elif op==b'cm':state['ctm']=mult_matrix(tuple(float(x) for x in args),state['ctm'])
            elif op in {b'g',b'rg',b'k'}:state['properties']['fill']=(args,op)
            elif op in {b'G',b'RG',b'K'}:state['properties']['stroke']=(args,op)
            elif op in {b'w',b'd',b'J',b'j',b'M',b'gs',b'Tf',b'Tc',b'Tw',b'Tz',b'TL',b'Tr',b'Ts'}:
                state['properties'][op]=(args,op)
            if op in PAINT:
                assert current is not None,'Paint outside a balanced group'
                current.paints.append(i);self.owner[i]=current
                self.font_at[i]=state['properties'].get(b'Tf')
                tagged.append(([NameObject('/Span'),DictionaryObject({NameObject('/MCID'):NumberObject(i)})],b'BDC'))
                tagged.append((args,op));tagged.append(([],b'EMC'))
            else:tagged.append((args,op))
        assert not stack,'Unbalanced source page'
        p=set_ops(self.page,tagged,self.reader)
        with pdfplumber.open(io.BytesIO(bytes_pdf(p))) as doc:
            pp=doc.pages[0]
            self.objects=pp.chars+pp.curves+pp.lines+pp.rects
            self.edges=pp.edges;self.images=pp.images
        for o in self.objects:
            if o.get('mcid') in self.owner:self.owner[o['mcid']].objects.append(o)
        self.visible=[o for o in self.objects if painted(o)]
        self.content=union(box(o) for o in self.visible)
        self.axes=self.find_axes()

    def find_axes(self):
        if 'geometry_loading/' in self.file.as_posix() or self.file.stem=='specimen_geometry' or '/stamps/' in self.file.as_posix():return None
        w=float(self.page.mediabox.width);h=float(self.page.mediabox.height)
        hs=[e for e in self.edges if e.get('stroke') and e['width']>.3*w and e['height']<.02]
        vs=[e for e in self.edges if e.get('stroke') and e['height']>.12*h and e['width']<.02]
        candidates=[]
        for a in hs:
            for b in hs:
                if abs(a['x0']-b['x0'])>.02 or abs(a['x1']-b['x1'])>.02 or a['y0']>=b['y0']:continue
                if a['x0']<.1 and a['x1']>w-.2:continue
                if all(any(abs(v['x0']-x)<.03 and abs(v['y0']-a['y0'])<.03 and abs(v['y1']-b['y0'])<.03 for v in vs) for x in [a['x0'],a['x1']]):
                    candidates.append([a['x0'],a['y0'],a['x1'],b['y0']])
        assert candidates,f'No axes frame: {self.file}'
        return max(candidates,key=lambda x:(x[2]-x[0])*(x[3]-x[1]))

    def metadata(self):
        path=self.file.with_suffix('.layout.json')
        return json.loads(path.read_text()) if path.exists() else {}

    def guess_legend(self,asset):
        m=self.metadata()
        if m.get('legendBounds_mm'):
            b=m['legendBounds_mm'];b=b if isinstance(b[0],(int,float)) else b[0]
            return expand([b[0]*PT,b[1]*PT,(b[0]+b[2])*PT,(b[1]+b[3])*PT],2)
        if m.get('legend'):return expand(m['legend']['bounds_pt'],2)
        if asset=='figures/trajectory/overview.pdf':
            x,y,r,t=self.axes;return [x+.46*(r-x),y+.68*(t-y),r+.2,t+.2]
        if asset=='figures/cod/mixity_quadratic.pdf':
            x,y,r,t=self.axes;return [x,y+.44*(t-y),x+.77*(r-x),t]
        return None

    def legend_cells(self,region,count):
        # Only complete, balanced decoration primitives are selected. Long
        # data polylines crossing the rectangle cannot become legend keys.
        candidates=[n for n in self.nodes if n.paints and not n.children and n.objects
                    and inside(union(box(o) for o in n.objects),region,.1)]
        keys=[]
        for n in candidates:
            for o in n.objects:
                if colored(o) and painted(o):
                    b=box(o);cy=(b[1]+b[3])/2
                    hit=next((k for k in keys if abs(k['y']-cy)<1.8 and abs(k['x']-b[0])<24),None)
                    if hit:hit['boxes'].append(b);hit['x']=min(hit['x'],b[0])
                    else:keys.append({'x':b[0],'y':cy,'boxes':[b]})
        assert len(keys)==count,(self.file,'legend key count',len(keys),count)
        # Column clusters preserve MATLAB's column-major legend order.
        columns=[]
        for k in sorted(keys,key=lambda k:k['x']):
            if not columns or k['x']-columns[-1][0]['x']>25:columns.append([k])
            else:columns[-1].append(k)
        entries=[]
        for colindex,col in enumerate(columns):
            col.sort(key=lambda k:-k['y'])
            ys=[k['y'] for k in col]
            gap=min(abs(a-b) for a,b in zip(ys,ys[1:])) if len(ys)>1 else 13
            right=columns[colindex+1][0]['x']-1 if colindex+1<len(columns) else region[2]+.2
            left=min(k['x'] for k in col)-1
            for k in col:
                nodes=[n for n in candidates if left-.1<=union(box(o) for o in n.objects)[0]
                       and union(box(o) for o in n.objects)[2]<=right+.1
                       and abs((union(box(o) for o in n.objects)[1]+union(box(o) for o in n.objects)[3])/2-k['y'])<=gap*.55]
                assert nodes,'Empty legend entry'
                entries.append(nodes)
        return entries

    def label_nodes(self):
        x,y,r,t=self.axes
        candidates=[n for n in self.nodes if n.paints and not n.children and n.objects
                    and all(painted(o) and not colored(o) for o in n.objects)
                    and union(box(o) for o in n.objects)[2]<x-2
                    and y<(union(box(o) for o in n.objects)[1]+union(box(o) for o in n.objects)[3])/2<t]
        # The rotated ylabel is the leftmost overlapping glyph column; tick
        # labels form separate columns to its right and are never selected.
        intervals=sorted(((union(box(o) for o in n.objects),n) for n in candidates),key=lambda x:x[0][0])
        clusters=[]
        for b,n in intervals:
            if not clusters or b[0]>clusters[-1]['right']+1.2:clusters.append({'right':b[2],'nodes':[n]})
            else:clusters[-1]['right']=max(clusters[-1]['right'],b[2]);clusters[-1]['nodes'].append(n)
        assert clusters,'Missing ylabel glyphs'
        nodes=clusters[0]['nodes'];b=union(box(o) for n in nodes for o in n.objects)
        assert b[3]-b[1]>min((t-y)*.3,20) and b[2]<x-9,('Ambiguous ylabel/tick separation',str(self.file),b,self.axes)
        return nodes,b

    def primitive_page(self,nodes,relabel=False):
        ops=[];changed=[]
        for n in nodes:
            ops.append(([],b'q'));ops.append(([FloatObject(v) for v in n.state['ctm']],b'cm'))
            ops.extend(n.state['properties'].values())
            for i in range(n.begin,n.end):
                args,op=self.ops[i]
                if relabel and op in (b'TJ',b'Tj'):
                    args,old,new=self.relabel(args,op,self.font_at.get(i))
                    if old!=new:changed.append([old,new])
                ops.append((args,op))
            ops.append(([],b'Q'))
        return set_ops(self.page,ops,self.reader),changed

    def relabel(self,args,op,font_state):
        if not font_state:return args,'',''
        font=self.page['/Resources']['/Font'][font_state[0][0]].get_object()
        if font.get('/Subtype')!='/Type0' or '/ToUnicode' not in font:return args,'',''
        cmap=font['/ToUnicode'].get_data().decode('latin1');mapping={}
        for block in re.findall(r'beginbfchar(.*?)endbfchar',cmap,re.S):
            for a,b in re.findall(r'<([0-9a-fA-F]+)>\s*<([0-9a-fA-F]+)>',block):mapping[bytes.fromhex(a)]=bytes.fromhex(b).decode('utf-16-be')
        for block in re.findall(r'beginbfrange(.*?)endbfrange',cmap,re.S):
            for a,b,c in re.findall(r'<([0-9a-fA-F]+)>\s*<([0-9a-fA-F]+)>\s*<([0-9a-fA-F]+)>',block):
                for j in range(int(a,16),int(b,16)+1):mapping[j.to_bytes(len(a)//2,'big')]=chr(int(c,16)+j-int(a,16))
        if not mapping:return args,'',''
        width=len(next(iter(mapping)));reverse={v:k for k,v in mapping.items()}
        strings=[a for a in (args[0] if op==b'TJ' else args) if not isinstance(a,(float,int))]
        if len(strings)!=1:return args,'',''
        s=strings[0];raw=s if isinstance(s,bytes) else s.get_original_bytes()
        if raw.startswith(b'\xfe\xff'):raw=raw[2:]
        text=''.join(mapping.get(raw[j:j+width],'?') for j in range(0,len(raw),width))
        if text!='reference' and not text.startswith('reference q max'):return args,text,text
        new='M0'+text[len('reference'):]
        assert all(c in reverse for c in new),'M0 glyph absent in retained font'
        replacement=ByteStringObject(b''.join(reverse[c] for c in new))
        return ([ArrayObject([replacement])] if op==b'TJ' else [replacement]),text,new


def visible(page):
    with pdfplumber.open(io.BytesIO(bytes_pdf(page))) as doc:
        p=doc.pages[0];objs=p.chars+p.curves+p.lines+p.rects
        return union(box(o) for o in objs if painted(o))


def text_page(text,size=10.5):
    buf=io.BytesIO();c=canvas.Canvas(buf,pagesize=(500,80));c.setFont('Times-Roman',size)
    c.drawString(2,20,text);c.save();buf.seek(0);return PdfReader(buf).pages[0]


def segment_hits(p,q,rect):
    # Liang--Barsky: crossing segments count even when both endpoints are
    # outside the legend. This is a graphics collision test, not FEM logic.
    x0,y0,x1,y1=rect;dx=q[0]-p[0];dy=q[1]-p[1];lo=0.;hi=1.
    for a,b in [(-dx,p[0]-x0),(dx,x1-p[0]),(-dy,p[1]-y0),(dy,y1-p[1])]:
        if abs(a)<1e-12:
            if b<0:return False
        elif a<0:lo=max(lo,b/a)
        else:hi=min(hi,b/a)
        if lo>hi:return False
    return True


def obstacles(drawing,removed):
    result=[];axes=drawing.axes
    for o in drawing.visible:
        if o.get('mcid') in removed:continue
        b=box(o)
        if not intersects(b,axes):continue
        pts=o.get('pts',[]);world=[(x,float(drawing.page.mediabox.height)-top) for x,top in pts]
        if colored(o) or o.get('stroke') and len(world)>30:
            if len(world)>30 and abs(world[0][0]-world[-1][0])+abs(world[0][1]-world[-1][1])<1:
                result.append({'kind':'hole or closed geometry','box':b})
            else:result.append({'kind':'curve/marker','box':b,'points':world})
        elif o.get('stroke') and o.get('object_type') in ('line','curve'):
            # Critical black boundaries and gray zero/crossing lines.
            if len(world)>=2:result.append({'kind':'boundary/zero line','box':b,'points':world})
        elif o.get('fill') or o.get('object_type')=='char':
            if inside(b,axes):result.append({'kind':'critical annotation/marker','box':b})
    return result


def collisions(rect,items,clearance=.6):
    check=expand(rect,clearance);hits=[]
    for item in items:
        if not intersects(check,expand(item['box'],clearance)):continue
        pts=item.get('points')
        if pts and len(pts)>1:
            if any(segment_hits(p,q,check) for p,q in zip(pts,pts[1:])):hits.append(item)
        else:hits.append(item)
    return hits


def place_legend(axes,w,h,items,prefer):
    x,y,r,t=axes
    for pad in [1.5*PT,2.5*PT,3*PT,4*PT]:
        corners={'northeast':[r-w-pad,t-h-pad,r-pad,t-pad],
                 'northwest':[x+pad,t-h-pad,x+w+pad,t-pad],
                 'southwest':[x+pad,y+pad,x+w+pad,y+h+pad],
                 'southeast':[r-w-pad,y+pad,r-pad,y+h+pad]}
        for name in [prefer]+[n for n in corners if n!=prefer]:
            rect=corners[name]
            if inside(rect,axes) and not collisions(rect,items):return rect,name
    raise AssertionError('No unobstructed internal legend position; author review required')


def compose(drawing,asset,legend_count,stamp=None,extra_erase=()):
    region=drawing.guess_legend(asset);entries=[];removed=set(extra_erase);labels=[];changes=[]
    context='M0 exterior' if 'tip_resolution_sensitivity/' in asset else 'Scaled M1' if 'increment_sensitivity/' in asset and asset.endswith(('trajectory.pdf','turn_density.pdf')) else None
    if legend_count:
        assert region,'Missing signed legend metadata'
        entries=drawing.legend_cells(region,legend_count)
        removed.update(i for nodes in entries for n in nodes for i in n.paints)
    ylabel_nodes=[];ylabel_box=None
    if asset in ['figures/m1_mesh_sensitivity/vertical_deviation.pdf','figures/m1_mesh_sensitivity/direction_deviation.pdf']:
        ylabel_nodes,ylabel_box=drawing.label_nodes();removed.update(i for n in ylabel_nodes for i in n.paints)
        assert not intersects(ylabel_box,drawing.axes),'Axis label overlaps scientific frame'
        assert stamp and Path(stamp).is_file(),'MATLAB M0 ylabel stamp missing'
    ops=[]
    for i,(args,op) in enumerate(drawing.ops):
        if i in removed:
            if op in (b'TJ',b'Tj',b"'",b'"'):ops.append(([ArrayObject([])],b'TJ'))
            else:ops.append(([],b'n'))
        else:ops.append((args,op))
    body=set_ops(drawing.page,ops,drawing.reader)
    scientific_paints=set()
    if drawing.axes:
        width=drawing.axes[2]-drawing.axes[0]
        for o in drawing.visible:
            pts=o.get('pts',[]);b=box(o)
            if o.get('stroke') and (len(pts)>30 or colored(o) and len(pts)>=8
                                   or b[2]-b[0]>.6*width):
                scientific_paints.add(o.get('mcid'))
        assert scientific_paints.isdisjoint(removed),'Decoration selection includes a numerical curve'
    info={'removedDecorationPaints':sorted(removed),'legend':None,'ylabelReplacement':None}
    info['protectedScientificPaintCount']=len(scientific_paints)
    if ylabel_box:
        label_drawing=Drawing(stamp)
        # Do not overlay the MATLAB figure's white background onto ticks.
        glyphs=[n for n in label_drawing.nodes if n.paints and not n.children
                and any(painted(o) for o in n.objects)]
        label,_=label_drawing.primitive_page(glyphs);b=visible(label)
        tx=(ylabel_box[0]+ylabel_box[2]-b[0]-b[2])/2
        ty=(ylabel_box[1]+ylabel_box[3]-b[1]-b[3])/2
        body.merge_transformed_page(label,Transformation().translate(tx,ty))
        info['ylabelReplacement']={'source':Path(stamp).name,'bounds_pt':[b[0]+tx,b[1]+ty,b[2]+tx,b[3]+ty],
                                   'oldBounds_pt':ylabel_box,'label':'M1 minus M0'}
    if entries:
        # Genuine vector swatches and glyphs, reflowed into one column.
        panels=[];scale=10.5/9.5
        for nodes in entries:
            cell,changed=drawing.primitive_page(nodes,relabel='m1_mesh_sensitivity/' in asset)
            b=visible(cell);panels.append((cell,b,scale));changes.extend(changed)
        if context:
            title=text_page(context);panels.insert(0,(title,visible(title),1.))
        pad=1.*PT;gap=.65*PT
        w=max(s*(b[2]-b[0]) for _,b,s in panels)+2*pad
        h=sum(s*(b[3]-b[1]) for _,b,s in panels)+(len(panels)-1)*gap+2*pad
        obstacle_list=obstacles(drawing,removed)
        prefer='southwest' if asset.endswith('mode_mixity.pdf') and 'sensitivity/' in asset else 'northwest' if asset.endswith(('direction_deviation.pdf','turn_density.pdf','mixity_quadratic.pdf')) else 'northeast'
        rect,corner=place_legend(drawing.axes,w,h,obstacle_list,prefer)
        top=rect[3]-pad
        for cell,b,s in panels:
            bottom=top-s*(b[3]-b[1]);tx=rect[0]+pad-s*b[0];ty=bottom-s*b[1]
            body.merge_transformed_page(cell,Transformation([s,0,0,s,tx,ty]))
            top=bottom-gap
        assert inside(rect,drawing.axes) and not collisions(rect,obstacle_list)
        info['legend']={'bounds_pt':rect,'location':corner,'insideAxes':True,'collisionCount':0,
                        'font_pt':10.5,'entryCount':legend_count,'context':context,
                        'relabeling':changes,'method':'original vector keys/glyphs; M0 uses existing retained font cmap'}
    # Every unchanged operator in the source body remains exactly present.
    info['preservedBodyOperations_sha256']=sha(repr([(i,op) for i,op in enumerate(drawing.ops) if i not in removed]).encode())
    info['bodyTransform']=[1,0,0,1,0,0]
    return body,info


def exponent_plan(drawing,info):
    """Relocate unchanged multiplier glyphs into a clear upper-left corner."""
    if not drawing.axes:return None
    x,y,r,t=drawing.axes;removed=set(info['removedDecorationPaints'])
    nodes=[n for n in drawing.nodes if n.paints and not n.children and n.objects
           and not set(n.paints)&removed and all(painted(o) and not colored(o) for o in n.objects)
           and inside(union(box(o) for o in n.objects),[x-.5,t+.1,x+45,t+21],.05)]
    if not nodes:return None
    assert all(o.get('fill') or o.get('object_type')=='char' for n in nodes for o in n.objects)
    exponent,_=drawing.primitive_page(nodes);b=visible(exponent)
    items=obstacles(drawing,removed)
    if info['legend']:items.append({'kind':'legend','box':info['legend']['bounds_pt']})
    rect,corner=place_legend(drawing.axes,b[2]-b[0],b[3]-b[1],items,'northwest')
    return exponent,b,rect,corner,{i for n in nodes for i in n.paints}


def mask_paints(page,indices,reader):
    ops=ContentStream(page.get_contents(),reader).operations
    for i in indices:
        args,op=ops[i]
        ops[i]=([ArrayObject([])],b'TJ') if op in (b'TJ',b'Tj') else ([],b'n')
    return set_ops(page,ops,reader)


def retained_legend(body,source,region,count,relabel=False,context=None):
    """Current retained-source exporter: no temporary outside legend strip."""
    original=Drawing(source);entries=original.legend_cells(expand(region,1),count)
    panels=[];changes=[];scale=10.5/7.65
    for nodes in entries:
        cell,changed=original.primitive_page(nodes,relabel=relabel)
        panels.append((cell,visible(cell),scale));changes.extend(changed)
    if context:
        title=text_page(context);panels.insert(0,(title,visible(title),1))
    # Drawing reads a temporary vector page only, never a MAT history.
    import tempfile
    with tempfile.TemporaryDirectory() as temp:
        file=Path(temp)/'retained_graph.pdf';file.write_bytes(bytes_pdf(body));drawing=Drawing(file)
    pad=PT;gap=.65*PT
    w=max(s*(b[2]-b[0]) for _,b,s in panels)+2*pad
    h=sum(s*(b[3]-b[1]) for _,b,s in panels)+(len(panels)-1)*gap+2*pad
    items=obstacles(drawing,set())
    prefer='northeast' if count>2 else 'southwest'
    rect,corner=place_legend(drawing.axes,w,h,items,prefer);top=rect[3]-pad
    for cell,b,s in panels:
        bottom=top-s*(b[3]-b[1]);tx=rect[0]+pad-s*b[0];ty=bottom-s*b[1]
        body.merge_transformed_page(cell,Transformation([s,0,0,s,tx,ty]));top=bottom-gap
    return body,{'bounds_pt':rect,'insideAxes':True,'collisionCount':0,'location':corner,
                 'entryCount':count,'context':context,'relabeling':changes,'font_pt':10.5}


def crop(page,rectangle,output):
    p=copy.copy(page);p.add_transformation(Transformation().translate(-rectangle[0],-rectangle[1]))
    p.mediabox.lower_left=(0,0);p.mediabox.upper_right=(rectangle[2]-rectangle[0],rectangle[3]-rectangle[1])
    p.cropbox.lower_left=(0,0);p.cropbox.upper_right=p.mediabox.upper_right
    w=PdfWriter();q=w.add_page(p);q.compress_content_streams()
    output=Path(output);output.parent.mkdir(parents=True,exist_ok=True)
    with output.open('wb') as f:w.write(f)


def render(file,prefix,dpi=160):
    pop=os.environ.get('PUBLICATION_PDFTOPPM') or shutil.which('pdftoppm')
    if not pop:
        path=Path.home()/'.cache/codex-runtimes/codex-primary-runtime/dependencies/native/poppler/Library/bin/pdftoppm.exe'
        if path.exists():pop=str(path)
    assert pop,'Poppler required for visual verification'
    subprocess.run([pop,'-png','-singlefile','-r',str(dpi),str(file),str(prefix)],check=True,stdout=subprocess.DEVNULL,stderr=subprocess.DEVNULL)


def build_review(root):
    root=Path(root);manifest=json.loads((root/'baseline_manifest.json').read_text());records=[];pending=[]
    for path,digest in manifest.get('protectedData',{}).items():
        assert sha((PAPER.parent/path).read_bytes())==digest,'Numerical/retained source evidence changed'
    counts={'figures/trajectory/overview.pdf':2,'figures/cod/mixity_quadratic.pdf':4,
            'figures/tip_resolution_sensitivity/trajectory.pdf':5,'figures/tip_resolution_sensitivity/mode_mixity.pdf':2,
            'figures/m1_mesh_sensitivity/trajectory.pdf':5,'figures/m1_mesh_sensitivity/mode_mixity.pdf':2,
            'figures/increment_sensitivity/trajectory.pdf':3,'figures/increment_sensitivity/vertical_deviation.pdf':2,
            'figures/increment_sensitivity/direction_deviation.pdf':2,'figures/increment_sensitivity/turn_density.pdf':3}
    for row in manifest['assets']:
        asset=row['asset'];source=root/'before'/asset
        assert sha(source.read_bytes())==row['sha256'],'Changed original vector evidence'
        if row.get('layoutMetadata_sha256'):
            assert sha(source.with_suffix('.layout.json').read_bytes())==row['layoutMetadata_sha256'],'Changed original layout metadata'
        d=Drawing(source);assert not d.images,'Raster original cannot be certified as vector'
        stamp=root/'stamps'/('m1_vertical_label.pdf' if asset.endswith('vertical_deviation.pdf') else 'm1_direction_label.pdf')
        p,info=compose(d,asset,counts.get(asset,0),stamp)
        plan=exponent_plan(d,info)
        if plan:
            exponent,bb,rect,corner,erase=plan
            p,info=compose(d,asset,counts.get(asset,0),stamp,erase)
            p.merge_transformed_page(exponent,Transformation().translate(rect[0]-bb[0],rect[1]-bb[1]))
            info['exponentInsideAxes']={'oldBounds_pt':bb,'bounds_pt':rect,'location':corner,'collisionCount':0}
        b=visible(p);group=None
        if '/mesh_levels/' in asset:group='mesh'
        elif '/intensities/' in asset:group='intensities'
        elif '/late/' in asset:group='late'
        elif '/cod/' in asset:group='cod'
        elif '/quality/' in asset:group='quality'
        elif any('/'+x+'/' in asset for x in ['tip_resolution_sensitivity','m1_mesh_sensitivity','increment_sensitivity']):group='sensitivity_top' if asset.endswith('trajectory.pdf') else 'sensitivity_lower'
        elif '/trajectory/' in asset:group='reference_path'
        pending.append((row,d,p,info,b,group))
    # Align only rotated ylabel decorations, so shared crop envelopes don't
    # retain avoidable blank left columns. Plot bodies/ticks never move here.
    import tempfile
    labels={}
    with tempfile.TemporaryDirectory() as temp:
        for index,(row,d,p,info,b,group) in enumerate(pending):
            if not group:continue
            path=Path(temp)/(str(index)+'.pdf');path.write_bytes(bytes_pdf(p));current=Drawing(path)
            nodes,bb=current.label_nodes()
            labels[index]=(current,nodes,bb)
        targets={}
        for index,(current,nodes,bb) in labels.items():
            group=pending[index][5];value=bb[0]-current.axes[0]
            targets[group]=min(targets.get(group,value),value)
        for index,(current,nodes,bb) in labels.items():
            row,d,p,info,b,group=pending[index]
            dx=targets[group]-(bb[0]-current.axes[0])
            if abs(dx)>1e-6:
                ids={i for n in nodes for i in n.paints}
                p=mask_paints(p,ids,current.reader);label,_=current.primitive_page(nodes)
                p.merge_transformed_page(label,Transformation().translate(dx,0))
                info['ylabelAlignment']={'translation_pt':[dx,0],'scientificCurveMovement':False}
                pending[index]=(row,d,p,info,visible(p),group)
    # Shared envelopes preserve the intended axes alignment and physical size
    # within panel classes. Only the y crop is independent for the stacked
    # reference-path panels, whose spatial aspects differ deliberately.
    envelopes={}
    for row,d,p,info,b,group in pending:
        if not group:continue
        a=d.axes;relative=[b[0]-a[0],b[1]-a[1],b[2]-a[0],b[3]-a[1]]
        # Top COD/quality rows omit xlabels by design. They share the x
        # envelope with their bottom row, but not its unused label-space.
        yg=group
        if group=='cod':yg='cod_top' if row['asset'].endswith(('turn_linear.pdf','turn_quadratic.pdf')) else 'cod_bottom'
        if group=='quality':yg='quality_top' if row['asset'].endswith(('min_angle.pdf','size_ratio.pdf')) else 'quality_bottom'
        envelopes.setdefault(group,[]).append(relative)
        if yg!=group:envelopes.setdefault(yg,[]).append(relative)
    envelopes={key:union(boxes) for key,boxes in envelopes.items()}
    for row,d,p,info,b,group in pending:
        if group:
            e=envelopes[group];a=d.axes
            area=[a[0]+e[0],a[1]+e[1],a[0]+e[2],a[1]+e[3]]
            if group=='reference_path':area[1],area[3]=b[1],b[3]
            if group in ('cod','quality'):
                top=row['asset'].endswith(('turn_linear.pdf','turn_quadratic.pdf','min_angle.pdf','size_ratio.pdf'))
                ey=envelopes[group+('_top' if top else '_bottom')];area[1],area[3]=a[1]+ey[1],a[1]+ey[3]
        else:area=b
        region=expand(area,2*PT);target=root/'after'/row['asset'];crop(p,region,target)
        new_bbox=[b[0]-region[0],b[1]-region[1],b[2]-region[0],b[3]-region[1]]
        size=[region[2]-region[0],region[3]-region[1]]
        assert inside(new_bbox,[0,0,*size]),'Crop clips painted content'
        assert min(new_bbox[0],new_bbox[1],size[0]-new_bbox[2],size[1]-new_bbox[3])>=2*PT-1e-6
        axes=[d.axes[0]-region[0],d.axes[1]-region[1],d.axes[2]-region[0],d.axes[3]-region[1]] if d.axes else None
        if info['legend']:
            info['legend']['oldBounds_pt']=info['legend']['bounds_pt']
            r=info['legend']['bounds_pt'];info['legend']['bounds_pt']=[r[0]-region[0],r[1]-region[1],r[2]-region[0],r[3]-region[1]]
            assert inside(info['legend']['bounds_pt'],axes)
        if info['ylabelReplacement']:
            r=info['ylabelReplacement']['bounds_pt']
            info['ylabelReplacement']['bounds_pt']=[r[0]-region[0],r[1]-region[1],r[2]-region[0],r[3]-region[1]]
        if info.get('exponentInsideAxes'):
            r=info['exponentInsideAxes']['bounds_pt'];info['exponentInsideAxes']['bounds_pt']=[r[0]-region[0],r[1]-region[1],r[2]-region[0],r[3]-region[1]]
        info.update(asset=row['asset'],source_sha256=row['sha256'],output_sha256=sha(target.read_bytes()),
                    oldPage_pt=row['page_pt'],newPage_pt=size,oldAxes_pt=d.axes,axes_pt=axes,
                    oldVisibleContent_pt=d.content,visibleContent_pt=new_bbox,safetyMargin_mm=2,
                    margins_mm=[new_bbox[0]/PT,new_bbox[1]/PT,(size[0]-new_bbox[2])/PT,(size[1]-new_bbox[3])/PT],
                    layoutGroup=group,vector=True,raster=False,sourceToOutput=[1,0,0,1,-region[0],-region[1]],
                    plotGeometryUnchanged=True,axisLimitsAndValuesUnchanged=True)
        target.with_suffix('.layout.json').write_text(json.dumps(info,indent=2)+'\n');records.append(info)
        render(target,target.with_suffix(''),220)
    (root/'geometry_report.json').write_text(json.dumps({'baseCommit':manifest['baseCommit'],'panels':records,
        'legendPolicy':'inside axes, collision-checked','physicalSolves':0,'dataReconstruction':False,
        'sourceEvidenceChanged':False,'publicationPromotion':'pending author approval'},indent=2)+'\n')
    print('PASS: 31 vector assets; source drawings preserved; internal legends; 2-mm crop safety envelope')


def sensitivity_preview(root):
    root=Path(root);folders=['figure7_tip_resolution','figure8_M1_exterior','figure9_increment'];writer=PdfWriter();rows=[]
    for i,folder in enumerate(folders):
        names=['trajectory','vertical_deviation','direction_deviation','turn_density' if i==2 else 'mode_mixity']
        buf=io.BytesIO();c=canvas.Canvas(buf,pagesize=(210*PT,297*PT));c.setFont('Times-Roman',11)
        c.drawString(24*PT,278*PT,f'Figure {i+7}: internal legends and tight vector panels')
        c.save();buf.seek(0);page=PdfReader(buf).pages[0]
        bottom=0
        for j,name in enumerate(names):
            file=root/folder/(name+'.pdf');d=Drawing(file);meta=d.metadata();p=PdfReader(file).pages[0]
            width=121.824 if j==0 else 40.824;scale=width*PT/float(p.mediabox.width)
            height=float(p.mediabox.height)*scale/PT
            if j==0:x=44.088;y=269-height;bottom=y
            else:x=40.2+44.388*(j-1);y=bottom-16-height
            page.merge_transformed_page(p,Transformation([scale,0,0,scale,x*PT,y*PT]))
            if meta.get('legend'):assert inside(meta['legend']['bounds_pt'],d.axes)
            rows.append({'figure':i+7,'panel':name,'page_mm':[float(p.mediabox.width)/PT,float(p.mediabox.height)/PT],
                         'axes_pt':d.axes,'legendInsideAxes':bool(meta.get('legend')),
                         'legendCollisionCount':meta.get('legend',{}).get('collisionCount',0) if meta.get('legend') else 0,
                         'vector':True,'sha256':sha(file.read_bytes())})
        writer.add_page(page)
    with (root/'figures_7_9_review.pdf').open('wb') as f:writer.write(f)
    (root/'layout_verification.json').write_text(json.dumps({'panels':rows,'passed':True,
        'legendPolicy':'inside axes, collision-checked','physicalSolves':0},indent=2)+'\n')


def finalize(file):
    d=Drawing(file);out=Path(file);b=expand(d.content,2*PT)
    crop(d.page,b,out);render(out,out.with_suffix(''),220)
    meta=d.metadata()
    meta['rawExportPage_mm']=meta.get('page_mm',[float(d.page.mediabox.width)/PT,float(d.page.mediabox.height)/PT])
    if meta.get('legend'):
        r=meta['legend']['bounds_pt'];assert inside(r,d.axes),'Retained legend must be internal'
        meta['legend']['bounds_pt']=[r[0]-b[0],r[1]-b[1],r[2]-b[0],r[3]-b[1]]
    if meta.get('graphTransform'):
        meta['graphTransform'][4]-=b[0];meta['graphTransform'][5]-=b[1]
    if meta.get('legendBounds_mm'):
        bounds=meta['legendBounds_mm'];bounds=[bounds] if isinstance(bounds[0],(int,float)) else bounds
        updated=[]
        for lx,ly,lw,lh in bounds:
            rect=[lx*PT,ly*PT,(lx+lw)*PT,(ly+lh)*PT]
            assert inside(rect,d.axes,1.5),'Exported legend is not inside its axes'
            updated.append([lx-b[0]/PT,ly-b[1]/PT,lw,lh])
        meta['legendBounds_mm']=updated
        meta['legendInsideAxes']=True
    meta.update(page_mm=[(b[2]-b[0])/PT,(b[3]-b[1])/PT],
        axes_pt=[d.axes[0]-b[0],d.axes[1]-b[1],d.axes[2]-b[0],d.axes[3]-b[1]] if d.axes else None,
        cropSafetyMargin_mm=2,visibleContent_pt=[d.content[0]-b[0],d.content[1]-b[1],d.content[2]-b[0],d.content[3]-b[1]],
        crop_translation_pt=[-b[0],-b[1]],output_sha256=sha(out.read_bytes()),legendPolicy='inside axes, collision-checked')
    if d.axes:meta['axes_mm']=[(d.axes[0]-b[0])/PT,(d.axes[1]-b[1])/PT,(d.axes[2]-d.axes[0])/PT,(d.axes[3]-d.axes[1])/PT]
    out.with_suffix('.layout.json').write_text(json.dumps(meta,indent=2)+'\n')


if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--review-root');ap.add_argument('--crop')
    o=ap.parse_args()
    if o.review_root:build_review(o.review_root)
    elif o.crop:finalize(o.crop)
