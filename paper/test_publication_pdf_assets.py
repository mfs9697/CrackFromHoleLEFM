"""Independent geometry/source checks for the Issue #78 vector candidates."""
import argparse
import hashlib
import itertools
import json
import tempfile
from pathlib import Path

from pypdf import PdfReader
import publication_pdf_assets as p


def images(page):
    def walk(resources):
        count=0
        for ref in resources.get('/XObject',{}).values():
            obj=ref.get_object()
            if obj.get('/Subtype')=='/Image':count+=1
            elif obj.get('/Subtype')=='/Form':count+=walk(obj.get('/Resources',{}))
        return count
    return walk(page.get('/Resources',{}))


def verify(root):
    root=Path(root);manifest=json.loads((root/'baseline_manifest.json').read_text())
    report=json.loads((root/'geometry_report.json').read_text());assert len(report['panels'])==31
    checks=0;geometry_objects=0
    for path,digest in manifest['protectedData'].items():
        assert p.sha((p.PAPER.parent/path).read_bytes())==digest,path;checks+=1
    for old,record in zip(manifest['assets'],report['panels']):
        asset=old['asset'];assert record['asset']==asset
        source=root/'before'/asset;target=root/'after'/asset
        assert p.sha(source.read_bytes())==old['sha256']
        assert p.sha(target.read_bytes())==record['output_sha256']
        assert images(PdfReader(target).pages[0])==0,'Raster replacement'
        before=p.Drawing(source);after=p.Drawing(target);dx,dy=record['sourceToOutput'][4:]
        assert record['sourceToOutput'][:4]==[1,0,0,1],'Only a rigid crop translation is permitted'
        if before.axes:
            expected=[before.axes[0]+dx,before.axes[1]+dy,before.axes[2]+dx,before.axes[3]+dy]
            assert max(abs(a-b) for a,b in zip(expected,after.axes))<2e-4,asset
        width=float(after.page.mediabox.width);height=float(after.page.mediabox.height)
        assert list(after.page.mediabox)==list(after.page.cropbox)
        assert p.inside(after.content,[0,0,width,height],2e-4)
        actual_margins=[after.content[0]/p.PT,after.content[1]/p.PT,(width-after.content[2])/p.PT,(height-after.content[3])/p.PT]
        assert min(actual_margins)>=1.99 and max(actual_margins)<=3.01,(asset,actual_margins)
        removed=set(record['removedDecorationPaints'])
        if record.get('ylabelAlignment'):
            nodes,_=before.label_nodes();removed.update(i for n in nodes for i in n.paints)
        # Match every unchanged painted primitive, including frames/ticks and
        # labels, under the declared translation. Extra legend glyphs are OK.
        buckets={}
        def color(o):return str(o.get('stroking_color'))+str(o.get('non_stroking_color'))
        for o in after.visible:
            b=p.box(o);key=(o['object_type'],color(o),round(b[0]*1000),round(b[1]*1000))
            buckets.setdefault(key,[]).append(b)
        for o in before.visible:
            if o.get('mcid') in removed:continue
            b=p.box(o);expected=[b[0]+dx,b[1]+dy,b[2]+dx,b[3]+dy]
            qx,qy=round(expected[0]*1000),round(expected[1]*1000);found=False
            for a,c in itertools.product([-1,0,1],repeat=2):
                key=(o['object_type'],color(o),qx+a,qy+c)
                if any(max(abs(x-y) for x,y in zip(expected,bb))<2e-4 for bb in buckets.get(key,[])):
                    found=True;break
            assert found,(asset,'Changed/lost original painted primitive',o.get('mcid'),expected)
            geometry_objects+=1
        legend=record['legend']
        if legend:
            assert legend['insideAxes'] and legend['collisionCount']==0
            assert p.inside(legend['bounds_pt'],after.axes)
            # Repeat the collision calculation from the original unmodified
            # graph primitives, not from a producer's boolean certification.
            old_rect=legend['oldBounds_pt'];obs=p.obstacles(before,removed)
            assert not p.collisions(old_rect,obs),(asset,'Legend crosses original scientific drawing')
        if asset.startswith('figures/m1_mesh_sensitivity/'):
            text=PdfReader(target).pages[0].extract_text()
            assert 'reference' not in text
            if legend:assert 'M0' in text and 'M1' in text
        checks+=1
    # The producer must refuse an altered source before writing candidates.
    with tempfile.TemporaryDirectory() as temp:
        temp=Path(temp);bad=json.loads(json.dumps(manifest));bad['protectedData']={}
        path=temp/'before'/bad['assets'][0]['asset'];path.parent.mkdir(parents=True)
        path.write_bytes(b'%PDF-corrupt source')
        (temp/'baseline_manifest.json').write_text(json.dumps(bad))
        try:p.build_review(temp)
        except AssertionError as e:assert 'Changed original' in str(e)
        else:raise AssertionError('Changed source accepted')
        assert not (temp/'after').exists()
    result={'passed':True,'panelCount':31,'checks':checks,'unchangedOriginalPaintedPrimitives':geometry_objects,
            'sourceMismatchRejected':True,'vectorPanels':31,'legendInsideAxes':True,'legendCollisions':0,
            'marginRange_mm':[2,3],'physicalSolves':0,'dataHistoriesReconstructed':False}
    (root/'verification.json').write_text(json.dumps(result,indent=2)+'\n')
    print(f'PASS: {checks} source/panel checks, {geometry_objects} original vector primitives preserved; 31 vector panels, internal legends, 2--3-mm margins')


if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('root');verify(ap.parse_args().root)
