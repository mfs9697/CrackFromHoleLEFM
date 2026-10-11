"""Lightweight layout/provenance tests; uses existing generated review PDFs."""
import argparse
import json
import shutil
import tempfile
from pathlib import Path

from pypdf import PdfReader
import normalize_retained_sensitivity_pdf as tool


def test(review_root):
    root = Path(review_root)
    profiles = json.loads((tool.SOURCES/'profiles.json').read_text())
    files = dict(m1_trajectory='figure8_M1_exterior/trajectory',
                 m1_vertical_deviation='figure8_M1_exterior/vertical_deviation',
                 m1_direction_deviation='figure8_M1_exterior/direction_deviation',
                 m1_mode_mixity='figure8_M1_exterior/mode_mixity',
                 increment_trajectory='figure9_increment/trajectory')
    for key, path in files.items():
        p = profiles[key]
        source = tool.SOURCES/(key+'.pdf')
        assert tool.sha(source.read_bytes()) == p['sha256']
        reader = PdfReader(source)
        graph, digest, count = tool.layer(reader,p['graph'])
        meta = json.loads((root/(path+'.layout.json')).read_text())
        assert meta['source'] == f'paper/data/retained_figure_7_9/{key}.pdf', 'Retained source must be repository-relative'
        assert meta['originalGraphOperations_sha256'] == digest
        assert meta['originalGraphOperationCount'] == count
        assert meta['numericalReconstruction'] is False and meta['rasterization'] is False
        if p['class'] == 'top':
            a,b,c,d,e,f = meta['graphTransform']
            assert a == d and b == c == 0, 'Retained spatial layer must have uniform scaling'
        if meta['legend']:
            spec=json.loads((tool.PAPER/'figure_7_9_layout.json').read_text())[p['class']]
            assert meta['legend']['bounds_pt'][1] > (spec['axes'][1]+spec['axes'][3])*tool.PT
        # The full selected source vector stream is a contiguous subsequence
        # of the wrapped output (resource names may be renamed by PDF merge).
        ops = tool.ContentStream(graph.get_contents(), reader).operations[3:-1]
        actual_reader=PdfReader(root/(path+'.pdf'))
        actual = tool.ContentStream(actual_reader.pages[0].get_contents(),actual_reader).operations
        # PDF float serialization can round tiny 8th-decimal device values;
        # compare normalized numerical operands rather than their repr.
        def signature(items):
            def norm(v):
                if isinstance(v,(int,float)):return round(float(v),5)
                if isinstance(v,(list,tuple)):return tuple(norm(x) for x in v)
                return str(v)
            return [(norm(args),op) for args,op in items if op != b'gs']
        wanted=signature(ops);found=signature(actual)
        starts=[i for i,x in enumerate(found) if x==wanted[0]]
        assert any(found[i:i+len(wanted)]==wanted for i in starts), 'Retained drawing operands changed'
    # A missing or changed retained source must fail before a plate/output is used.
    original=tool.SOURCES
    with tempfile.TemporaryDirectory() as temp:
        tool.SOURCES=Path(temp)
        shutil.copy(original/'profiles.json',tool.SOURCES/'profiles.json')
        (tool.SOURCES/'m1_trajectory.pdf').write_bytes(b'%PDF-Changed-source')
        try:
            tool.normalize('m1_trajectory','missing-plate.pdf',Path(temp)/'forbidden.pdf')
        except AssertionError as e:
            assert 'identity' in str(e)
        else: raise AssertionError('Changed source accepted')
        assert not (Path(temp)/'forbidden.pdf').exists()
    tool.SOURCES=original
    audit=json.loads((root/'layout_verification.json').read_text())
    assert audit['passed'] and len(audit['panels'])==12
    assert len(PdfReader(root/'figures_7_9_review.pdf').pages)==3
    print('PASS: 12 panel geometries; five retained vector streams unchanged; source mismatch rejected; three-page review')


if __name__ == '__main__':
    ap=argparse.ArgumentParser();ap.add_argument('review_root');test(ap.parse_args().review_root)
