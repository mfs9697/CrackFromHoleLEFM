"""Place immutable retained vector graphics onto common MATLAB axes plates.

No coordinates of physical states are extracted, sampled, or reconstructed.
The signed profiles select whole original PDF drawing commands, preserving
the curves, markers, event markers and legends. Only PDF layout transforms
are added. Source axes decorations are replaced by a blank MATLAB axes plate.
Spatial layers are scaled uniformly, calibrated by the retained circular hole;
comparison layers are resized to the common (non-spatial) plotting rectangle.
"""
import argparse
import copy
import hashlib
import io
import json
import os
import shutil
import subprocess
from pathlib import Path

import pdfplumber
from pypdf import PdfReader, PdfWriter, Transformation
from pypdf.generic import ContentStream, NameObject, DecodedStreamObject
from reportlab.pdfgen import canvas

PAPER = Path(__file__).resolve().parent
SOURCES = PAPER / 'data' / 'retained_figure_7_9'
PT = 72 / 25.4


def sha(data):
    return hashlib.sha256(data).hexdigest()


def layer(reader, span):
    original = reader.pages[0]
    stream = ContentStream(original.get_contents(), reader)
    # All five signed files use the same two outer device transforms.
    assert [op for _, op in stream.operations[:3]] == [b'q', b'cm', b'cm']
    selected = stream.operations[span[0]:span[1]]
    assert not any(op in [b'BT', b'TJ', b'Tj'] for _, op in selected) or span[0] > 1800
    out = copy.copy(original)
    payload = ContentStream(original.get_contents(), reader)
    payload.operations = stream.operations[:3] + selected + [([], b'Q')]
    contents = DecodedStreamObject()
    contents.set_data(payload.get_data())
    out[NameObject('/Contents')] = contents
    return out, sha(repr(selected).encode()), len(selected)


def visible_bounds(page):
    writer = PdfWriter(); writer.add_page(page)
    buffer = io.BytesIO(); writer.write(buffer); buffer.seek(0)
    with pdfplumber.open(buffer) as doc:
        p = doc.pages[0]
        objs = p.chars + p.curves + p.lines + p.rects
        assert objs, 'Empty retained legend layer'
        return [min(o['x0'] for o in objs), p.height - max(o['bottom'] for o in objs),
                max(o['x1'] for o in objs), p.height - min(o['top'] for o in objs)]


def normalize(key, axes_plate, output):
    profile = json.loads((SOURCES / 'profiles.json').read_text())[key]
    layout = json.loads((PAPER / 'figure_7_9_layout.json').read_text())
    spec = layout[profile['class']]
    source = SOURCES / (key + '.pdf')
    assert sha(source.read_bytes()) == profile['sha256'], 'Retained source identity changed'
    reader = PdfReader(source)
    graph, graph_sha, n = layer(reader, profile['graph'])
    dest = PdfReader(axes_plate).pages[0]
    plate_info = json.loads(Path(axes_plate).with_suffix('.layout.json').read_text())
    assert plate_info['xlim'] == ([135,300] if profile['class']=='top' else [0,100])
    assert plate_info['ylim'] == profile['ylim'], 'Axes plate does not match retained units/range'
    expected = [v * PT for v in spec['page']]
    # MATLAB print rounds its MediaBox to whole PDF points. The wrapper sets
    # exact boxes; the inner axes coordinates remain physical MATLAB sizes.
    assert all(abs(float(a) - b) < 1.1 for a, b in zip(dest.mediabox[2:], expected))
    dest.mediabox.upper_right = expected; dest.cropbox.upper_right = expected
    x, y, w, h = [v * PT for v in spec['axes']]
    if profile['class'] == 'top':
        a, b, c, d = profile['circle_bbox_topdown_pt']
        assert abs((c - a) - (d - b)) < 1e-4, 'Source hole is not circular'
        s = (w / 165) / ((c - a) / 60)
        cx = (a + c) / 2; cy = float(reader.pages[0].mediabox.height) - (b + d) / 2
        target_cx = x + w * (170 - 135) / 165
        target_cy = y + h * (-20 + 55) / 70
        transform = [s, 0, 0, s, target_cx - s * cx, target_cy - s * cy]
    else:
        a, b, c, d = profile['axes_bbox_topdown_pt']
        bottom = float(reader.pages[0].mediabox.height) - d
        sx, sy = w / (c - a), h / (d - b)
        transform = [sx, 0, 0, sy, x - sx * a, y - sy * bottom]
    dest.merge_transformed_page(graph, Transformation(transform), over=True)
    legend_info = None
    if profile['legend']:
        from publication_pdf_assets import retained_legend
        original_legend, _, _ = layer(reader, profile['legend'])
        count = 5 if key == 'm1_trajectory' else 3 if key == 'increment_trajectory' else 2
        dest, legend_info = retained_legend(dest, source, visible_bounds(original_legend), count,
            relabel=key.startswith('m1_'), context='Scaled M1' if key=='increment_trajectory' else None)
    writer = PdfWriter(); writer.add_page(dest)
    writer.add_metadata({'/Subject': 'Layout-only wrapping of retained accepted vector evidence; no numerical reconstruction'})
    output = Path(output); output.parent.mkdir(parents=True, exist_ok=True)
    with output.open('wb') as f: writer.write(f)
    report = dict(source=source.relative_to(PAPER.parent).as_posix(), source_sha256=profile['sha256'], panelClass=profile['class'],
                  page_mm=spec['page'], axes_mm=spec['axes'], graphTransform=transform,
                  xlim=[135,300] if profile['class']=='top' else [0,100], ylim=profile['ylim'],
                  originalGraphOperations_sha256=graph_sha, originalGraphOperationCount=n,
                  legend=legend_info, numericalReconstruction=False, rasterization=False,
                  sourceTypography='Original keys and glyphs retained; M0 relabeled in existing font; internal collision checks')
    output.with_suffix('.layout.json').write_text(json.dumps(report, indent=2) + '\n')
    from publication_pdf_assets import finalize
    finalize(output)
    report=json.loads(output.with_suffix('.layout.json').read_text())
    return report


def exact_box(file, kind):
    spec = json.loads((PAPER / 'figure_7_9_layout.json').read_text())[kind]
    target = [v * PT for v in spec['axes']]
    with pdfplumber.open(file) as doc:
        actual = frame_box(doc.pages[0], target)
    sx, sy = target[2] / actual[2], target[3] / actual[3]
    assert abs(sx-1) < .005 and abs(sy-1) < .005, 'Unexpected MATLAB print scaling'
    reader = PdfReader(file); p = reader.pages[0]
    # Correct MATLAB/painters' rounded device viewport, rather than simply
    # declaring a different MediaBox. This makes the physical axes exact too.
    p.add_transformation(Transformation([sx,0,0,sy,target[0]-sx*actual[0],target[1]-sy*actual[1]]))
    p.mediabox.upper_right = [v * PT for v in spec['page']]
    p.cropbox.upper_right = p.mediabox.upper_right
    writer = PdfWriter(); writer.add_page(p)
    with open(file, 'wb') as f: writer.write(f)


def frame_box(page, expected):
    """Find the exported axes frame from its four existing straight edges."""
    horizontal = [e for e in page.edges if e['width'] > expected[2] * .5 and e['height'] < .01]
    vertical = [e for e in page.edges if e['height'] > expected[3] * .5 and e['width'] < .01]
    candidates = []
    for a in horizontal:
        for b in horizontal:
            if abs(a['x0'] - b['x0']) > .01 or abs(a['x1'] - b['x1']) > .01 or a['y0'] >= b['y0']:
                continue
            left = any(abs(v['x0'] - a['x0']) < .01 and abs(v['y0'] - a['y0']) < .01
                       and abs(v['y1'] - b['y0']) < .01 for v in vertical)
            right = any(abs(v['x0'] - a['x1']) < .01 and abs(v['y0'] - a['y0']) < .01
                        and abs(v['y1'] - b['y0']) < .01 for v in vertical)
            if left and right:
                candidates.append([a['x0'], a['y0'], a['width'], b['y0'] - a['y0']])
    assert candidates, 'No axes frame found'
    return min(candidates, key=lambda c: sum(abs(x - y) for x, y in zip(c, expected)))


def review(root):
    from publication_pdf_assets import sensitivity_preview
    sensitivity_preview(root)


if __name__ == '__main__':
    ap = argparse.ArgumentParser()
    ap.add_argument('--key'); ap.add_argument('--plate'); ap.add_argument('--output')
    ap.add_argument('--exact-box'); ap.add_argument('--class', dest='kind', choices=['top', 'lower'])
    ap.add_argument('--review-root'); ap.add_argument('--asset-review-root'); ap.add_argument('--crop')
    opt = ap.parse_args()
    if opt.asset_review_root:
        from publication_pdf_assets import build_review
        build_review(opt.asset_review_root)
        from test_publication_pdf_assets import verify
        verify(opt.asset_review_root)
        from prepare_issue78_review import prepare
        prepare(opt.asset_review_root)
    elif opt.crop:
        from publication_pdf_assets import finalize
        finalize(opt.crop)
    elif opt.review_root: review(opt.review_root)
    elif opt.exact_box: exact_box(opt.exact_box, opt.kind)
    else: normalize(opt.key, opt.plate, opt.output)
