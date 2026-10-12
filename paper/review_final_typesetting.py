"""Create a vector three-column comparison proof and warning/provenance record.

Requires pypdf, reportlab, Pillow. No plots, numerical calculations, or PDF assets
are changed. Pass the author's supplied PDF explicitly; never use it as a figure
source. The middle column is the committed Issue #78 reviewed manuscript.
"""
from pathlib import Path
import argparse
import hashlib
import io
import json
import re
from pypdf import PdfReader, PdfWriter, Transformation
from reportlab.pdfgen import canvas
from PIL import Image, ImageDraw

PAPER = Path(__file__).resolve().parent


def final_pass_warnings(path):
    record = json.loads(path.read_text(encoding='utf-8-sig'))
    assert record['exitCode'] == 0 and record['pdfExists']
    log = record['log']
    # Tectonic performs its own bibliography/reference reruns. Report convergence
    # warnings separately from expected first-pass unresolved labels/citations.
    assert '(main.tex' in log, 'Expected a main.tex compilation log'
    final = log[log.rfind('(main.tex'):]
    assert 'Output written on main.xdv' in final, 'Incomplete final TeX pass'
    return {
        'compiler': record['compiler'],
        'final_pass_underfull_boxes': re.findall(r'^Underfull .*$', final, re.M),
        'final_pass_overfull_boxes': re.findall(r'^Overfull .*$', final, re.M),
        'final_pass_latex_warnings': re.findall(r'^LaTeX Warning:.*$', final, re.M),
        'transient_undefined_citation_warnings': len(re.findall(r"LaTeX Warning: Citation .*?undefined", log, re.S)),
    }


def review(supplied, before_log, after_log, render_dir=None):
    folder = PAPER / 'reviews/final_typesetting'
    folder.mkdir(parents=True, exist_ok=True)
    paths = [supplied, PAPER / 'reviews/issue78/compiled/manuscript_after.pdf', PAPER / 'main.pdf']
    readers = [PdfReader(p) for p in paths]
    rows = [
        ('Introduction compound-word overflow', 2, 2, 2),
        ('Figure 2 and interaction-integral flow', 7, 7, 7),
        ('Material-domain notation and adjacent definitions', 8, 8, 8),
        ('Characteristic states / late-response figure', 13, 13, 13),
        ('COD figure flow', 14, 14, 14),
        ('Tip-resolution discussion and table', 15, 15, 15),
        ('Figure 7 captions and complete exterior-mesh discussion', 16, 16, 16),
        ('Figure 8 approved M0/M1 asset and captions', 17, 17, 17),
        ('Figure 9 captions and definitions', 18, 18, 18),
        ('Section 5.7 ending', 19, 19, 19),
        ('Discussion begins without a detached Section 5.7 sentence', 20, 20, 20),
        ('Conclusion / appendix ordering', 21, 21, 21),
        ('Figure S1 precedes References', 24, 24, 22),
        ('References conclude the article', 23, 23, 23),
    ]
    writer = PdfWriter()
    width, height = 1825, 912
    column_names = ['Author-supplied PDF (local stale Figure 8d)', 'Issue #78 approved review', 'Final typesetting / PR #80 graphics']
    for title, *page_numbers in rows:
        stream = io.BytesIO()
        c = canvas.Canvas(stream, pagesize=(width, height))
        c.setFont('Helvetica-Bold', 15)
        c.drawString(15, height-22, title)
        c.setFont('Helvetica', 12)
        for col, number in enumerate(page_numbers):
            c.drawString(15+col*605, height-45, f'{column_names[col]} - page {number}')
        c.setFont('Helvetica', 10)
        c.drawString(15, 12, 'Native vector manuscript pages; comparison only. No plot geometry or numerical evidence edited.')
        c.save()
        page = writer.add_page(PdfReader(stream).pages[0])
        for col, number in enumerate(page_numbers):
            source = readers[col].pages[number-1]
            page.merge_transformed_page(source, Transformation().translate(15+col*605, 28))
    writer.add_metadata({'/Title': 'Final typesetting: supplied / Issue 78 / corrected comparisons', '/Author': 'Danylo Selivanov'})
    proof = folder / 'affected_pages_comparison.pdf'
    with proof.open('wb') as f:
        writer.write(f)
    manifest = {
        'baseline_commit': '799f2e190fa90dd41b7fee178763d66ffcaada29',
        'inputs': [{'role': name, 'sha256': hashlib.sha256(p.read_bytes()).hexdigest(), 'pages': len(r.pages)}
                   for name, p, r in zip(column_names, paths, readers)],
        'comparison_rows': [{'defect': row[0], 'supplied_page': row[1], 'reviewed_page': row[2], 'final_page': row[3]} for row in rows],
        'baseline_warnings': final_pass_warnings(before_log),
        'final_warnings': final_pass_warnings(after_log),
        'comparison_sha256': hashlib.sha256(proof.read_bytes()).hexdigest(),
        'physical_calculations': 0,
    }
    (folder / 'review_manifest.json').write_text(json.dumps(manifest, indent=2) + '\n')
    if render_dir:
        pages = sorted(render_dir.glob('page-*.png'))
        assert len(pages) == len(readers[2].pages), 'Render every final page before making contact sheets'
        for k in range(0, len(pages), 4):
            sheet = Image.new('RGB', (1860, 2720), '#e8e8e8')
            draw = ImageDraw.Draw(sheet)
            for j, path in enumerate(pages[k:k+4]):
                image = Image.open(path).convert('RGB')
                image.thumbnail((900, 1310))
                x, y = 20+(j % 2)*930, 40+(j // 2)*1360
                sheet.paste(image, (x, y))
                draw.text((x, y-22), path.name, fill='black')
            sheet.save(render_dir / f'contact-{k//4+1}.png')
    return manifest


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--supplied', type=Path, required=True)
    parser.add_argument('--before-log', type=Path, required=True)
    parser.add_argument('--after-log', type=Path, required=True)
    parser.add_argument('--render-dir', type=Path)
    args = parser.parse_args()
    result = review(args.supplied, args.before_log, args.after_log, args.render_dir)
    print(json.dumps({'pages': result['inputs'][2]['pages'], 'baseline': result['baseline_warnings'], 'final': result['final_warnings']}, indent=2))
