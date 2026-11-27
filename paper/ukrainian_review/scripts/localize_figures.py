"""Replace only label glyphs in the frozen vector PDFs. No plotting or FEM."""

import argparse
from collections import Counter
import hashlib
from io import BytesIO
import json
from pathlib import Path
import re

import fitz
from pypdf import PdfReader

EDITION = Path(__file__).resolve().parents[1]
# Rectangles are in unscaled source-PDF points (top-left origin).
# Unlisted labels contain only mathematical symbols, SI units, or abbreviations.
# Each record supplies source wording, replacement, bounding region, and angle.
LABELS = {
    1: [
        ("Free", "Вільний край", (45, 81, 70, 94), 0),
        ("Free", "Вільний край", (243, 81, 270, 94), 0),
    ],
    5: [
        ("Accepted", "Прийнята", (245, 17, 295.7, 30.5), 0),
        ("Qualified, unsolved extension", "Перевірене продовження без розв’язку", (245, 30, 401, 43), 0),
    ],
    6: [("unsolved", "без розв’язку", (312, 54, 358, 66), 0)],
    9: [("[deg]", "[град]", (4, 70, 18, 94), 90)],
    11: [("[deg]", "[град]", (4, 61, 18, 84), 90)],
    12: [("COD–EDI turn [deg]", "Кут COD–EDI [град]", (4, 12, 18, 116), 90)],
    13: [("COD–EDI turn [deg]", "Кут COD–EDI [град]", (4, 12, 18, 116), 90)],
    14: [
        ("mode mixity", "змішаність мод", (4, 8.4, 20, 69.5), 90),
        ("Crack length", "Довжина тріщини", (70, 150.5, 136, 164), 0),
    ],
    15: [
        ("mode mixity", "змішаність мод", (4, 8.4, 20, 69.5), 90),
        ("Crack length", "Довжина тріщини", (70, 150.5, 136, 164), 0),
    ],
    16: [
        ("M0 exterior", "Зовнішня сітка M0", (267, 13, 319, 29), 0),
        ("local symmetry: 84.146 mm", "локальна симетрія: 84.146 mm", (282, 60.0, 412, 71.8), 0),
        ("last accepted", "останній прийнятий", (299, 72, 412, 83), 0),
    ],
    18: [("[mdeg]", "[мград]", (4, 35, 18, 68), 90)],
    19: [("M0 exterior", "Зовнішня M0", (51, 97.5, 103, 113), 0)],
    20: [
        ("M0 q max: P17", "M0, макс. q: P17", (291, 38, 360, 54), 0),
        ("local symmetry: 84.146 mm", "локальна симетрія: 84.146 mm", (290, 53, 413, 68), 0),
        ("P23: last accepted", "P23: останній прийнятий", (290, 65, 413, 81), 0),
    ],
    22: [("[mdeg]", "[мград]", (4, 34, 18, 67), 90)],
    24: [("Scaled M1", "Масштабована M1", (311, 13, 359, 29), 0)],
    25: [("Successive", "Послідовні", (4, 67.5, 18, 120), 90)],
    26: [("Direction difference [mdeg]", "Різниця напрямів [мград]", (4, 9, 18, 140), 90)],
    27: [
        ("Scaled M1", "Масшт. M1", (51, 13, 99, 29), 0),
        ("deg", "град", (4, 56, 18, 76), 90),
    ],
    28: [("Minimum angle [deg]", "Мінімальний кут [град]", (4, 13, 18, 118), 90)],
    29: [("Adjacent size ratio", "Відн. розмірів сусідніх ел.", (4, 17, 18, 112), 90)],
    30: [
        ("PCG iterations", "Ітерації PCG", (4, 26, 18, 104), 90),
        ("Crack length", "Довжина тріщини", (80, 139, 146, 153), 0),
    ],
    31: [
        ("True residual / gate", "Фактична нев’язка / поріг", (4, 14, 18, 115), 90),
        ("Crack length", "Довжина тріщини", (80, 139, 146, 153), 0),
    ],
}


def digest(data):
    return hashlib.sha256(data).hexdigest()


def outline_key(path):
    """Compare visible vector geometry independent of PDF object IDs."""
    def normalize(obj):
        if isinstance(obj, float):
            return round(obj, 4)
        if isinstance(obj, dict):
            return {k: normalize(v) for k, v in obj.items() if k not in ("seqno", "layer")}
        if isinstance(obj, (list, tuple, fitz.Point, fitz.Rect, fitz.Quad)):
            return [normalize(v) for v in obj]
        return obj

    return json.dumps(normalize(path), sort_keys=True)


def union(rects):
    rect = fitz.Rect()
    for part in rects:
        rect |= part
    return rect


def glyphs_in_region(drawings, spans, box):
    outlines = [path for path in drawings if path["type"] == "f" and path["rect"] in box]
    native = [span for span in spans if (fitz.Rect(span["bbox"]) & box).get_area() /
              fitz.Rect(span["bbox"]).get_area() >= 0.90]
    return outlines, native


def remove_label_paint(doc, page, source, boxes, native_boxes):
    """Suppress only the paint operator for each bounded fill path."""
    reader = PdfReader(BytesIO(source))
    content = reader.pages[0].get_contents()
    stack = []
    matrix = fitz.Matrix(1, 0, 0, 1, 0, 0)
    points = []
    clipping = False
    new_operations = []
    traces = iter(page.get_texttrace())
    removed_text = 0
    painting = {b"S", b"s", b"f", b"F", b"f*", b"B", b"B*", b"b", b"b*", b"n"}
    for operands, operator in content.operations:
        if operator == b"q":
            stack.append(matrix)
        elif operator == b"Q":
            matrix = stack.pop()
        elif operator == b"cm":
            matrix = fitz.Matrix(*(float(x) for x in operands)) * matrix
        elif operator in (b"m", b"l", b"c", b"v", b"y"):
            points.extend((float(x), float(y)) for x, y in zip(operands[::2], operands[1::2]))
        elif operator == b"re":
            x, y, w, h = map(float, operands)
            points.extend(((x, y), (x+w, y), (x+w, y+h), (x, y+h)))
        elif operator in (b"W", b"W*"):
            clipping = True
        elif operator in painting:
            if operator in (b"f", b"F", b"f*") and points and not clipping:
                transform = matrix * page.transformation_matrix
                mapped = [fitz.Point(*pt) * transform for pt in points]
                x, y = zip(*mapped)
                bounds = fitz.Rect(min(x), min(y), max(x), max(y))
                if any(bounds in box for box in boxes):
                    operator = b"n"
            points = []
            clipping = False
        if operator in (b"Tj", b"TJ") and any(operands):
            trace = next(traces)
            if any(fitz.Rect(trace["bbox"]) in box for box in native_boxes):
                removed_text += 1
                continue
        new_operations.append((operands, operator))
    assert next(traces, None) is None
    assert removed_text == len(native_boxes)
    content.operations = new_operations
    xref = doc.get_new_xref()
    doc.update_object(xref, "<<>>")
    doc.update_stream(xref, content.get_data())
    doc.xref_set_key(page.xref, "Contents", f"{xref} 0 R")


def draw_label(page, font, bounds, box, label, angle, color, figure):
    if angle == 90:
        room = box.height
        size = min(10.8, room / font.text_length(label, fontsize=1))
        length = font.text_length(label, fontsize=size)
        point = (bounds.x0 + 0.78 * bounds.width,
                 (box.y0 + box.y1 + length) / 2)
    else:
        size = min(10.8, 0.96 * box.width / font.text_length(label, fontsize=1))
        if figure in (1, 16, 19, 24, 27):
            size = 9.0
        if figure == 1:
            size = 7.0
        if figure in (19, 27):
            size = min(size, 8.4)
        if figure in (14, 15, 30, 31):
            size = min(size, 9.3)
        x = bounds.x0
        if figure == 1:
            x -= 3 if x < 100 else 19
        point = (x, (bounds.y0 + bounds.y1) / 2 + size * 0.32)
    page.insert_text(point, label, fontname="UkrainianLabel", fontsize=size,
                     color=color, rotate=angle, overlay=True)
    return {"font_size_pt": round(size, 3), "text_origin_pt": list(point)}


def localize(source, output, changes, font_path, figure):
    original = source.read_bytes()
    with fitz.open(stream=original, filetype="pdf") as doc:
        page = doc[0]
        original_boxes = (page.mediabox, page.cropbox)
        original_images = len(page.get_images())
        old_drawings = page.get_drawings()
        spans = [s for block in page.get_text("dict")["blocks"]
                 for line in block.get("lines", ()) for s in line["spans"]]
        removed = Counter()
        records = []
        for old, new, region, angle in changes:
            box = fitz.Rect(region)
            outlines, native = glyphs_in_region(old_drawings, spans, box)
            if not outlines and not native:
                raise ValueError(f"{source.name}: no glyphs found for {old}")
            removed.update(outline_key(path) for path in outlines)
            bounds = union([p["rect"] for p in outlines] + [fitz.Rect(s["bbox"]) for s in native])
            color = outlines[0]["fill"] if outlines else (0, 0, 0)
            # Remove outlined glyphs individually. Keep the background and images.
            records.append((old, new, angle, box, bounds, color, outlines, native))

        expected = Counter(outline_key(path) for path in old_drawings) - removed
        info = []
        if records:
            outline_boxes = [r[3] for r in records if r[6]]
            native_boxes = [fitz.Rect(span["bbox"]) + (-0.01, -0.01, 0.01, 0.01)
                            for record in records for span in record[7]]
            remove_label_paint(doc, page, original, outline_boxes, native_boxes)
            page = doc.reload_page(page)
            actual = Counter(outline_key(path) for path in page.get_drawings())
            if actual != expected:
                raise ValueError(f"{figure} {source.name}: unexpected vector changes "
                                 f"({sum((expected-actual).values())} missing, "
                                 f"{sum((actual-expected).values())} extra, "
                                 f"{sum(removed.values())} selected)")
            font = fitz.Font(fontfile=str(font_path))
            page.insert_font(fontname="UkrainianLabel", fontfile=str(font_path))
            for old, new, angle, box, bounds, color, paths, texts in records:
                position = draw_label(page, font, bounds, box, new, angle, color, figure)
                info.append(dict(original=old, ukrainian=new, source_region_pt=list(box),
                                 outlined_glyphs=len(paths), native_spans=len(texts), angle=angle, **position))
            doc.subset_fonts()
            doc.set_metadata({})
            output.write_bytes(doc.tobytes(garbage=4, deflate=True))
        else:
            output.write_bytes(original)

        with fitz.open(output) as built:
            copied = built[0]
            assert (copied.mediabox, copied.cropbox) == original_boxes
            assert Counter(outline_key(path) for path in copied.get_drawings()) == expected
            assert len(copied.get_images()) == original_images
            text = copied.get_text()
            for entry in info:
                assert entry["ukrainian"] in text, (figure, entry["ukrainian"], text)
        return {
            "sha256_english": digest(original), "sha256_ukrainian": digest(output.read_bytes()),
            "vector_paths_preserved": sum(expected.values()),
            "vector_geometry_check": "PASS",
            "status": "localized" if records else "symbolic-labels-and-international-units",
            "replacements": info,
        }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--paper-dir", type=Path, default=EDITION.parent)
    parser.add_argument("--font", type=Path, required=True, help="Noto Serif Regular TTF")
    args = parser.parse_args()
    manuscript = (args.paper_dir / "submissions/aam/manuscript.tex").read_text(encoding="utf-8")
    paths = re.findall(r"\\includegraphics\[[^]]*\]\{([^}]+)\}", manuscript)
    audit = {}
    for i, path in enumerate(paths, 1):
        output = EDITION / path
        output.parent.mkdir(parents=True, exist_ok=True)
        audit[path] = localize(args.paper_dir / path, output, LABELS.get(i, ()), args.font, i)
        print(f"{i:2} {audit[path]['status']:39} {path}", flush=True)
    report = {
        "method": "Ukrainian text replaces glyphs in the existing vector exports. "
                  "Original scientific curves, markers, meshes, ticks, axes, graphics boxes "
                  "and all non-label drawing paths are preserved.",
        "font_sha256": digest(args.font.read_bytes()),
        "pymupdf": fitz.VersionBind,
        "assets": audit,
    }
    (EDITION / "audit/FIGURE_PROVENANCE.json").write_text(
        json.dumps(report, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")


if __name__ == "__main__":
    main()
