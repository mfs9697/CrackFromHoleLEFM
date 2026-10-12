"""Ensure the committed PDF corresponds to the current manuscript source.

After rebuilding paper/main.pdf:
    python paper/verify_pdf_source_match.py --record
    python paper/verify_pdf_source_match.py

Requires pypdf. No FE solves or figure generation are performed.
The manifest checks byte identity; the textual checks prevent a stale PDF
from being accepted simply because hashes were refreshed.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import re
import unicodedata
from pathlib import Path

from pypdf import PdfReader

PAPER = Path(__file__).resolve().parent
DEFAULT_MANIFEST = PAPER / "reviews/final_typesetting/pdf_source_manifest.json"

# Deliberately anchored to the actual authorized final manuscript wording.
# Do not weaken these checks to make an obsolete PDF pass.
SENTINELS = {
    "abstract_coupled_refinement_scope":
        "coupled refinement reveals trends in mode mixity and "
        "turning density without isolating the effect of the increment",
    "accumulated_length_definition":
        "The accumulated crack-path length measured from the initiation point",
    "accumulated_length_not_horizontal":
        "denotes this length rather than the horizontal coordinate of the crack tip",
    "software_neutral_solver":
        "the PCG iteration may stagnate at a residual level slightly above "
        "the prescribed tolerance because of finite-precision effects",
    "single_restart_qualification":
        "one Krylov restart from the returned iterate is permitted",
    "tip_resolution_scope":
        "The h0 and 2h0 paths share the interaction-domain radii "
        "and the M0 exterior sizing law",
    "exterior_mesh_scope":
        "The M0 and M1 trajectories are propagated separately",
}


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def tokens(text: str) -> str:
    """Normalize TeX/PDF typography, hyphenation and math/text spacing.

    PDF extractors can return "Theh0" for the rendered "The $h_0$".
    Match long, distinctive source/PDF sentinels after collapsing these
    typographic boundaries. Exact file hashes are checked separately.
    """
    text = unicodedata.normalize("NFKC", text).lower()
    # TeX subscript notation h_0 tokenizes in the same way as PDF h0.
    text = text.replace("_", "")
    return "".join(re.findall(r"[a-z0-9]+", text))


def validate(tex_path: Path, pdf_path: Path) -> tuple[dict, list[str]]:
    tex_path, pdf_path = Path(tex_path), Path(pdf_path)
    tex = tex_path.read_text(encoding="utf-8")
    pdf = PdfReader(str(pdf_path))
    rendered = "\n".join(page.extract_text() or "" for page in pdf.pages)
    src, dst = tokens(tex), tokens(rendered)

    failures = []
    for key, snippet in SENTINELS.items():
        test = tokens(snippet)
        if test not in src:
            failures.append(f"{key}: required wording missing from main.tex")
        if test not in dst:
            failures.append(f"{key}: required wording absent from compiled PDF")
    if re.search(r"\bMATLAB\b", tex, flags=re.I):
        failures.append("The manuscript source still mentions MATLAB")
    if re.search(r"\bMATLAB\b", rendered, flags=re.I):
        failures.append("Compiled PDF still mentions MATLAB")
    if pdf.metadata.get("/Author") != "Danylo Selivanov":
        failures.append("Unexpected PDF author metadata")
    if pdf.metadata.get("/Title") != (
        "Incremental LEFM Prediction of a Crack Trajectory from a Circular Hole"
    ):
        failures.append("Unexpected PDF title metadata")

    figure_paths = re.findall(
        r"\\includegraphics\[[^\]]*\]\{([^}]+)\}", tex
    )
    if len(figure_paths) != 31 or len(set(figure_paths)) != 31:
        failures.append(f"Expected 31 distinct approved graphic inclusions, got {len(figure_paths)}")
    files = {}
    for rel in figure_paths:
        path = tex_path.parent / rel
        if not path.is_file():
            failures.append(f"Missing approved graphic: {rel}")
        else:
            files[rel] = sha256(path)

    manifest = {
        "version": 1,
        "main_tex_sha256": sha256(tex_path),
        "references_bib_sha256": sha256(tex_path.parent / "references.bib"),
        "main_pdf_sha256": sha256(pdf_path),
        "pdf_page_count": len(pdf.pages),
        "included_graphic_pdf_sha256": files,
        "sentinel_names": list(SENTINELS),
    }
    return manifest, failures


def check(
    tex_path: Path = PAPER / "main.tex",
    pdf_path: Path = PAPER / "main.pdf",
    manifest_path: Path = DEFAULT_MANIFEST,
    record: bool = False,
) -> dict:
    result, failures = validate(tex_path, pdf_path)
    if failures:
        raise ValueError("PDF/source mismatch:\n- " + "\n- ".join(failures))
    if record:
        manifest_path.parent.mkdir(parents=True, exist_ok=True)
        manifest_path.write_text(
            json.dumps(result, indent=2, ensure_ascii=False) + "\n",
            encoding="utf-8",
        )
    else:
        if not manifest_path.is_file():
            raise ValueError(
                f"Missing PDF/source manifest {manifest_path}; "
                "rebuild the PDF and run --record before merging."
            )
        recorded = json.loads(manifest_path.read_text(encoding="utf-8"))
        if recorded != result:
            changes = sorted(set(recorded) | set(result))
            mismatched = [key for key in changes if recorded.get(key) != result.get(key)]
            raise ValueError(
                "PDF/source manifest differs from files: " + ", ".join(mismatched)
            )
    return result


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--tex", type=Path, default=PAPER / "main.tex")
    parser.add_argument("--pdf", type=Path, default=PAPER / "main.pdf")
    parser.add_argument("--manifest", type=Path, default=DEFAULT_MANIFEST)
    parser.add_argument(
        "--record", action="store_true",
        help="Write manifest only after PDF text and source checks have passed.",
    )
    args = parser.parse_args()
    try:
        result = check(args.tex, args.pdf, args.manifest, args.record)
    except (ValueError, FileNotFoundError) as exc:
        parser.exit(1, f"FAIL: {exc}\n")
    print(
        f"PASS: {result['pdf_page_count']} PDF pages; "
        "current wording and all 31 included graphics bound to source."
    )
