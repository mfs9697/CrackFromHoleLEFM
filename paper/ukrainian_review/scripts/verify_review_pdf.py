"""Read the final XeLaTeX recorder/log, BibTeX output, PDF links and page text."""

import argparse
import json
from pathlib import Path
import re

import fitz

from verify_preservation import EDITION, hash_bytes

JOB = "AAM_UA_Bogdanov_review"
PAGE_MARKERS = (
    "0.011415", "0.36648", "84.14941", "214.931", "83.268112", "0.06734",
)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("build_dir", type=Path, help="Directory containing the final XeLaTeX outputs")
    parser.add_argument("--write-audit", action="store_true")
    args = parser.parse_args()
    manifest = json.loads((EDITION / "SOURCE_MANIFEST.json").read_text(encoding="utf-8"))
    paths = {i["path"].removeprefix("paper/") for i in manifest["graphic_assets"]}
    pdf = args.build_dir / f"{JOB}.pdf"
    log = (args.build_dir / f"{JOB}.log").read_text(encoding="utf-8")
    bbl = (args.build_dir / f"{JOB}.bbl").read_text(encoding="utf-8")
    aux = (args.build_dir / f"{JOB}.aux").read_text(encoding="utf-8")
    blg = (args.build_dir / f"{JOB}.blg").read_text(encoding="utf-8")
    fls = (args.build_dir / f"{JOB}.fls").read_text(encoding="utf-8").splitlines()
    bibkeys = re.findall(r"\\bibcite\{([^}]+)\}", aux)
    graphic_inputs = {}
    working_dir = Path(next(line[4:] for line in fls if line.startswith("PWD ")))
    for line in fls:
        if line.startswith("INPUT ") and "figures/" in line and line.endswith(".pdf"):
            path = (working_dir / line[6:]).resolve()
            graphic_inputs[path.relative_to(EDITION).as_posix()] = hash_bytes(path.read_bytes())
    packaged_figures = {p: hash_bytes((EDITION / p).read_bytes()) for p in paths}
    expected_dois = {x for x in manifest["cited_bibliography"].values() if x}
    with fitz.open(pdf) as document:
        uris, internal, page_text, bad_pages = set(), 0, [], []
        names = document.resolve_names()
        broken_destinations = []
        for page in document:
            page_text.append(page.get_text())
            for link in page.get_links():
                if "uri" in link:
                    uris.add(link["uri"])
                elif link["kind"] in (fitz.LINK_GOTO, fitz.LINK_NAMED):
                    named = link.get("nameddest")
                    if not 0 <= link.get("page", -1) < len(document) or (named and named not in names):
                        broken_destinations.append(named or link.get("page"))
                    internal += 1
            if abs(page.rect.width - 595.28) > 0.05 or abs(page.rect.height - 841.89) > 0.05:
                bad_pages.append(page.number + 1)
        text = "\n".join(page_text)
        checks = {
            "main_source_recorded": any(line.startswith("INPUT ") and
                (working_dir / line[6:]).resolve() == (EDITION / "main_ua.tex").resolve() for line in fls),
            "compiled_all_localized_figure_paths_and_hashes": graphic_inputs == packaged_figures,
            "reference_keys_in_original_citation_order": bibkeys == list(manifest["cited_bibliography"]),
            "bibliography_entry_count": len(re.findall(r"\\bibitem\[", bbl)) == len(bibkeys),
            "all_doi_urls_hyperlinked": expected_dois <= uris,
            "internal_links_resolve": internal > 0 and not broken_destinations,
            "a4_page_dimensions": not bad_pages,
            "at_least_one_unique_numeric_marker_per_table": all(marker in text for marker in PAGE_MARKERS),
            "ukrainian_front_matter_and_appendix": all(s in text for s in
                ("Анотація", "Ключові слова.", "Додаткові показники чисельної прийнятності", "Література")),
            "all_main_figures_and_tables_rendered": all(f"Рис. {i}:" in text for i in range(1, 11)) and
                all(f"Табл. {i}:" in text for i in range(1, 7)),
            "sequential_equation_numbers": all(re.search(r"(?m)^\(" + str(i) + r"\)$", text)
                for i in range(1, 29)),
            "no_missing_reference_font_or_hyphenation_warnings": not re.search(
                r"undefined (?:references|citations)|multiply defined|Missing character|"
                r"No hyphenation patterns|font .* not found|mapping .* not found", log, re.I),
            "no_overfull_boxes": not re.search(r"Overfull \\[hv]box", log),
            "bibtex_without_warnings_or_errors": not re.search(r"(?im)^(?:warning--|error[: -]|\(there (?:was|were))", blg),
        }
        result = {
            "programme": "XeLaTeX → BibTeX → XeLaTeX → XeLaTeX",
            "engine": log.splitlines()[0],
            "paper_size": "A4, 12 pt, Noto Serif, Ukrainian hyphenation",
            "pdf_sha256": hash_bytes(pdf.read_bytes()), "pages": len(document),
            "source_sha256": hash_bytes((EDITION / "main_ua.tex").read_bytes()),
            "biblio_doi_links": len(expected_dois & uris),
            "pdf_internal_links": internal,
            "checks": checks,
        }
    print(json.dumps(result, ensure_ascii=False, indent=2))
    if args.write_audit:
        (EDITION / "audit/BUILD_VERIFICATION.json").write_text(
            json.dumps(result, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    if not all(checks.values()):
        raise SystemExit(1)


if __name__ == "__main__":
    main()
