"""Compare the Ukrainian edition against the pinned AAM Git objects."""

import argparse
import hashlib
import json
from pathlib import Path
import re
import subprocess

EDITION = Path(__file__).resolve().parents[1]
REPO = EDITION.parents[1]
BASE = "7597cba1a378e6f436a3bfab304f85702902ee4e"
EN = "paper/submissions/aam/manuscript.tex"
EN_BLOB = "3f418857f7d54b351b1d62e4a003c199c3c7bc9d"
BIB = "paper/references.bib"
NUMBER = re.compile(r"(?<![\w.])[+-]?\d+(?:\.\d+)?(?!\w|\.\d)")
DISPLAY = re.compile(
    r"\\begin\{(equation\*?|align\*?)\}(.*?)\\end\{\1\}|\\\[(.*?)\\\]", re.S)
MATH = re.compile(
    r"\\begin\{(?:equation\*?|align\*?)\}.*?\\end\{(?:equation\*?|align\*?)\}"
    r"|\\\[.*?\\\]|(?<!\\)\$(?:\\.|[^$])*?\$", re.S)
TOKEN = re.compile(r"\\[A-Za-z]+|\\[^\sA-Za-z]|\d+(?:\.\d+)?|[^\s]")
GRAPHIC = re.compile(r"\\includegraphics\[[^\]]*\]\{([^}]+)\}")
REFERENCE = re.compile(r"\\(eqref|ref|cite|label)\{([^}]+)\}")
METADATA = {
    "title": "Incremental LEFM Prediction of a Crack Trajectory from a Circular Hole",
    "author": "Danylo Selivanov",
    "email": "fract@inmech.kyiv.ua",
    "affiliation": "S.P. Timoshenko Institute of Mechanics, National Academy of Sciences of Ukraine, Kyiv, Ukraine",
}
FIGURE_PROVENANCE = (
    "paper/ISSUE78_FIGURE_PIPELINE.md",
    "paper/FINAL_FIGURE_WORKFLOW.md",
    "paper/FIGURE_CONSISTENCY_AUDIT.md",
    "paper/FIGURES_7_9_LAYOUT_REPORT.md",
    "paper/EVIDENCE_AND_PROVENANCE.md",
    "paper/FINAL_TYPESETTING_REVIEW.md",
)
# Explicit localizations inside math (English form: Ukrainian form). Everything else -- operators, signs,
# indices and numbers -- must remain literally identical. AAM-UA7 (WP5, findings T-01..T-06, T-14) adds the
# Cyrillic unit symbols, one-to-one, and the sentence-final period of one display after which the Ukrainian
# text starts a new sentence (the English continues with "and with ..."). Negative controls below check that
# a different unit is still detected.
MATH_TEXT = {
    r"\text{in }": r"\text{в }",
    r"\mathrm{mm}": r"\text{мм}",
    r"\mu\mathrm{m}": r"\text{мкм}",
    r"\mathrm{MPa}": r"\text{МПа}",
    r"\mathrm m": r"\text{м}",
    r"[$\mu$m]": "[мкм]",   # table headers (T-05); the unit itself is checked by table_units
    "\\frac{L_{\\mathrm{tr}}}{\\Delta a}=1.00,\n\\]": "\\frac{L_{\\mathrm{tr}}}{\\Delta a}=1.00.\n\\]",
}
UNIT_TEXT = {
    "mm": "mm", "мм": "mm", "deg": "deg", "град": "deg", "mdeg": "mdeg",
    "мград": "mdeg", "m": "m", "м": "m", "GPa": "GPa", "ГПа": "GPa",
    "MPa": "MPa", "МПа": "MPa", "мкм": "um",
}


def git(*args):
    return subprocess.check_output(["git", "-C", str(REPO), *args])


def blob(path):
    return git("rev-parse", f"{BASE}:{path}").decode().strip()


def pinned(path):
    return git("show", f"{BASE}:{path}")


def hash_bytes(data):
    return hashlib.sha256(data).hexdigest()


def identity(path):
    return {"path": path, "git_blob_sha": blob(path), "sha256": hash_bytes(pinned(path))}


def comments_removed(text):
    return re.sub(r"(?<!\\)%[^\n]*", "", text)


def numerical_macros(text):
    return text.split("% BEGIN AUDITED DATA:", 1)[1].split("% END AUDITED DATA", 1)[0]


def scientific_body(text):
    return re.split(r"\\section\{", text, maxsplit=1)[1].split(r"\bibliography", 1)[0]


def math_tokens(text):
    for en, ua in MATH_TEXT.items():
        text = text.replace(ua, en)
    return [TOKEN.findall(item.group()) for item in MATH.finditer(text)]


def table_sources(text):
    return list(re.finditer(r"\\begin\{table\}.*?\\end\{table\}", text, re.S))


def plain_numbers(text):
    # Ranges have a typesetting dash. Do not treat the second endpoint as negative.
    text = re.sub(r"(\d)--(?=\d)", r"\1–", text)
    text = re.sub(r"(\d)-(?=\d)", r"\1–", text)
    return NUMBER.findall(text)


def table_cells(text):
    rows = []
    for table in table_sources(text):
        matched = re.search(r"\\midrule\s*(.*?)\\bottomrule", table.group(), re.S)
        content = matched[1].replace(r"\CharacteristicRows", re.search(
            r"\\newcommand\{\\CharacteristicRows\}\{(.*?)(?=\n% END AUDITED DATA)",
            text, re.S)[1][:-1])
        records = [row.strip() for row in content.split(r"\\") if row.strip()]
        rows.append([[TOKEN.findall(f) for f in MATH.findall(cell)] +
                     plain_numbers(MATH.sub("", cell)) +
                     re.findall(r"\\(?:LSLengthMM)\b", cell)
                     for row in records for cell in row.split("&")])
    return rows


def units_in_tables(text):
    headers = [t.group().split(r"\midrule")[0].split(r"\toprule")[-1]
               for t in table_sources(text)]
    return [re.findall(r"\[(mm|мм|deg|град|mdeg|мград|мкм|\$\\mu\$[mм])(?=[\]$\/])", h)
            for h in headers]


def normalized_units_in_tables(text):
    # micrometre: English "$\mu$m" and the Ukrainian symbol "мкм" (AAM-UA7, T-05) are the same unit "um"
    return [["um" if x.startswith("$") else UNIT_TEXT[x] for x in t]
            for t in units_in_tables(text)]


def section_shape(text):
    return re.findall(r"\\(section\*?|subsection|appendix)\b(\*)?", text)


def compare(en, ua):
    en_body, ua_body = [scientific_body(comments_removed(t)) for t in (en, ua)]
    en_no_layout, ua_no_layout = [re.sub(
        r"\\begin\{(?:minipage|subfigure)\}(?:\[t\])?\{[^}]+\}|"
        r"\\(?:setlength|Needspace)\{[^}]+\}(?:\{[^}]+\})?|"
        r"\\includegraphics\[[^\]]*\]\{[^}]+\}|\\begin\{tabular\*?\}(?:\{[^}]+\}){1,2}|"
        r"\\begin\{figure\}\[[^\]]+\]|(?<=панелі~)\([bcd]\)",
        "", body) for body in (en_body, ua_body)]
    checks = {
        "audited_macro_definitions": numerical_macros(en) == numerical_macros(ua),
        "all_math_tokens": math_tokens(en_body) == math_tokens(ua_body),
        "display_math_blocks": len(DISPLAY.findall(en_body)) == len(DISPLAY.findall(ua_body)),
        "table_numerical_cells_and_signs": table_cells(en) == table_cells(ua),
        "table_units": normalized_units_in_tables(en) == normalized_units_in_tables(ua),
        "all_labels_and_references": REFERENCE.findall(en_body) == REFERENCE.findall(ua_body),
        "citation_groups": re.findall(r"\\cite\{([^}]+)\}", en_body) == re.findall(r"\\cite\{([^}]+)\}", ua_body),
        "graphics_identity_and_order": GRAPHIC.findall(en_body) == GRAPHIC.findall(ua_body),
        "table_and_figure_environments": re.findall(r"\\begin\{(table|figure|subfigure)\}", en_body) ==
            re.findall(r"\\begin\{(table|figure|subfigure)\}", ua_body),
        "section_hierarchy": section_shape(en) == section_shape(ua),
        "appendix_structure": REFERENCE.findall(en_body.split(r"\appendix")[1]) ==
            REFERENCE.findall(ua_body.split(r"\appendix")[1]),
        "explicit_numerical_lexemes": plain_numbers(en_no_layout) == plain_numbers(ua_no_layout),
    }
    en_abstract = re.search(r"\\abstract\{(.*?)\n\}", en, re.S)[1]
    ua_abstract = re.search(r"\\begin\{abstract\}(.*?)\\end\{abstract\}", ua, re.S)[1]
    checks["abstract_numbers"] = plain_numbers(en_abstract) == plain_numbers(ua_abstract)
    checks["abstract_math_and_citations"] = (
        math_tokens(en_abstract) == math_tokens(ua_abstract) and
        REFERENCE.findall(en_abstract) == REFERENCE.findall(ua_abstract)
    )
    return checks


def negative_controls(en, ua):
    mutations = {
        "numerical_value": ("0.00675308135", "0.00675308136", "all_math_tokens"),
        "math_sign": (r"-3K_{\mathrm{II}}", r"+3K_{\mathrm{II}}", "all_math_tokens"),
        "missing_citation": (r"\cite{Bowie1956}", "", "citation_groups"),
        "missing_figure": (GRAPHIC.search(ua)[0], "", "graphics_identity_and_order"),
        # AAM-UA7: the localized units remain checked one-to-one
        "wrong_unit_in_display": (r"2h_0=0.054025\ \text{мм}", r"2h_0=0.054025\ \text{м}", "all_math_tokens"),
        "wrong_unit_in_text": (r"$3.04~\text{мкм}$", r"$3.04~\text{мм}$", "all_math_tokens"),
        "wrong_unit_in_table_header": ("[мкм]", "[мм]", "table_units"),
        "display_punctuation": (r"\Delta a}=1.00.", r"\Delta a}=1.01.", "all_math_tokens"),
    }
    results = {}
    for name, (old, new, check) in mutations.items():
        assert old in ua
        results[name] = {"detected": not compare(en, ua.replace(old, new, 1))[check],
                         "check": check}
    return results


def inventory(en):
    clean = comments_removed(en)
    citations = re.findall(r"\\cite\{([^}]+)\}", clean)
    keys = list(dict.fromkeys(key.strip() for group in citations for key in group.split(",")))
    bib = pinned(BIB).decode()
    dois = {}
    for key in keys:
        record = re.search(r"@[A-Za-z]+\{" + re.escape(key) + r",(.*?)(?=\n\})", bib, re.S)
        assert record, key
        doi = re.search(r"\bdoi\s*=\s*\{([^}]+)\}", record[1], re.I)
        dois[key] = "https://doi.org/" + doi[1] if doi else None

    def environments(regex):
        return [{"start_line": clean.count("\n", 0, m.start()) + 1, "latex": m.group()}
                for m in re.finditer(regex, clean, re.S)]

    numbered = [m for m in DISPLAY.finditer(clean) if m[1]]
    return {
        "base_commit": BASE,
        "english_manuscript": identity(EN),
        "author_metadata": METADATA,
        "companion_sources": [identity(f) for f in
                              ("paper/submissions/aam/README.md", "paper/submissions/aam/check_scientific_freeze.py",
                               "paper/main.tex", BIB, *FIGURE_PROVENANCE)],
        "audited_macros_latex": numerical_macros(en).split("\n", 3)[-1].strip(),
        "sections": [{"command": m[1], "title": m[2]} for m in re.finditer(
            r"\\(section\*?|subsection)\{([^}]+)\}", clean)],
        "display_equations": [{"numbered_environment": m[1], "latex": m[0],
                               "labels": re.findall(r"\\label\{([^}]+)\}", m[0])}
                              for m in DISPLAY.finditer(clean)],
        "equation_counts": {
            "numbered_environments": len(numbered),
            "numbered_equations": sum(1 + (len(re.findall(r"\\\\", m[2])) if m[1] == "align" else 0) for m in numbered),
            "unnumbered_display_equations": sum(m[1] is None for m in DISPLAY.finditer(clean)),
            "equation_labels": len(re.findall(r"\\label\{eq:[^}]+\}", clean)),
        },
        "tables": environments(r"\\begin\{table\}.*?\\end\{table\}"),
        "table_numerical_cells": table_cells(en),
        "figures": environments(r"\\begin\{figure\}.*?\\end\{figure\}"),
        "graphic_assets": [identity(f"paper/{f}") for f in GRAPHIC.findall(en)],
        "citation_groups": citations,
        "cited_bibliography": dois,
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--write-audit", action="store_true")
    parser.add_argument("--source-only", action="store_true",
                        help="Skip verifying localized figure outputs while developing figures")
    args = parser.parse_args()
    en = pinned(EN).decode()
    ua = (EDITION / "main_ua.tex").read_text(encoding="utf-8")
    assert blob(EN) == EN_BLOB
    report = compare(en, ua)
    report["bibliography_bytes"] = pinned(BIB) == (EDITION / "references.bib").read_bytes()
    report["protected_paths_unchanged"] = not git("diff", "--name-only", BASE, "--", ".",
                                                  ":(exclude)paper/ukrainian_review/**").strip()
    manifest = inventory(en)
    if not args.source_only:
        # AAM-UA7: every graphic is a reconstruction (or retained vector layer) recorded in FIGURE_SOURCE_MANIFEST.json;
        # the English asset must be the pinned one, the packaged Ukrainian file must match the manifest, and the
        # independent data verification (FIGURE_DATA_VERIFICATION.json, incl. negative controls) must have passed.
        figures = {a["path"].removeprefix("paper/ukrainian_review/"): a for a in json.loads(
            (EDITION / "audit/FIGURE_SOURCE_MANIFEST.json").read_text(encoding="utf-8"))["assets"]}
        verification = json.loads((EDITION / "audit/FIGURE_DATA_VERIFICATION.json").read_text(encoding="utf-8"))
        report["all_graphic_hashes"] = verification["pass"] and len(figures) == len(manifest["graphic_assets"]) and all(
            f["sha256"] == figures[path]["englishAsset"]["sha256"] and
            hash_bytes((EDITION / path).read_bytes()) == figures[path]["ukrainianSHA256"] and
            figures[path]["independentVerification"] == "PASS" and
            figures[path]["panelId"] in verification["baseline"]["panels"]
            for f in manifest["graphic_assets"] for path in [f["path"].removeprefix("paper/")]
        )
    controls = negative_controls(en, ua)
    result = {"checks": report, "negative_controls": controls,
              "allowed_translations_inside_math": MATH_TEXT,
              "number_of_inline_and_display_math_spans": len(math_tokens(scientific_body(comments_removed(en)))),
              "body_citation_groups": len(manifest["citation_groups"]),
              "number_of_bibliography_entries": len(manifest["cited_bibliography"]),
              "number_of_tables": len(manifest["tables"]),
              "number_of_figures": len(manifest["figures"]),
              "number_of_graphics": len(manifest["graphic_assets"]),
              **manifest["equation_counts"]}
    print(json.dumps(result, ensure_ascii=False, indent=2))
    if args.write_audit:
        for path, data in [(EDITION / "SOURCE_MANIFEST.json", manifest),
                           (EDITION / "audit/PRESERVATION_CHECKS.json", result)]:
            path.write_text(json.dumps(data, ensure_ascii=False, indent=2) + "\n", encoding="utf-8", newline="\n")
    if not all(report.values()) or not all(c["detected"] for c in controls.values()):
        raise SystemExit(1)


if __name__ == "__main__":
    main()
