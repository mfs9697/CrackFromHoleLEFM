"""Check that journal formatting has not altered the approved scientific paper.

Run from any folder:
    python paper/submissions/aam/check_scientific_freeze.py

This read-only check compares the full scientific body, audited macro block,
abstract, figures and citation keys against paper/main.tex. Exact explicit
editorial differences (bibliography, Appendix Fig. S1 -> Fig. 10) are
reversible and hence do not waive checks on numeric or mechanical statements.
"""
from pathlib import Path
import re

HERE = Path(__file__).resolve().parent
PAPER = HERE.parents[1]
SOURCE = PAPER / "main.tex"
CANDIDATE = HERE / "manuscript.tex"


def between(s: str, left: str, right: str) -> str:
    if s.count(left) != 1 or s.count(right) != 1:
        raise ValueError(f"Expected exactly one {left!r} / {right!r}")
    return s.split(left, 1)[1].split(right, 1)[0]


def scientific_body(s: str, is_candidate: bool) -> str:
    start = r"\section{Introduction}"
    if s.count(start) != 1:
        raise AssertionError(f"Missing/unexpected body marker: {start}")
    body = start + s.split(start, 1)[1]
    if is_candidate:
        end = "\n% EDITORIAL HOLD:"
        if body.count(end) != 1:
            raise AssertionError("Candidate missing editorial-hold terminator")
        body = body.split(end, 1)[0]
        # These are the only permitted journal-specific changes in the body.
        swaps = (
            (r"\bibliographystyle{sn-basic}", r"\bibliographystyle{unsrt}"),
            (r"\bibliography{../../references}", r"\bibliography{references}"),
            (r"\section{Additional numerical diagnostics}",
             r"\section{Supplementary numerical diagnostics}"),
        )
        for new, old in swaps:
            if body.count(new) != 1:
                raise AssertionError(f"Unexpected occurrence of {new}")
            body = body.replace(new, old)
        anchor = r"\section{Supplementary numerical diagnostics}"
        appendix_figure_reset = (
            "\\renewcommand{\\thefigure}{S\\arabic{figure}}\n"
            "\\renewcommand{\\theHfigure}{S\\arabic{figure}}\n"
            "\\setcounter{figure}{0}\n"
        )
        body = body.replace(anchor, anchor + "\n" + appendix_figure_reset, 1)
    else:
        body = body.split(r"\end{document}", 1)[0]
    return body.strip()


def check() -> dict:
    src, dst = SOURCE.read_text(encoding="utf-8"), CANDIDATE.read_text(encoding="utf-8")
    assert scientific_body(src, False) == scientific_body(dst, True), (
        "Scientific body has changed beyond explicit journal formatting differences"
    )

    start, end = r"\newcommand{\qmix}", "% END AUDITED DATA"
    block = start + between(src, start, end) + end
    assert block in dst, "Audited numerical macro block differs"
    abstract = between(src, r"\begin{abstract}", r"\end{abstract}").strip()
    assert abstract in dst, "Abstract has been modified"

    graphic_pattern = r"\\includegraphics(?:\[[^\]]*\])?\{([^}]+)\}"
    sf = re.findall(graphic_pattern, src)
    df = re.findall(graphic_pattern, dst)
    assert sf == df and len(sf) == 31 and len(set(sf)) == 31, (
        "The 31 approved PDF figure inclusions differ"
    )
    for f in df:
        if not (PAPER / f).is_file():
            raise AssertionError("Approved figure missing: " + f)

    keyword_match = re.search(r"\\keywords\{([^}]*)\}", dst, re.DOTALL)
    assert keyword_match is not None, "Missing journal keyword field"
    keywords = keyword_match.group(1).replace("\n", " ")
    assert len([x.strip() for x in keywords.split(";") if x.strip()]) == 6, (
        "Archive of Applied Mechanics requires 4–6 keywords"
    )
    assert r"\author*[1]{\fnm{Danylo} \sur{Selivanov}}" in dst
    assert "fract@inmech.kyiv.ua" in dst
    assert r"\documentclass[pdflatex,sn-mathphys-num]{sn-jnl}" in dst
    assert r"\bibliographystyle{sn-basic}" in dst
    assert src.count(r"\cite{") == dst.count(r"\cite{"), "Citation calls differ"
    return {"graphics": len(df), "keywords": 6, "scientific_body": "exact",
            "numerical_macros": "exact", "abstract": "exact"}


if __name__ == "__main__":
    print("PASS:", check())
