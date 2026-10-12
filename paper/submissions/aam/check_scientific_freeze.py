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
        # Restore only the precise, author-confirmed funding addition.
        funding = (
            "% BEGIN AAM AUTHOR-CONFIRMED FUNDING\n"
            "\\section*{Declarations}\n"
            "\\noindent\\textbf{Funding.}\n"
            "This research was supported by the regular institutional funding\n"
            "of the S.P. Timoshenko Institute of Mechanics,\n"
            "National Academy of Sciences of Ukraine.\n"
            "No specific research grant was awarded for this work.\n"
            "% END AAM AUTHOR-CONFIRMED FUNDING\n"
        )
        assert body.count(funding) == 1, "Funding statement edited without authorization"
        body = body.replace(funding, "", 1)
        # Springer 3.1's table wrapper (center + threeparttable) cannot
        # close reliably with the float package's forced [H] placement.
        # Restore the exact two original [H] choices solely for comparison.
        placements = (
            (
                "\\begin{table}[htbp]\n\\centering\n\\small\n"
                "\\caption{One-step directional map for the centered-hole symmetric control.",
                "\\begin{table}[H]\n\\centering\n\\small\n"
                "\\caption{One-step directional map for the centered-hole symmetric control.",
            ),
            (
                "\\begin{table}[htbp]\n\\centering\n\\small\n"
                "\\setlength{\\tabcolsep}{4.0pt}\n"
                "\\caption{Successive differences between independently propagated paths in",
                "\\begin{table}[H]\n\\centering\n\\small\n"
                "\\setlength{\\tabcolsep}{4.0pt}\n"
                "\\caption{Successive differences between independently propagated paths in",
            ),
        )
        for journal_placement, original_placement in placements:
            assert body.count(journal_placement) == 1, "Approved table placement changed"
            body = body.replace(journal_placement, original_placement, 1)
        assert r"\begin{table}[H]" not in s, (
            "Springer Nature class does not support [H] on its wrapped tables"
        )
        # Reverse the single Springer table-caption width correction.
        # Springer wraps table captions in threeparttable; the natural
        # narrow tabular width otherwise produces a cramped caption.
        # These two substitutions change no caption, header, or data.
        table_one_layout = (
            (
                r"\begin{tabular*}{0.64\linewidth}{@{\extracolsep{\fill}}rr}"
                + "\n" + r"\toprule" + "\n"
                + r"$\theta_2$ [deg] & $F(\theta_2)=\theta_3$ [deg]",
                r"\begin{tabular}{rr}"
                + "\n" + r"\toprule" + "\n"
                + r"$\theta_2$ [deg] & $F(\theta_2)=\theta_3$ [deg]",
            ),
            (
                r"\bottomrule" + "\n" + r"\end{tabular*}"
                + "\n" + r"\end{table}"
                + "\n\nThe opening-mode intensity",
                r"\bottomrule" + "\n" + r"\end{tabular}"
                + "\n" + r"\end{table}"
                + "\n\nThe opening-mode intensity",
            ),
        )
        for journal_layout, original_layout in table_one_layout:
            assert body.count(journal_layout) == 1, (
                "Unexpected Table 1 layout change in AAM candidate"
            )
            body = body.replace(journal_layout, original_layout, 1)
        # Reverse three exact table-header line breaks used to bring
        # Tables 5 and 6 inside the journal text column. No header
        # content, units, or tabulated data may change.
        table_header_layout = (
            (
                r"\shortstack[c]{$q_{\mathrm{K},\max}/\Delta a$\\{}[mm$^{-1}$]}",
                r"$q_{\mathrm{K},\max}/\Delta a$ [mm$^{-1}$]",
            ),
            (
                r"\shortstack[c]{$(\Delta\theta/\Delta a)_{q_{\max}}$\\{}[deg/mm]}",
                r"$(\Delta\theta/\Delta a)_{q_{\max}}$ [deg/mm]",
            ),
            (
                r"\shortstack[c]{$\max|\delta\theta|$\\{}[mdeg]}",
                r"$\max|\delta\theta|$ [mdeg]",
            ),
        )
        for journal_header, original_header in table_header_layout:
            assert body.count(journal_header) == 1, (
                "Unexpected Table 5/6 heading change in AAM candidate"
            )
            body = body.replace(journal_header, original_header, 1)
        # Revert a 0.2-pt per-side padding reduction in Table 3
        # needed for a 1.292-pt Springer text-column overrun.
        tip_table_spacing = (
            r"\setlength{\tabcolsep}{2.6pt}" + "\n"
            + r"\caption{Structured-tip sensitivity at fixed",
            r"\setlength{\tabcolsep}{2.8pt}" + "\n"
            + r"\caption{Structured-tip sensitivity at fixed",
        )
        assert body.count(tip_table_spacing[0]) == 1, (
            "Unexpected Table 3 padding change in AAM candidate"
        )
        body = body.replace(*tip_table_spacing, 1)
        # These are the only permitted journal-specific changes in the body.
        swaps = (
            ("% The sn-mathphys-num class option selects sn-mathphys-num.bst automatically.\n"
             + r"\bibliography{../../references}",
             r"\bibliographystyle{unsrt}" + "\n" + r"\bibliography{references}"),
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
        body = body.replace(anchor + "\n", anchor + "\n" + appendix_figure_reset, 1)
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
    assert len([x.strip() for x in keywords.split(",") if x.strip()]) == 6, (
        "Archive of Applied Mechanics requires 4–6 keywords"
    )
    assert r"\author*[1]{\fnm{Danylo} \sur{Selivanov}}" in dst
    assert "fract@inmech.kyiv.ua" in dst
    assert r"\usepackage{amsmath,amssymb,bm,booktabs}" in dst, "Wrong/missing amssymb package"
    assert "pdfauthor={Danylo Selivanov}" in dst, "Missing PDF author metadata"
    assert "pdftitle={Incremental LEFM Prediction of a Crack Trajectory from a Circular Hole}" in dst, "Missing PDF title metadata"
    assert r"\documentclass[pdflatex,sn-mathphys-num]{sn-jnl}" in dst
    assert r"\bibliographystyle{sn-basic}" not in dst
    assert r"\bibliographystyle{sn-mathphys-num}" not in dst
    assert r"\documentclass[pdflatex,sn-mathphys-num]{sn-jnl}" in dst
    assert src.count(r"\cite{") == dst.count(r"\cite{"), "Citation calls differ"
    return {"graphics": len(df), "keywords": 6, "scientific_body": "exact",
            "numerical_macros": "exact", "abstract": "exact"}


if __name__ == "__main__":
    print("PASS:", check())
