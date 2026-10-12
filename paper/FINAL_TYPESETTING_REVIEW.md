# Final typesetting and notation review

Based on PR #80 head `799f2e190fa90dd41b7fee178763d66ffcaada29` on
`paper-issue78-promote-reviewed-figures`. The implementation branch is
`paper-final-typesetting-notation`; its draft PR targets that PR #80 branch.
No PR was merged. No physical calculations, plotting, figure promotion,
mesh reconstruction, or changes to numerical evidence were performed.

## Additional academic-language review (2026-10-12)

Following the initial PR #81 typesetting commit, a separate, author-requested
language and terminology pass was made **on the same draft branch**:

- Audited all 24 occurrences of *reference* in the typesetting manuscript:
  18 ambiguous numerical-baseline/comparison occurrences were replaced by
  **M0**, **production**, **baseline**, or **$h_0$**, according to meaning.
  The remaining six occurrences describe the **unit reference traction/loading**,
  not the M0 exterior-mesh law, and are intentionally retained.
- Figure 7 now names the **$h_0$ path** when discussing the interpolated
  crossing; Figure 8 names the **M0 path**. M0/M1 is reserved for
  alternative **exterior mesh-sizing laws**, not tip-resolution labels.
- Checked all seven uses of *finite-element*. They are appropriate
  hyphenated adjectives (mesh, displacement field, method, etc.);
  the standalone abbreviation **FE** is not used or needed. No
  abbreviations were introduced merely to shorten the manuscript.
- Tightened the Abstract and selected introductory and results
  passages. In particular, the Abstract now avoids implying that
  coupled crack-increment/mesh refinement isolates increment effects.
  No numerical value, formula, table, citation or accepted figure was changed.
- The reproducibility guard `verify_final_typesetting.py` now carries
  **25 exact, uniquely matched, authorized prose substitutions**; after
  reversing those and the previously approved material-domain correction,
  all remaining nonfigure prose must exactly match the original PR #80
  baseline. Its numerical token checker no longer treats a numeric
  **symbolic subscript** (such as `h_0`) as a measurement. With this
  refinement, both baseline and revised manuscripts contain the same
  **708 numerical tokens in the same order**; this is not a change to
  any scientific datum.

**Verification status:** The source-level equivalence of all 22 prose
substitutions and the 708-token numerical sequence was checked before
updating the preservation script. The original compiled `paper/main.pdf`
and the review artifacts below correspond to the **initial PR #81
typesetting version** and have **not** been recompiled for these later
editorial commits. A fresh local LaTeX compilation, visual check and
rerun of the preservation script are **required before merging**.
Do not interpret the historical zero-warning PDF log below as a compile
of the latest manuscript source.

## Follow-up: accumulated crack length and software-neutral PCG description

The author identified a remaining `MATLAB` mention after locally compiling
the Codex worktree. On the active PR #81 branch, Section 4.6 now describes
finite-precision PCG stagnation and the permitted single Krylov restart
without referring to a specific software package. Section 3.3 no longer
redefines `\Delta a`, which was introduced in Section 2.2, and now
states that `a_k=k\Delta a` is the accumulated crack-path length from
`P_0`, not the horizontal tip coordinate.

All three paragraph-level changes are included in the strict 25-entry
authorized-prose whitelist. Comparison against the PR #80 scientific
baseline establishes identical display equations, numerical tables,
citations, graphic paths, and the entire numeric-token sequence. This
follow-up has **not yet been compiled into a new local PDF**; the
committed `paper/main.pdf` is a historical typesetting proof, not
a proof of the current source. A fresh compilation is required before
merge. The choice of computational implementation can remain in the
public reproducibility documentation rather than in this methods paragraph.

## Corrections

| Defect | Result |
| --- | --- |
| Figures 7–9 narrow captions | Local ragged-right, 9-pt subcaptions, shorter wording. Mathematical definitions moved to full-width captions where needed. Panel arrangements, widths, PDF paths, and all graphics remain unchanged. |
| Figure 8(d) says “reference” in supplied PDF | The primary checkout contains a locally modified, stale panel. The new worktree uses the correct PR #80 asset, already identical to the approved Issue #78 candidate. All four active panels distinguish M0/M1. Neither local modified file was overwritten. |
| Section 4.4 redefines the material domain | It now refers to the existing material domain and confines the nonzero gradient of the radial weight to the interaction annulus. Equation (20), its normalization, and the description of nodal T6 interpolation/quadrature are unchanged. |
| Figure 2 interrupts the integral explanation | A section-boundary float barrier puts Figure 2 before Section 4.4; definitions and Equation (20) are adjacent on page 8. |
| Table 2 / Figure 5 whitespace | Both permit placement near their references on page 13, with a narrow barrier before the COD discussion. Figure 6 follows its discussion on page 14 without separating a sentence across an intervening figure. |
| Exterior-mesh sentence interrupted by figures | Figure 7 is submitted near its first trajectory comparison; a local space reservation keeps the complete exterior-mesh discussion on page 16, before Figure 8 on page 17. |
| Section 5.7 final sentence on next page | Its final paragraph stays together on page 19; the Discussion starts on page 20. A short reservation protects the Discussion opening, and a barrier after Table 6 avoids another float interrupting its following paragraph. |
| References before diagnostics | Appendix A / Figure S1 now precedes References on page 22. References conclude on page 23. Figure S1 and its unique hyperlink anchor are preserved. |
| Introduction overfull line | A paragraph-local 1-em emergency stretch removes the 1.32713-pt overflow without changing a word or suppressing warning thresholds. |

The supplied PDF and the Issue #78 review each have 24 pages; the corrected
article has 23. Figure numbers, table numbers, equation numbers, and labels
are unchanged. Figure 8 remains on a dedicated figure page; its balanced
whitespace is intentional. No arbitrary image trims, negative spacing,
blanket page clears, or new forced `[H]` placements were introduced.

Figure 7 compares two structured-tip resolutions under M0; Figure 8 compares
M0/M1 exterior laws; Figure 9 remains the coupled, dimensionless scaled-M1
4/2/1-mm family. The caption changes retain absolute-direction differences,
local MTS turning density, exact-common-length comparisons, native-grid
histories, and the distinction between interpolation and solved states.

## Graphics diagnosis and preservation

The stale local Figure 8(d) SHA-256 is
`08ffd494bff504438306fc43649570a35e5b9ad2c0fa5cc8835e499e55f27536`.
The approved/active PR #80 asset SHA-256 is
`dd1f8b5dc91ef05b97567f14f4e41b04a4eefd0a877a27fd188550d43b695098`.
This is a rebuild from already promoted graphics, not an additional asset edit.
The primary checkout's modified `paper/main.pdf` and modified panel are retained.

`reviews/final_typesetting/verification.json` records:

- 123 protected files in `paper/data` and `paper/figures`, byte-identical to PR #80.
- All 31 active graphic PDFs byte-identical to their approved Issue #78 candidates.
- 29 displayed mathematical blocks, six complete tables (including their captions),
  and the audited numerical macros unchanged.
- 773 numerical literals outside figure captions preserved in order.
- All nonfigure prose unchanged except the specifically authorized material-domain
  correction; no SIF, MTS, COD, interpolation, gate, interpretation, or conclusion edited.
- Citation keys/order, 81 LaTeX labels, graphics paths/order, and panel geometry unchanged.
- Sole author, affiliation, email, title/author PDF metadata unchanged.
- All 70 internal PDF hyperlinks resolve; Figure S1 uses its preserved unique anchor.
- Five in-memory regression mutations (equation, numerical macro, table value,
  citation, historical graphic path) are rejected. No test alters a file.

## Build, warnings, and visual review

The complete project was compiled with the existing LaTeX skill workflow and
Tectonic 0.17.0, including BibTeX/reference reruns and all promoted graphics.
`reviews/final_typesetting/review_manifest.json` contains input/output hashes,
comparison-page mapping, and final-pass warning records.

| Warning | PR #80 baseline | Final |
| --- | ---: | ---: |
| Narrow-caption underfull boxes | 11 | 0 |
| Introduction overfull box | 1 (1.32713 pt) | 0 |
| Converged LaTeX warnings / unresolved citations or references | 0 | 0 |

First-pass unresolved citation messages during Tectonic's automatic reruns are
transient; none remain after convergence. The report separates these from final
warnings. Poppler emitted environment font fallback diagnostics for Symbol and
ArialUnicode while rendering; no missing or clipped glyphs were observed, and
the verified plot/font assets were not changed.

All 23 final pages were rendered with Poppler at 110 dpi and visually inspected,
including individual detail checks of the rearranged results pages. No axes,
legends, subcaptions, main captions, mathematical displays, tables, or text overlap
or clipping was found. Character bounding-box inspection also found zero
characters outside the page. All significant figures and accepted/unsolved
state qualifications were retained. Font size, subfigure labels, notation, colors,
and table alignment were checked against the approved review.
`reviews/final_typesetting/visual_review.json` binds the final PDF and rendered
page hashes to that review. A fresh render after the last compile was
pixel-identical to all 23 reviewed page images.

The 14-page vector proof `reviews/final_typesetting/affected_pages_comparison.pdf`
places affected pages from the supplied PDF, the committed Issue #78 approved
review, and the final article side by side. It records the stale Figure 8(d)
legend, caption spacing, mathematical-flow corrections, pagination, and end matter.
The page mapping distinguishes the relocated appendix from the References pages.

## Reproduction

From the repository root, with Tectonic installed:

```text
python -c "from pathlib import Path; Path('build/final-typesetting/reproduced').mkdir(parents=True, exist_ok=True)"
cd paper
tectonic -X compile --print --outdir ../build/final-typesetting/reproduced main.tex
cd ..
python paper/verify_final_typesetting.py --pdf build/final-typesetting/reproduced/main.pdf --self-test
```

Verification requires `pypdf`; creating the comparison proof additionally requires
`reportlab` and `Pillow`. These are available in the bundled runtime. The article
and its checks need no investigator-local MAT archives or additional calculations.

The proof workflow accepts the author's PDF explicitly, never as a figure source:

```text
python paper/review_final_typesetting.py --supplied PATH_TO_SUPPLIED_PDF --before-log PATH_TO_BASELINE_COMPILE_JSON --after-log PATH_TO_FINAL_COMPILE_JSON --render-dir PATH_TO_FINAL_PAGE_PNGS
```

The compile JSON files follow the existing LaTeX skill's `--json` format. Baseline
compilation uses PR #80's unchanged `main.tex`; the middle proof column uses the
committed `reviews/issue78/compiled/manuscript_after.pdf`. Scratch builds/renders
are ignored; the final PDF, vector proof, and verification records are committed.

## Remaining limitations / author judgment

The existing `verify_manuscript.py` chain still fails at
`verify_isolated_tip_correction.py:71`, which requires the exact obsolete phrases
`$\eta_e=1$` and “connectivity need not be identical.” An in-memory comparison
confirms the identical failure at the unmodified PR #80 baseline and this branch;
see `reviews/final_typesetting/legacy_verifier_status.json`. Numerical checks
before that prose assertion pass. The legacy test was left unchanged; the new
typesetting safeguards independently protect scientific values and evidence.

No unresolved typesetting defect was identified. Whether diagnostics should
eventually be a separate journal supplement remains an author/journal decision.
For this revision they remain Appendix A with Figure S1, before the bibliography;
no separate publication structure has been created.
