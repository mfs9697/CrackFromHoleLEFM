# Issue #78: self-contained vector graphics and tight publication crops

Base: `paper-scientific-batch6`, commit `2e4cf7e39693a45470545f27e08678606ac84b80`.
The feature branch is independent; no upstream PR was merged. All scientific
data, accepted source vectors, solver/qualification algorithms, gates, tables,
conclusions, authorship and affiliation remain unchanged.

## Review outputs and approval boundary

The candidate PDFs are in `paper/reviews/issue78/after/figures/`, using the same
31 asset names as `paper/main.tex`. The original publication PDFs remain in
`paper/figures/`; byte-identical signed copies are in the review's `before/`
directory. **Nothing was promoted into the accepted publication folder.**
This follows the issue's explicit requirement: “The author must inspect the
new publication PDFs and approve before promotion/merge.”

Review files:

- `comparison-1.png` through `comparison-4.png`: all 31 before/after panels.
- `geometry_report.json`: old/new page boxes, axes/content/legend bounds,
  margins, vector status, hashes, translation matrices and decoration changes.
- `verification.json`: independent source/geometry checks.
- `compiled/manuscript_before.pdf` and `compiled/manuscript_after.pdf`:
  complete manuscripts using the corresponding vector assets.
- `compiled/figures_7_9_before.pdf` and `compiled/figures_7_9_after.pdf`:
  three-page, native manuscript-size LaTeX proofs, with the real captions.
- `manuscript-contact-1.png` through `manuscript-contact-3.png`: all final pages.

## One reproducible route

From the repository root in MATLAB:

```matlab
addpath(genpath(pwd));
R = review_issue78_figures();
% redraw_all_figures() uses this same review-only route.
% Optional separate destination: review_issue78_figures('OutputRoot',folder)
```

This generates only two MATLAB vector ylabel stamps (M1 minus M0), then applies
the shared Python vector transformer, checks sources/geometry, creates PNGs
and contact sheets, and prepares TeX proofs. It does not read/reconstruct
numerical histories, assemble a mesh, run FEM, reevaluate EDI, fit COD, or
recompute MTS/coordinates/SIFs. Python needs `pypdf`, `pdfplumber`, `reportlab`
and Pillow; Poppler supplies PNG rendering. Existing runtime overrides
`PUBLICATION_PYTHON` and `PUBLICATION_PDFTOPPM` are supported.

Compile the prepared TeX with an installed LaTeX engine, for example from
`paper/reviews/issue78`:

```text
tectonic manuscript_after.tex --outdir compiled
tectonic figures_7_9_after.tex --outdir compiled
```

The bundled Tectonic compiler was used for the delivered proofs. The isolated
`review_figures_7_9()` command also uses the signed vector route and keeps
outputs outside the accepted publication folder.

After inspecting and approving the candidates, the author can explicitly run
`promote_issue78_figures(R,'AuthorApproved',true)`. The promotion helpers reject
calls without that approval argument, preflight all reviewed hashes/files, and
copy graphics/metadata only. No promotion or merge was performed for this PR.

## Source map and preservation

| Panels | Original producing routine / numerical provenance | Action here |
|---|---|---|
| Figure 1 | `plot_specimen_geometry.m`; approved vector composition | Copy vector operators; crop full content. |
| Figure 2 (3) | `main_tip_core_mesh_level_comparison.m`; accepted mesh illustrations | Crop existing vectors only; no mesh generator called. |
| Figures 3--6 and S1 (15) | `plot_existing_manuscript_figures.m`; saved `evidence_exact.mat` records | Existing vector operators retained, legends/decorations normalized, tight crops. No extraction/refitting/recalculation. |
| Figure 7 (4) | `plot_tip2h0_vs_reference_publication.m`; corrected portable isolated-core evidence | Accepted vectors retained. Legends identify h0/2h0 under **M0 exterior**. |
| Figure 8 (4) | Original retained M1 vector evidence; original full M1 history unavailable | No history fabricated. Genuine legend keys/glyphs reflowed; displayed `reference` becomes **M0** using the retained font's existing M and 0 glyphs. Only the two display ylabels are typeset as M1 minus M0. |
| Figure 9 (4) | Accepted native 4/2/1-mm trajectory vector; companions from saved increment plotting records | All accepted vectors retained. Trajectory and density legends identify increments under **Scaled M1**. |

`paper/data/retained_figure_7_9/` was not changed. Neither were the numerical
snapshots, old historical CSVs, paths, provenance keys, or MATLAB reference
variable names. M0 is a display name for the requested reference exterior law,
not a new curve/data family. M0 does not identify only one Figure 7 curve.

## Legend and crop implementation

`publication_pdf_assets.py` tags individual PDF paint operations during
analysis. It extracts balanced decoration primitives while preserving the
original line/marker keys, paths, glyphs, resources and graphics states. Old
legend paints are disabled individually: there is no rectangular white mask
that could erase a curve behind a legend. M0 substitution edits only the
derived legend strings in the original embedded font. New label stamps contain
only their visible glyphs; their MATLAB page backgrounds are not overlaid on
ticks or data.

Legends use clear internal corners, with original swatches and uniformly
scaled glyphs. Segment/rectangle checks catch crossing curves even when both
segment endpoints are outside a proposed legend. Hole regions, markers,
physical boundaries, zero/crossing lines and critical annotations are obstacles.
The exporter refuses a layout if no unobstructed internal position exists.

All final PDF crops use a 2-mm safety envelope. Shared envelopes maintain
companion axes sizes/alignment. Rotated ylabel decorations are aligned to the
same left column; five unchanged scientific multipliers are moved into clear
internal corners, retaining their factors and tick values. COD/quality top
rows without xlabels do not inherit their lower row's empty xlabel band.
The resulting individual exterior margins are **2.00--2.99 mm**.

Only page translations are applied to scientific plot bodies: no anisotropic
stretch, changed axes limits, internal-axes crop, digitization, or raster
replacement. The 31 final PDFs have matching MediaBox/CropBox, and their PNGs
and `.layout.json` describe those final assets.

The live MATLAB export helpers now use `publication_place_internal_legend.m`
and the same crop/PNG route. `normalize_retained_sensitivity_pdf.py` no longer
enforces an outside-strip policy. The layout configuration, verification JSON
and regression checks certify internal placement instead of `legendOutsideData`.

## Executed verification

- 54 independent source/panel checks; 31 vector PDFs; zero legend collisions.
- 14,860 unchanged painted primitives matched to their originals under crop
  translations, including curves, frames, ticks and unaffected labels.
- All 23 numerical/retained-source files in the signed data manifest remain
  byte-identical. Changed source evidence is rejected before candidate output.
- MATLAB exercised internal legend placement, viewport preservation, tight
  PDF/PNG export and the retained M0/M1 exporter successfully.
- Full before/after manuscripts compile (25 / 24 pages); both three-page
  Figures 7--9 proofs compile with no overfull boxes/oversized floats.
- All 31 before/after panels and all 24 final manuscript pages were inspected.
  Figure 9 retains natural subfigure heights and its separating paragraph/gap;
  its main caption does not overlap subcaptions.

One **pre-existing** 1.327-pt introductory-text overfull hbox occurs in both
whole-manuscript builds. It is outside this figure task and was not altered.
There are no new figure-related overflow warnings. The missing raw M1 history
still limits independent numerical replay; this route reproduces its graphical
evidence without claiming to recover that history.

The older `verify_isolated_tip_correction.py` prose assertion requiring the
literal text `eta_e=1` / “connectivity need not be identical” also fails on
the unmodified batch-6 manuscript. This pre-existing prose/test mismatch is
reported rather than weakening its numerical checks. The independent source,
vector-geometry and unchanged-manuscript-scope checks above pass.

Manuscript edits are limited to M0 display terminology in Figure 8 and removal
of the obsolete fixed-height containers in Figures 7/8. Panel order, labels,
arrangements, all numerical values and scientific conclusions are preserved.
