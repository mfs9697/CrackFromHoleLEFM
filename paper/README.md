# Pilot manuscript: crack trajectory from a circular hole

## Manuscript figure and problem-statement conventions (2026-10-10)

All ten `figure` environments, including subfigure captions and labels, are
maintained **directly in `paper/main.tex`**. External files under
`paper/figures/` hold graphical PDF assets only; no external
`figure.tex` wrappers are included or required. The approved specimen
schematic is displayed at `0.78\\textwidth`, selected against the
11-point manuscript text in a PDFLaTeX proof. Other graph-panel widths
retain their audited arrangements pending separate font-scale checks.

The **Problem formulation** introduces dimensions, hole coordinates,
material, loads, initiation threshold and crack-increment variables
symbolically. Physical numerical values are consolidated under
**Numerical parameters and initiation state** in the implementation
section; the first-tip SIFs and turn are reported in results. The
schematic measures the center ordinate from the bottom edge, while the
computational frame is shifted upward by `B/2`. The right-edge
clearance is consistently `2A-x_{23}`.


The author-approved specimen/loading schematic is included as the first figure in
`main.tex` from `figures/geometry_loading/specimen_geometry.pdf`. It is a
PDF-wrapped black-and-white rendering of the supplied original image, not a
newly generated mechanics result. The drawing measures $y_c$ from the lower
plate edge; the numerical frame shifts the vertical origin up by $B/2$.
Thus the same physical hole is at (170,80) mm in the sketch and
(170,-20) mm in the numerical frame. The curved crack is illustrative.
The old procedural figure generator has been retired.


Current correction (2026-10-09): the manuscript retains its existing structure
and now uses the completed isolated-core tip-resolution study. See
[ISOLATED_TIP_REPLACEMENT.md](ISOLATED_TIP_REPLACEMENT.md) for portable reproduction
and [FIGURE_CONSISTENCY_AUDIT.md](FIGURE_CONSISTENCY_AUDIT.md) for all Figures 1--8
and S1. Run `python paper/verify_isolated_tip_correction.py` for the current
portable data/table/figure checks. The older `verify_manuscript.py` is a
historical full-local-archive audit with obsolete pilot layout assumptions;
its reference manifests are preserved, not overwritten by this correction.
The historical pilot notes below are retained as development history.

Manuscript layer prepared on `paper-pilot-crack-trajectory`, based on the
current plotting/development snapshot
`cfdf0f110010d688a4e3c48f6d88a00fd17dc698` (2026-10-06).
The pilot manuscript layer is based on audited numerical archives. The
reference solver/extractor implementation is unchanged; subsequent branches
have added a qualified coarser-exterior M1 trajectory and publication plots
derived from the saved reference and M1 runs.

## Title alternatives, followed by the selected title

1. **Incremental LEFM Prediction of a Crack Trajectory from a Circular Hole**
2. **Mode-Mixity Reversal Along a Crack Trajectory Approaching a Free Boundary**
3. **Current-Tip Direction Selection for a Hole-Initiated Crack**

Title 1 is selected. It states the mechanics problem directly without
promoting the boundary interpretation into an established causal result.
The reversal is the central result in the abstract and discussion.

## Contents

- `main.tex`: complete single-column article with governing equations,
  implemented MTS rule, characteristic states, five main figures (including
  the independent M1 exterior-mesh sensitivity figure), and one supplementary
  quality figure.
- `references.bib`: explicit literature TODOs; no fabricated entries.
- `FIGURE_PLAN.md`: compact figure selection and source mapping.
- `EVIDENCE_AND_PROVENANCE.md`: numerical claim-to-source mapping and limits.
- `REMAINING_WORK.md`: actual gaps before journal submission.
- `data/`: audited CSV/JSON, exact compact MAT snapshot, original accepted
  Stage-I source, derived metrics, and source SHA256 manifest.
- `export_pilot_evidence.m`: read-only manuscript data extraction/audit.
- `plot_existing_manuscript_figures.m`: redraws Figures 1--4 and S1 in
  MATLAB as separate vector panels for LaTeX subcaptions; no FE solve or
  extraction replay is performed.
- `verification/crack_path/plot_m1_vs_reference_publication.m`: generates
  the separate vector panels for the M1 exterior-mesh sensitivity figure.
- `sync_manuscript_data.py`: updates only embedded tables, coordinates, and
  numerical macros in the TeX source from the audited data snapshot.

The article distinguishes P1's previously accepted seed from the 22 accepted
P2–P23 physical files. P24 appears only as dashed qualified-unsolved geometry.
The linear zero estimate is always an interpolation, never a solved state.
The conclusion concerns a fixed 4-mm sequence; no path-convergence claim is made.

## Build and review

The manuscript uses external vector PDF panels under `paper/figures/`, with
panel letters and panel captions supplied by LaTeX `subcaption`. The plotted
numerical evidence remains in `paper/data/`. `references.bib` is intentionally
inactive until entries have been verified. Before compiling after a clean
checkout, generate the MATLAB panels from the saved evidence archives.

From the repository root, generate the figure panels first:

```matlab
addpath(genpath(pwd));
plot_existing_manuscript_figures();
plot_m1_vs_reference_publication();
```

Then, from `paper/`, compile `main.tex` twice with pdfLaTeX or XeLaTeX.
Alternatively:

```text
tectonic -X compile main.tex
```

This pilot was compiled with the bundled Tectonic compiler. The built-in
preview compiler returned `Unable to find standard directories for platform`;
the source was preserved and the standalone PDF was compiled separately.
The local MiKTeX latexmk wrapper also lacks its Perl runtime; neither global
configuration nor installed software was changed to work around it.
Rendered PDF pages were inspected for layout and figure legibility.
Paper-local line-ending rules preserve the audited text hashes across
Windows checkouts without changing any scientific source outside this layer.

To re-audit the investigator archive in MATLAB, without any physical solve:

```matlab
addpath(genpath(pwd));
E = export_pilot_evidence( ...
    'RunDir',fullfile(pwd,'verification','crack_path','final_clean_run'), ...
    'FrozenStateFile',fullfile(pwd,'paper','data','accepted_stage1_source.mat'));
```

Then run `python paper/sync_manuscript_data.py`. If the evidence changes,
review the narrative as well as the refreshed numbers; this is not a tool
for automatically authoring a different physical study. The reader stops
on archive/state disagreement rather than silently reconciling it.

## Candid assessment

1. **Strongest contribution:** a reproducible current-tip-driven trajectory
   showing a mode-II sign reversal while KI continues to grow. All stored
   COD definitions corroborate the P21–P22 sign-change bracket.
2. **Weakest point:** the trajectory still has only one crack-increment
   size. The independent M1 run now demonstrates negligible sensitivity to
   substantial coarsening of the exterior mesh, but it does not test local
   tip-core/EDI resolution or increment-size convergence.
3. **Likeliest reviewer scrutiny:** the location/robustness of the local-
   symmetry crossing and the causal attribution to the free boundary.
   The 84.15-mm value is interpolation; boundary distance and curvature
   co-evolve without a separating control.
4. **Most valuable next calculation:** a matched physical 4/2/1-mm
   increment study, with independently controlled mesh/extraction resolution
   and comparison at common crack lengths. Do not couple every numerical
   change into one refinement and call it increment convergence.
5. **Submission readiness:** the evidence is sufficient for a serious
   internal fracture-mechanics pilot and an informative draft. It is not
   yet ready for journal submission: verified literature, increment-size
   sensitivity, and local tip/EDI-domain checks remain necessary. No formal
   journal “pilot submission” category is assumed.

The full prioritized list is in `REMAINING_WORK.md`.
