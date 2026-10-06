# Pilot manuscript: crack trajectory from a circular hole

Manuscript layer prepared on `paper-pilot-crack-trajectory`, based on the
current plotting/development snapshot
`cfdf0f110010d688a4e3c48f6d88a00fd17dc698` (2026-10-06).
The numerical solver, mesher, extractors, MTS routine, and plotting routines
are unchanged. No physical solve or mesh generation was performed.

## Title alternatives, followed by the selected title

1. **Incremental LEFM Prediction of a Crack Trajectory from a Circular Hole**
2. **Mode-Mixity Reversal Along a Crack Trajectory Approaching a Free Boundary**
3. **Current-Tip Direction Selection for a Hole-Initiated Crack**

Title 1 is selected. It states the mechanics problem directly without
promoting the boundary interpretation into an established causal result.
The reversal is the central result in the abstract and discussion.

## Contents

- `main.tex`: complete single-column article, 193-word provisional abstract,
  seven keywords, governing equations, implemented MTS rule, characteristic
  states, four main vector figures, and one supplementary quality figure.
- `references.bib`: explicit literature TODOs; no fabricated entries.
- `FIGURE_PLAN.md`: compact figure selection and source mapping.
- `EVIDENCE_AND_PROVENANCE.md`: numerical claim-to-source mapping and limits.
- `REMAINING_WORK.md`: actual gaps before journal submission.
- `data/`: audited CSV/JSON, exact compact MAT snapshot, original accepted
  Stage-I source, derived metrics, and source SHA256 manifest.
- `export_pilot_evidence.m`: read-only manuscript data extraction/audit.
- `sync_manuscript_data.py`: updates only embedded tables, coordinates, and
  numerical macros in the TeX source from the audited data snapshot.

The article distinguishes P1's previously accepted seed from the 22 accepted
P2–P23 physical files. P24 appears only as dashed qualified-unsolved geometry.
The linear zero estimate is always an interpolation, never a solved state.
The conclusion concerns a fixed 4-mm sequence; no path-convergence claim is made.

## Build and review

The TeX is standalone: plots and data are embedded, so there are no external
image or table inputs. `references.bib` is intentionally inactive until
entries have been verified. The source can be opened in the built-in editor.

For a regular TeX installation, compile `main.tex` twice with pdfLaTeX or
XeLaTeX. Alternatively:

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
2. **Weakest point:** the trajectory has only one crack-increment size.
   Tight solver and regression tolerances do not establish path accuracy.
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
   yet ready for journal submission: verified literature, trajectory
   convergence, and late-tip extraction checks remain necessary. No formal
   journal “pilot submission” category is assumed.

The full prioritized list is in `REMAINING_WORK.md`.
