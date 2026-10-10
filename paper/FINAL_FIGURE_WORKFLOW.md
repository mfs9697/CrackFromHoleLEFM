# Final publication figure workflow

Work continues on PR #69's source branch so the approved specimen drawing,
symbolic formulation and inline figure environments are preserved. All figure
environments remain in `main.tex`; no `figure.tex` wrapper is reintroduced.

From the repository root:

```matlab
addpath(genpath(pwd));
R = redraw_all_figures();
```

This command is plot-only: no stiffness assembly, FEM solve, mesh generation,
SIF evaluation or COD refit. `redraw_manifest.json` explicitly identifies
regenerated and retained graphics and records their sources and hashes.

## Source map

| Figure | Source | Action |
|---|---|---|
| 1, geometry/loading | `plot_specimen_geometry.m`; author's approved composition retained as `specimen_geometry_approved_reference.pdf` | MATLAB vector drawing, symbolic dimensions only, unsigned outward sigma arrows, origin O and illustrative crack. Drawing origin is lower-left; numerical y is shifted by -B/2. |
| 2, mesh hierarchy | `verification/crack_path/main_tip_core_mesh_level_comparison.m` and audited `build_stage2_scaled_audited_core.m` | Regenerated separately in mesh-only mode; exact element/COD fingerprints are checked. Excluded deliberately from strictly plot-only `redraw_all_figures()`. |
| 3, reference path | `evidence_exact.mat`, `plot_existing_manuscript_figures.m` | Regenerated. Solid P0-P23 and dashed qualified-unsolved P24; interpolation diamond not a solved state. |
| 4, KI/KII/absolute direction | Same saved reference table | Regenerated; theta is absolute in the frozen frame, not a local MTS turn. |
| 5, late mode mixity/MTS turn | Same saved reference table | Regenerated; Delta theta belongs to the current accepted field and predicts the next segment. |
| 6, COD comparisons | Same 176 saved fits | Regenerated from stored differences; no fitting/extraction repeated. |
| 7, isolated-core sensitivity | Validated portable isolated dataset and unchanged reference | Regenerated with existing plotting routine. Core multiplier eta_c varies, requested exterior multiplier eta_e=1; connectivity may change. |
| 8, M1 exterior sensitivity | Verified existing four vector panels | Retained: raw M1 history remains unavailable. No replacement history is invented. |
| 9, increment sensitivity | Saved accepted plotting records recovered from the executed 2026-10-09 figure audit | Three companion panels regenerated from exact saved metrics/native records. Complete native vertex histories remain unavailable: the previously verified wide trajectory export from commit 81f1608 is retained, with separate provenance. Common-length coordinates are not interpolated into missing native trajectories. |
| S1, numerical acceptance | Saved reference solver and mesh summaries | Regenerated; gates unchanged, no P24 failed physical field plotted. |

The compact `increment_plot_records.mat` is a copy of the existing plotting
result, not a new numerical result. Its 23 exact-common-length comparison rows
and 161 native q/turn-density rows were checked against the published figures
and tables during the previous accepted-archive audit. The original archive
paths remain in its provenance; those physical archives are now absent.

## Figure 2: approved mesh panels (separate mesh-only workflow)

Figure 2 has three reviewed local PDF panels. The publication source uses
`figures/mesh_levels/tip_mesh_H2.pdf`, `tip_mesh_H1.pdf`, and
`tip_mesh_H0.pdf` with concise subcaptions `4h_0`, `2h_0`, and `h_0`.
These are **mesh illustrations**; the paper distinguishes the physically
qualified H1/H0 levels from the illustration-only H2 level (COD counts
10/15/12/9 fail the 12-native-point gate).

To promote **the already reviewed** PDF and PNG panels from their scratch
folder into the manuscript folder, run from the repository root:

```matlab
addpath(genpath(pwd));
P = promote_figure2_review();
```

The promotion helper checks the original T3 and T6 node counts, all four COD
point counts, the three tip sizes and the expected H2/H1/H0 qualification
statuses **before copying any files**. It does not run a solve or even a mesh
builder. The output is six locally updated graphic files
(`tip_mesh_H2/H1/H0.pdf` and `.png` in `paper/figures/mesh_levels`).
Commit these with GitHub Desktop to make the approved binaries available on
GitHub; merely pulling the PR cannot upload MATLAB-generated binary files.

If generation must be repeated, it remains a separate **mesh-only**
diagnostic, not part of `redraw_all_figures`:

```matlab
addpath(genpath(pwd));
R = main_tip_core_mesh_level_comparison( ...
    'OutputDir', fullfile(pwd,'paper','figures','mesh_levels_review'));
disp(R.summary);
```

Examine the new PDFs before running the promotion helper again.

## Physical sizing and typography

`publication_style.m` defines the manuscript's 162-mm text width and matching
0.94/0.485/0.315 export widths. `publication_export_axis.m` prints a fixed-size
vector PDF page rather than tightly cropping axes and subsequently shrinking
unknown margins. At the corresponding inline LaTeX width, regenerated labels,
legends and annotations are 11 pt and ticks are 9.5 pt. Mathematical text and
ticks explicitly use MATLAB's LaTeX interpreter. MATLAB's Computer Modern
math/text rendering is compatible with the manuscript's Latin Modern layout.

Figures 7, 8 and 9 share a wide trajectory above three companion panels.
A two-by-two layout was considered but rejected because the retained portrait
M1 panels would either exceed the page or require substantially smaller text.
The common arrangement preserves quantities and spatial aspect ratios while
allowing the available-data panels to be rendered at their final physical size.

**Retained-source limitation:** original M1 and increment-trajectory graphics
cannot be fully regenerated without their raw source coordinates. Their
embedded fonts are retained, not represented as having passed the new
interpreter/font checks. Figure 2 is no longer in this category: its audited
deterministic mesh-only generator is available and has passed the exact T3,
T6-node and four-window COD-sampling fingerprints.
The wide increment asset is a verified graphical variant of the same accepted
4/2/1-mm histories, not a substituted numerical family. This limitation is
recorded for review rather than concealed by numerical reconstruction.

Scientific values, accepted evidence snapshots, numerical solver sources and
acceptance gates are unchanged. Panel letters are generated by LaTeX, and all
existing figure/table/equation labels are retained.
