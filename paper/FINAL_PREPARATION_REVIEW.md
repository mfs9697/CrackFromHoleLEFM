# Final scientific and graphical preparation, 2026-10-10

PR #69 was open when work began. Its source branch
`paper-use-approved-specimen-figure` was used, retaining the author's approved
specimen composition, the symbolic problem statement, all inline figure
environments and all 71 pre-existing labels. No branch was merged.

## Executed work

- MATLAB vector schematic: rectangular 2A-by-B plate, generic circular hole,
  lower-left origin O, explicit axes, unsigned outward tensile sigma arrows,
  free lateral edges and an explicitly illustrative curvilinear crack.
  The approved reference PDF is preserved separately. No numerical values or
  computed trajectory are inserted into the sketch.
- Unified `redraw_all_figures` plot-only driver and shared physical export/style
  helpers. Regenerated labels/legends/annotations are set at 11 pt and ticks at
  9.5 pt on PDF pages matching inline LaTeX widths. Interpreters are explicitly
  LaTeX, including ticks. The workflow distinguishes retained assets from
  regenerated data-driven panels; it does not invoke meshing or numerical
  extraction routines.
- Figures 7/8/9 share one wide trajectory and three companion panels. The
  retained portrait M1 source prevents a legible two-by-two arrangement;
  the chosen common topology preserves the quantities and spatial scales.
- Dark-blue unboxed hyperlinks, a 190-word number-free abstract, an Introduction
  contribution/organization paragraph, and a 47-family notation register.
  Local MTS turns, frozen-frame absolute directions, global directions and
  between-study differences are distinguished. Eta multipliers replace option
  names in the article; software options remain unchanged.
- Presentation precision is reduced only in LaTeX displays. Full-precision
  evidence and numerical acceptance checks are retained. Scientific prose
  emphasizes recursive direction selection, numerical sensitivity and model
  limits rather than the chronology of archive audits.

## Verification actually executed

1. MATLAB `redraw_all_figures` executed, then executed again into a scratch
   reproduction directory. Its reports record zero physical solves, zero mesh
   studies, zero EDI evaluations and zero COD fits. It retained eight source-
   limited panels explicitly. Subsequent plot-only margin/annotation repairs
   were executed and visually checked; they do not change numerical arrays.
2. MATLAB `test_isolated_tip_publication` executed: nine input/provenance guards
   passed, including missing, historical, incomplete and unaccepted inputs.
3. `verify_isolated_tip_correction.py` executed: 59 saved-data/table checks and
   all ten inline environments/31 graphic references passed. Display rounding
   is checked against the unrounded saved values, not by patching those values.
4. `verify_manuscript.py` and `verify_final_presentation.py` executed: the
   abstract, notation, unique labels, vectors, inline references and physical
   export dimensions passed. The compiled PDF has 63 unboxed functioning link
   annotations and 77 unique source labels, including six new section labels.
5. The optional strict local-archive check was executed and correctly stopped:
   the original full reference archive is unavailable. Portable evidence checks
   passed; raw physical field revalidation is **not** claimed.
6. Tectonic compiled the complete manuscript. Final PDF pages were rendered and
   visually inspected, with close checks of the specimen, COD legends, P24
   annotation and sensitivity layouts. There are no overfull text boxes or
   unresolved references in the final PDF. Narrow-caption underfull advisories
   are non-fatal and do not clip content.
7. Baseline hashes verify 244 original evidence/solver files unchanged. Original
   bibliography bytes and citation order are preserved. Numerical results,
   solver algorithms and physical acceptance gates are not modified.

## Unresolved, explicitly retained limitations

The raw M1 history, mesh-hierarchy coordinates and complete native increment
trajectories remain unavailable. Accepted plotting records survive for the
increment's three comparison panels; their native trajectory is not rebuilt
from common-length points. The retained wide vector is a previously verified
export of the same accepted increment family, with separate hash provenance.

Eight retained panels cannot be fully restyled to the exact 11/9.5-pt target or
certified for their original MATLAB interpreters without raw sources. Their
original embedded fonts remain a documented exception. No numerical data were
invented to disguise that limitation. Supplying the missing saved histories
would allow exact typography regeneration without a new physical solve.

The model remains quasi-static LEFM direction selection with prescribed
increments. Formal continuum error bounds, an exact local-symmetry point and
a uniquely isolated boundary mechanism are not asserted. Stronger versions of
those interpretations would require evidence outside this task and were not
implemented.

## Commands and reports

```matlab
addpath(genpath(pwd));
R = redraw_all_figures();
```

```text
python paper/verify_manuscript.py
python paper/verify_isolated_tip_correction.py
python paper/verify_final_presentation.py
```

See `FINAL_FIGURE_WORKFLOW.md`, `NOTATION_AND_EDITORIAL_REVIEW.md`,
`notation_register.csv`, `figures/redraw_manifest.json` and
`final_preparation_verification.json` for detailed provenance and scope.
