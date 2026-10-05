# SIF audit harness


## Current audit status — closed

The asymmetric tiny-Mode-II forensic audit is scientifically closed on the
audit branch. The accepted result for the audited problem is

\[
K_{II}/K_I \approx 1.066\times10^{-4}.
\]

The final evidence includes exact reflection-paired symmetry controls,
prescribed Williams-field extraction controls, Level-0 COD/EDI
cross-extraction, a Level-0-qualified SGS-PCG solver, and a three-scale
structured C03 convergence family at \(s=1,\ 1/\sqrt2,\ 1/2\).

See **[AUDIT_CLOSURE.md](AUDIT_CLOSURE.md)** for the complete evidence chain,
numerical recommendation, limitations, production guidance, and merge record.
See **[STEP_INDEX.md](STEP_INDEX.md)** for the canonical numbering map,
including the historical 39--43, 48/50/51, and restored Step61 records.

The exact pre-cleanup repository state is frozen on
`archive/sif-audit-closed-2026-10-05`.

The material below documents the historical verification harness and earlier
gates. Statements there describing EDI as a prototype or the audit as
unfinished should be read in their original chronological context.

This folder is the verification layer for comparing the historical
mirror-based circular J/mode-separation extractor with the interaction
equivalent-domain integral (EDI) extractor.

The central rule is simple: **solve the FEM problem once, then pass the
identical mesh and displacement vector to both SIF extractors.**

## First control case

`cfg_crack_path_two_leg_control.m` reproduces the current non-perforated
Crack-Path geometry/material/mesh scale as a two-leg *elastic* crack
control. The second leg is traction free here; this is intentionally not a
CZM solve.

Run the first gate directly:

```matlab
addpath(genpath(pwd));
R = main_step1_same_field_compare();
```

or configure the control explicitly:

```matlab
C = cfg_crack_path_two_leg_control(2.0);
R = run_crack_path_old_vs_edi(C);
```

The output contains:

- `R.solution`: the single FEM field used by both methods;
- `R.old`: `SIF_LEFM_circle2_debug` result and diagnostics;
- `R.edi`: `SIF_LEFM_interaction_EDI` result and diagnostics;
- `R.difference`: signed method-to-method differences.

## Historical status at the first control stage

At the beginning of the audit, the interaction EDI implementation was still a
prototype whose normalization, Mode-II sign, and auxiliary-field derivatives
required independent verification. The material below preserves that
chronological state. Those verification gates were subsequently completed;
for the final scientific status use `AUDIT_CLOSURE.md` rather than this early
stage description.

The older published two-segment benchmark remained a separate reproduction
gate at this point in the chronology.

## Literal-lattice crack geometry gate (Step 4C)

The approved parent is `build_literal_ring_lattice.m`: 64 equal angular
segments on **every** ring, alternating angular phase 0 and half a sector,
r0=0.005, r1=0.20, 46 radial intervals and q approximately 1.08349620.
That builder is unchanged. The earlier graded-ring builders and Step 4
experiment remain available as audit history; their alternative parent
topologies are not the approved parent for this gate.

Run only the geometry gate:

```matlab
addpath(genpath(pwd));
G = main_step4c_preview_literal_crack_cut();
```

For headless artifact export:

```matlab
G = main_step4c_preview_literal_crack_cut( ...
    'Visible','off','OutputDir',fullfile(pwd,'geometry_preview'));
```

The driver produces a full mesh, a negative-x seam zoom, and an inner-boundary
crack-entry zoom. Blue circles and orange crosses identify coincident upper
and lower T6 face nodes without displacing either face. Split children are
colored; their parent triangle IDs remain available in `G.info.cut`.
With `OutputDir`, it writes three PNG figures, a MAT file containing the
parent and final meshes, and a JSON audit summary. With no output directory,
it returns the data and figure handles without writing files.

`build_literal_ring_crack_cut_mesh.m` cuts x2=0, x1<0 in the completed parent
T3 mesh. Edge intersections are shared before separate upper/lower IDs are
created. Every new point lies on an original triangle edge; points on a
staggered-ring chord are not projected onto a circle. Only intersected
triangles are subdivided. Original on-seam vertices have their tiny
`sin(pi)` residual snapped to exact zero in the cut mesh only; the returned
parent preserves the original coordinates. Unsplit triangles touching the
lower crack face require only substitution of duplicate face IDs.

The builder requires the T3 audit to pass before converting to T6, and then
audits the T6 result. Checks cover unchanged original nodes and connectivity
outside the cut neighborhood, correct parent/child subdivision, positive
areas, cut-neighborhood quality and angles, exact face coordinates,
distinct coincident upper/lower IDs, and the absence of elements or shared
edges connecting across the crack. The quality gates require Q >= 0.70
and minimum angle >= 25 degrees near the cut, where
Q = 4 sqrt(3) area / (a^2+b^2+c^2).

Run the regression checks (including deliberately corrupted meshes) with:

```matlab
report = test_literal_crack_cut();
assert(report.passed);
```

The geometry driver passed in MATLAB R2023a on 2026-09-27: 46 parent
triangles split, 5,750 outside-neighborhood rows unchanged, 3,078 T3 nodes,
5,934 triangles and 12,089 T6 nodes. There are 47 T3 and 93 T6 coincident
upper/lower node pairs. Minimum near-cut quality is 0.751266 and minimum
angle is 30.083984 degrees. All triangle areas are positive, with maximum
relative parent/child area discrepancy 2.18e-16. The local bisections are
expected to have lower quality than the uncut near-equilateral triangles;
the old experimental builder's 0.8 seam threshold does not apply here.
The regression also passed, including rejection of all 15 deliberately
corrupted inputs (seam bridges, welded faces, slivers and invalid T6 data).

No FEM solve or SIF extraction is called by this driver. The SIF experiment
must remain untouched and unrun until this geometry gate passes.
