# SIF extraction on asymmetric FEM meshes: audit log

## Purpose

The audit asks a narrow question: how much does the historical
mirror-based mode separation distort (K_I) and (K_{II}) when the
discrete FEM field is not mirror symmetric, and does an interaction
equivalent-domain integral provide a stable alternative on arbitrary
asymmetric meshes?

No attempt will be made to preserve or justify old numerical values.  Old
results are retained as historical/control values until independently
verified.

## Source-code audit

Two repositories are involved:

- **mfs9697/Crack-Path** supplies the established Crack-Path geometry,
  T6 FEM, loading, and the historical SIF workflow.
- **mfs9697/CrackFromHoleLEFM** contains both the historical
  `SIF_LEFM_circle2` family and the prototype
  `SIF_LEFM_interaction_EDI`.

The old `SIF_LEFM_circle2.m` present in both repositories implements the
same conceptual procedure: mirrored contour points (P,Q) are evaluated
independently in the discrete T6 field, symmetric/antisymmetric fields are
formed, separate (J_I,J_{II}) values are integrated, and those are
converted to SIFs.

The key numerical vulnerability is therefore not continuum mode separation
itself. It is the assumption that the two independently interpolated FEM
fields at mirrored locations have the parity properties of the continuum
solution. An asymmetric triangulation need not satisfy that assumption.

## Step 1 implemented on branch `sif-asymmetric-mesh-audit`

The first code-extraction step deliberately leaves both SIF extractors
unchanged.

A verification-only pipeline has been added:

1. `cfg_crack_path_two_leg_control.m` defines a two-leg elastic control
   case using the current non-perforated Crack-Path scales.
2. `build_crack_path_polyline_LEFM_mesh.m` adapts the pencil-channel
   geometry from Crack-Path and collapses *all* crack segments to a
   traction-free polyline while retaining distinct crack-face topology.
3. `solve_crack_path_polyline_field.m` solves one T6 displacement field
   with the Crack-Path loading and rigid-body constraints.
4. `run_crack_path_old_vs_edi.m` sends that exact same FEM field to
   `SIF_LEFM_circle2_debug` and `SIF_LEFM_interaction_EDI`.

The helper `subdivide_last_leg.m`, which was used by existing geometry
code but absent from CrackFromHoleLEFM, has also been restored from
Crack-Path.

Two repository-level issues found during extraction are intentionally kept
outside the production path for now:

- `kinking_LEFM_1leg.m` in CrackFromHoleLEFM calls
  `geom_pencil_1leg.m`, but that geometry file is absent from the target
  repository;
- the existing `mesh_pencil_domain.m` computes its `Hmin` from an
  already subdivided last segment and then divides by `ncoh` again,
  effectively introducing an extra factor of `ncoh` in the local target
  size. The verification builder therefore implements the Crack-Path
  `h_last = L_last/ncoh` rule directly instead of using that routine.

Neither production issue is changed in this PR because the first goal is to
isolate the SIF comparison from unrelated refactoring.

## What this step does *not* establish

This first control case is **not yet the historical published two-segment
benchmark**. The exact mapping of the older publication geometry and its
plane-state convention is still to be closed before a published numerical
baseline is declared reproduced.

Likewise, `SIF_LEFM_interaction_EDI.m` remains a prototype. At this stage
old-vs-EDI differences are method-to-method differences, not EDI-based
errors.

## Verification gates

The next gates are:

- independently verify the EDI normalization and mode-II sign using exact
  analytical crack-tip fields;
- establish a deliberately mirror-symmetric FEM mesh and require old/EDI
  agreement there;
- construct controlled asymmetric mesh families while keeping geometry,
  material, loading, and crack path fixed;
- sweep old contour radius and EDI annulus radii;
- correlate old-method changes with explicit mesh-asymmetry measures;
- only after those gates, decide which extractor is suitable for production
  crack-path calculations.

Incremental crack growth is outside the scope of this audit until these
gates are closed.


## First same-field numerical result

The first local MATLAB run of `main_step1_same_field_compare` produced:

| Method | KI | KII |
|---|---:|---:|
| old mirror/J | 0.33973 | 0.0033961 |
| interaction EDI | 0.68007 | 0.0073326 |

with the same FEM field, old contour radius (8.0\times10^{-4}), and EDI
annulus ([8.0\times10^{-5},8.0\times10^{-4}]).

The mode-I ratio is approximately (0.68007/0.33973 \approx 2.0018). This
is a strong diagnostic signature of a possible factor-of-two normalization
issue in the EDI conversion, but it is not by itself sufficient to alter the
production implementation. The mode-II value is small and correspondingly
more sensitive; its ratio is not used as the normalization diagnostic.

## Step 2: synthetic Williams-field normalization test

The branch now includes:

- `validate_EDI_Williams_fields.m`;
- `main_step2_validate_edi_normalization.m`.

This test constructs an independent polar annulus with duplicated crack-face
nodes, prescribes exact leading-order Williams mode-I, mode-II, and mixed-mode
displacements, and asks the unchanged EDI implementation to recover the
imposed SIFs. It reports both raw/input and (0.5\,\mathrm{raw}/\mathrm{input})
ratios across mesh refinements.

Decision rule:

- if raw/input tends to 1, the current EDI normalization is consistent;
- if raw/input tends to 2 while (0.5\,\mathrm{raw}/\mathrm{input}) tends
  to 1, the conversion requires an explicit factor (1/2);
- pure-mode cross leakage and the sign of recovered mode II are checked at
  the same time.

No EDI source-code correction has been made before this gate is run.


## Step 2 result: normalization error identified and corrected

The independent exact-field test was run locally on 2026-09-27.

For the coarse polar annulus mesh (Nr=8, Nth=64), the pre-correction
implementation returned approximately:

- pure mode I: KI_recovered/KI_input = 1.9722;
- pure mode II: KII_recovered/KII_input = 1.9722;
- mixed mode: both modal recovery ratios = 1.9722.

For the refined mesh (Nr=16, Nth=128), the ratios were:

- pure mode I: 2.0088;
- pure mode II: 2.0088;
- mixed mode: 2.0088 for both modes.

The non-imposed modal components converged toward zero, and pure mode II
was recovered with the same sign as the imposed field.

This establishes that the interaction integral itself was producing the
standard interaction quantity

\[
I^{(1,2)} = \frac{2}{E'}\left(
K_I^{(1)}K_I^{(2)} + K_{II}^{(1)}K_{II}^{(2)}
\right),
\]

whereas the conversion in `SIF_LEFM_interaction_EDI.m` had omitted the
factor (1/2).

The conversion has therefore been corrected to

\[
K = \frac{E'}{2K^{aux}} I.
\]

This correction is based on the independent analytical Williams-field
experiment, not on fitting the EDI result to the historical mirror/J value.

The Williams-field driver has now been converted into a regression test:
after the correction, recovered/input must approach one, pure-mode cross
leakage must approach zero, and mode-II sign must remain positive.


## Post-correction Williams-field regression

The corrected EDI implementation was rerun locally on 2026-09-27.

Fine mesh (Nr=16, Nth=128):

- pure mode I: KI_recovered/KI_input = 1.0044;
- pure mode II: KII_recovered/KII_input = 1.0044;
- mixed mode: both modal recovery ratios = 1.0044;
- pure-I cross leakage: KII = -2.2813e-7;
- pure-II cross leakage: KI = 6.8440e-7.

The coarse mesh (Nr=8, Nth=64) recovered approximately 0.9861 for both
modes, so the refinement trend brackets unity and the fine-mesh error is
about 0.44%.

The post-fix regression therefore passes the configured 2% tolerance.
Mode-II sign is confirmed independently by the pure-II exact field.

Status of the EDI normalization/sign gate: **PASSED**.
