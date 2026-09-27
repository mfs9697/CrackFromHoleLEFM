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
