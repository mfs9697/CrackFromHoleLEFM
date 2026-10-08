# First-segment local-symmetry study

## Research status

This directory starts a **new research sequence** after closure of the
asymmetric SIF audit. It is intentionally not numbered Step71.

The closed audit remains on `sif-asymmetric-mesh-audit`; its exact pre-cleanup
snapshot is preserved on `archive/sif-audit-closed-2026-10-05`.

The new scientific objective is:

\[
\text{hole-only initiation}
\;\longrightarrow\;
K_{II}(\theta_1)
\;\longrightarrow\;
K_{II}(\theta_1^*)=0
\;\longrightarrow\;
\text{subsequent MTS crack path}.
\]

No first-segment direction has yet been selected.

## Which starting geometry is authoritative?

Three configurations occur in repository history and must not be conflated.

1. **Historical early prototype:** hole center `[0.150,-0.060]` m. This is
   the configuration described in the older `Documents/CrPathLEFM.tex`.
2. **Centered verification configuration:** the generic
   `cfg_hole_initiation.m` currently places the hole at `[0.5*A,0]`.
   This was introduced for centered symmetry/parity benchmarking.
3. **Accepted asymmetric single-crack benchmark:** hole center
   `[0.17,-0.02]` m, radius `0.03` m in the `0.30 x 0.20` m plate.
   This is the configuration validated in audit Steps 19 onward and is the
   physical starting point for the new crack-path study.

The new wrapper `cfg_first_segment_asymmetric.m` freezes item 3 explicitly so
future changes to the generic centered configuration cannot silently redefine
the crack-path problem.

## Frozen Stage-I configuration

The Stage-I baseline uses:

- plate: `A=0.30 m`, `B=0.10 m`;
- circular hole: center `[0.17,-0.02] m`, radius `0.03 m`;
- `Npoly=480`;
- plane strain, `E=210e3`, `nu=0.30`;
- unit remote tension in `y`;
- minimal anchoring;
- strength threshold `sigma_c=300`;
- production boundary estimator `boundary_extrapolated_t6`;
- material-side offsets `[0.05,0.10,0.25] h_hole`;
- linear radial extrapolation to the hole boundary;
- mesh-scaled local quadratic angular fit with factor `c=3`.

The initial short-crack length reserved for the next stage is

\[
a_0=4\ {\rm mm}.
\]

This value does not affect the Stage-I hole-only solve.

## Historical fingerprint

The accepted fine Stage-I result from audit Steps 19/22/23 is

\[
\phi_*=-1.5606127781^\circ,
\]

with

\[
\mathbf{x}_*
=
[0.19998887220,\,-0.020817033905]\ {\rm m}.
\]

For the circular hole this corresponds to

\[
\mathbf n_{\rm mat}
\approx
[0.9996290732,\,-0.0272344635],
\]

\[
\mathbf t
\approx
[0.0272344635,\,0.9996290732].
\]

The new Stage-I run uses these numbers only as a broad configuration
fingerprint; the newly computed values are saved as the authoritative starting
state.

The same audit found a physically separate secondary tensile maximum about
`1.34%` below the dominant maximum. Therefore the present study is a
single-crack first-initiation example; a later multi-crack problem would need
fresh stress redistribution after the first crack grows.

## Stage-I freeze run

The guarded driver is:

`verification/crack_path/main_stage1_freeze_starting_state.m`

Run:

```matlab
addpath(genpath(pwd));

R0 = main_stage1_freeze_starting_state( ...
    'AllowSolve', true);

disp(R0.summary);
disp(R0.peaks);
disp(R0.gates);
```

The driver performs exactly one **hole-only** FEM solve. It does not build a
crack or evaluate an SIF.

Persistent compact outputs are written to:

- `verification/crack_path/stage1_starting_state.mat`;
- `verification/crack_path/stage1_starting_state.csv`;
- `verification/crack_path/stage1_boundary_stress.csv`.

## Gate after Stage I

If the Stage-I freeze passes, the next task is **not** an angle sweep yet.

The next task is to construct the `theta=0` 4-mm trial crack and qualify an
interaction-EDI domain whose complete support remains inside a controlled
reflection-paired tip/core region and clear of the hole mouth.

Only after that 4-mm EDI geometry/extraction gate passes will we compute the
physical curve

\[
\theta \mapsto K_{II}(\theta)
\]

and search for the first-segment local-symmetry root.


## Near-tip h/2 refinement

The controlled fixed-geometry L0/L1 near-tip refinement protocol is documented
in [TIP_REFINEMENT_H2.md](TIP_REFINEMENT_H2.md). It halves the structured tip
scale while keeping the accepted crack geometry, EDI radii, paired-core radius,
reference exterior law, physics, and solver gates fixed.


## Near-tip coarsening sensitivity

See [TIP_COARSENING_H1.md](TIP_COARSENING_H1.md) for the guarded fixed-geometry H1 = 2 hTip qualification experiment around the P21--P22 local-symmetry bracket.


## Independent local tip-resolution trajectory

See [TIP_2H0_INDEPENDENT_TRAJECTORY.md](TIP_2H0_INDEPENDENT_TRAJECTORY.md) for the guarded independent `2h0` trajectory experiment with the reference exterior mesh law.


## Symmetric path-stability control

See [SYMMETRIC_PATH_STABILITY.md](SYMMETRIC_PATH_STABILITY.md) for the centered-hole five-probe one-step stability test of the EDI--MTS recurrence.
