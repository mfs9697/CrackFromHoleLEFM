# Independent 2h0 tip-resolution trajectory

## Purpose

This run tests **accumulated crack-path sensitivity to the local structured
tip resolution** while retaining the reference exterior mesh law.

It is the local-resolution counterpart of the independent M1 exterior-mesh
trajectory experiment.

## Scientific isolation

Unchanged:

- frozen Stage-I initiation point and local frame;
- prescribed first segment, `theta_1 = 0`;
- fixed crack increment, `Delta a = 4 mm`;
- reference exterior size law;
- paired-core radius `rCore = 3 mm`;
- EDI radii `rInner = 0.4 mm`, `rOuter = 2.6 mm`;
- material, loading and plane strain;
- FE-nodal EDI and 16-point quadrature;
- MTS direction rule;
- solver tolerances and physical acceptance gates.

Changed only:

- production tip scale `h0 = 0.0270123254 mm`;
- alternative tip scale `2h0 = 0.0540246508 mm`;
- paired-core T3 fingerprint `12678 -> 3318`;
- production EDI support fingerprint `11316 -> 2976`;
- native COD sampling `38/55/44/34 -> 19/28/23/18`.

The 2h0 core has already passed fixed-geometry qualification and physical
tests at P21 and P22.

## Independence rule

The accepted H0 trajectory is not reused after the common prescribed first
segment.

1. P1 is rebuilt with `CoreScale = 2` and the reference exterior law.
2. P1 is physically solved on that mesh.
3. Its own `KI(P1)` and `KII(P1)` determine `theta_2`.
4. Every later point is generated recursively from the immediately preceding
   2h0 SIFs and MTS turn.
5. Reference vertices are used only for a post-run comparison.

## Phase 1 — qualification-only P1 preflight

From the repository root:

```matlab
addpath(genpath(pwd));
close all

M0 = main_tip2h0_independent_trajectory();

disp(M0.P1Qualification);
```

This performs **no physical solve**.

Proceed only if all P1 structural and prescribed-Williams gates pass and the
summary reports the 2h0 fingerprints.

## Phase 2 — independent P1--P23 trajectory

After accepting the preflight:

```matlab
M = main_tip2h0_independent_trajectory( ...
    'AllowPhysicalSolves', true, ...
    'MaxSegments', 23, ...
    'PlotEachStep', false);

disp(M.Path.stepTable);
```

Output:

```text
verification/crack_path/tip_2h0_independent_run/
    p1_seed/
    trajectory/
    tip2h0_independent_run_summary.mat
    tip2h0_vs_reference.csv
```

The comparison CSV is generated only after propagation and has no influence
on the alternative trajectory.

## Acceptance questions

The principal trajectory-level checks will be:

- whether the positive mode-mixity maximum remains at P17;
- whether P21 remains positive and P22 negative in `KII/KI`;
- whether the turning-sense reversal remains bracketed by P21--P22;
- maximum accumulated vertical tip-coordinate difference from H0;
- maximum Euclidean tip separation;
- maximum absolute direction difference;
- maximum absolute mode-mixity difference;
- the visualization-only interpolated local-symmetry location.

If the complete 2h0 path passes, it can be plotted with the same four-panel
logic as the M1 exterior-mesh sensitivity figure.
