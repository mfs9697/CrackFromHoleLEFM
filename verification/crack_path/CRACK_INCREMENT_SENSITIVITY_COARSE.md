# Crack-increment sensitivity on a coarse dimensionless mesh family

## Scientific question

The crack-increment study uses the three-level family

```text
Delta a = 4 mm, 2 mm, 1 mm.
```

The 4-mm and 2-mm trajectories are completed. The 1-mm level is the final
refinement used to test whether the trajectory, turning density, and late
curvature-reversal location approach stable limits.

The comparison should not use the expensive production meshes.  Instead, both
increments use the same **dimensionless coarse family**, combining the two
already qualified coarsening ideas:

- local structured core: `CoreScale = 2`;
- exterior mesh: M1.

This is a convergence-family test of the whole fixed-increment algorithm.

## Dimensionless numerical family

For all three increments,

```text
hTip / Delta a = 2 * 0.00675308135
r_i  / Delta a = 0.10
r_o  / Delta a = 0.65
r_c  / Delta a = 0.75

M1 far cap / Delta a    = 1.25
M1 transition / Delta a = 1.00
M1 far slope            = 0.15
M1 boundary growth      = 0.35
```

Therefore the absolute scales are:

| quantity | Delta a = 4 mm | Delta a = 2 mm | Delta a = 1 mm |
|---|---:|---:|---:|
| hTip | 0.0540247 mm | 0.0270123 mm | 0.0135062 mm |
| r_i | 0.4 mm | 0.2 mm | 0.1 mm |
| r_o | 2.6 mm | 1.3 mm | 0.65 mm |
| r_c | 3.0 mm | 1.5 mm | 0.75 mm |
| M1 far cap | 5.0 mm | 2.5 mm | 1.25 mm |
| M1 transition | 4.0 mm | 2.0 mm | 1.0 mm |

The same nondimensional COD windows and qualification gates are retained.

Because the physical mesh scales with `Delta a`, this experiment should be
described as fixed-increment **scheme convergence**, not as a pure
single-parameter experiment with all absolute FE lengths held fixed.

## Common physical initialization

Stage I remains the accepted hole-only solution.  The initiation point and
frozen material frame are unchanged.

Only the reserved Stage-II increment stored in the frozen summary is
overridden in memory.  No Stage-I physical quantity is recomputed or changed.

For each increment:

1. build and qualify its own P1 first-segment mesh on the M1 + CoreScale=2 family;
2. physically solve P1;
3. use that run's own P1 SIFs to determine theta_2;
4. propagate recursively with no reference-path vertices imposed;
5. target 92 mm total crack length:
   - 23 segments for 4 mm;
   - 46 segments for 2 mm;
   - 92 segments for 1 mm.

## Guarded workflow

Qualification only:

```matlab
addpath(genpath(pwd));
close all

Q4 = main_increment_sensitivity_coarse('IncrementMM',4);
Q2 = main_increment_sensitivity_coarse('IncrementMM',2);
Q1 = main_increment_sensitivity_coarse('IncrementMM',1);

disp(Q4.P1Qualification);
disp(Q2.P1Qualification);
disp(Q1.P1Qualification);
```

After both P1 preflights pass, run the 4-mm coarse trajectory:

```matlab
M4 = main_increment_sensitivity_coarse( ...
    'IncrementMM',4, ...
    'TargetLengthMM',92, ...
    'AllowPhysicalSolves',true);
```

Then run the 2-mm trajectory:

```matlab
M2 = main_increment_sensitivity_coarse( ...
    'IncrementMM',2, ...
    'TargetLengthMM',92, ...
    'AllowPhysicalSolves',true);
```

The final refinement is:

```matlab
M1 = main_increment_sensitivity_coarse( ...
    'IncrementMM',1, ...
    'TargetLengthMM',92, ...
    'AllowPhysicalSolves',true);
```

The 1-mm run uses its own output directory and supports resume from the
saved `path_run_state.mat`.

The driver supports resume from its saved `path_run_state.mat`.

## Comparison

After both runs finish:

```matlab
Sda = plot_increment_sensitivity_coarse_publication();
disp(Sda);
```

The comparison uses exact common crack lengths

```text
4, 8, 12, ..., 92 mm
```

for the coordinate and direction differences. No trajectory interpolation is
used there. The publication convergence panel uses the native turning density
`Delta theta / Delta a` for all three levels. Raw `q_K=KII/KI` is retained as
a diagnostic, together with `q_K/Delta a`, its positive maximum, and its
linearly interpolated zero crossing.

The publication diagnostics are:

- maximum vertical-coordinate difference at matched crack lengths;
- maximum Euclidean corresponding-tip separation;
- maximum absolute direction difference;
- maximum absolute mode-mixity difference at matched lengths;
- crack length of the positive mode-mixity maximum for each increment;
- P21/P22-type sign-change bracket for each native history;
- linear visualization-only local-symmetry estimate for each increment.

## Interpretation rule

The three-level 4/2/1-mm family is intended to test whether successive path
differences contract and whether the turning-density and curvature-reversal
observables approach stable limits **within this scaled coarse numerical
family**.

Even three levels should be interpreted as a numerical convergence study,
not as proof of a formal asymptotic order. It complements, rather than
replaces, the separate exterior-mesh and local-tip-resolution tests.
