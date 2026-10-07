# Independent M1 trajectory

## Purpose

This run tests **accumulated trajectory sensitivity** to the alternative M1
exterior mesh family.  It is intentionally different from the preceding
fixed-geometry diagnostic.

The reference path is not reused after the common prescribed first segment.
Instead, the entire directional recurrence is recomputed on M1.

## What makes the run independent

The common physical initialization is retained:

- frozen Stage-I initiation point and local frame;
- prescribed first segment `theta_1=0`;
- fixed increment `Delta a=4 mm`;
- same material, loading, plane strain, and MTS rule.

The numerical trajectory is then independent:

1. the P1 mesh is rebuilt with the M1 exterior law;
2. P1 is physically solved on M1;
3. `KI(P1)` and `KII(P1)` from that M1 solve determine `theta_2`;
4. P2 is therefore generated from the M1 P1 field;
5. every later vertex is generated only from the immediately preceding M1
   physical SIFs and MTS turn.

No accepted reference vertex after P1 is imposed.

This distinction matters: the earlier fixed-geometry comparison established
that the exterior coarsening perturbs SIFs only slightly at P17/P22.  The
present experiment allows those small differences to accumulate recursively
into a genuinely different trajectory.

## M1 mesh family

The local numerical machinery is unchanged:

- reflection-paired current-tip core radius: 3 mm;
- tip size: the accepted `hTip/Delta a=0.00675308135`;
- EDI radii: 0.4 and 2.6 mm;
- FE-nodal EDI weight;
- 16-point quadrature;
- the same COD windows;
- all structural/synthetic/solver gates.

Only the exterior size law differs from the reference:

- far cap: `1.25 Delta a = 5 mm`;
- transition length: `1.00 Delta a = 4 mm`;
- far slope: `0.15`;
- boundary-metric growth: `0.35`.

## Run protocol

Switch to branch:

```powershell
git fetch origin
git switch m1-independent-trajectory
git pull
```

From the MATLAB repository root:

### Preflight — no physical solve

```matlab
addpath(genpath(pwd));
close all
M0 = main_m1_independent_trajectory();
```

This qualifies the M1 P1 geometry only.  It does not solve elasticity.

Proceed only if the P1 structural and prescribed-Williams gates all pass.

### Full independent P1--P23 trajectory

```matlab
M = main_m1_independent_trajectory( ...
    'AllowPhysicalSolves',true, ...
    'MaxSegments',23, ...
    'PlotEachStep',false);
```

The run intentionally stops at accepted physical P23.  It does not attempt
P24.

Expected output:

```text
verification/crack_path/m1_independent_run/
    p1_seed/
    trajectory/
    M1_independent_run_summary.mat
    M1_vs_reference.csv
```

The reference comparison CSV is created **after** the M1 propagation and has
no role in determining the M1 path.

## Interpretation rule

Do not judge the M1 run by whether it reproduces every reference vertex.
Path divergence is the observable.

The principal comparisons are:

- `x(a), y(a)`;
- absolute local direction `theta(a)`;
- `KI(a)`, `KII(a)`, and `KII/KI(a)`;
- MTS turn;
- crack length of the positive mode-mixity maximum;
- bracket/location of the mode-mixity sign change;
- onset and magnitude of the turning reversal.

A trajectory-level conclusion is warranted only after the independent M1
history has been completed and audited.


## MATLAB full-run result — PASS

The independent M1 trajectory was completed through accepted physical
segment P23 with `stopReason=max_segments_reached`.

The run retained the prescribed Stage-I/P1 physical initialization, but P1
was recomputed on M1 and no accepted reference trajectory vertex after P1
was imposed.

### Key trajectory events

The M1 trajectory reproduces the same event sequence as the reference run:

- the positive mode-mixity maximum remains at **P17 (68 mm)**;
- at P17:
  - `KI = 0.786789676082`;
  - `KII = +0.00254832743269`;
  - `KII/KI = +0.00323889281997`;
  - next MTS turn `= -0.371140693296 deg`;
- P21 (84 mm) remains slightly positive:
  - `KII/KI = +8.11529313483e-05`;
  - next MTS turn `= -0.00929944077984 deg`;
- P22 (88 mm) is negative:
  - `KII/KI = -0.00214148164742`;
  - next MTS turn `= +0.245393094791 deg`;
- therefore the mode-mixity sign change and turning-sense reversal remain
  bracketed by **P21--P22**;
- the linear visualization-only local-symmetry estimate from the M1
  `q=KII/KI` values is
  `a_LS = 84.1460481757 mm`;
- P23 remains strongly negative in mode mixity and predicts a larger
  positive next turn:
  - `KI = 1.2443571689`;
  - `KII = -0.00678037468375`;
  - `KII/KI = -0.00544889751367`;
  - next MTS turn `= +0.624354409785 deg`;
  - predicted `theta_24 = -2.7761761404 deg`.

### Comparison with the accepted reference trajectory

At the most sensitive late states, the independently propagated M1 path is
extremely close to the reference:

| state | quantity | reference | M1 | M1-reference |
|---|---:|---:|---:|---:|
| P17 | theta (deg) | -2.44443486845 | -2.44441112407 | +2.374e-5 |
| P17 | KII/KI | +0.003238930803 | +0.003238892820 | -3.798e-8 |
| P21 | theta (deg) | -3.63664264573 | -3.63662420420 | +1.844e-5 |
| P21 | KII/KI | +8.13426737e-5 | +8.11529313e-5 | -1.897e-7 |
| P22 | theta (deg) | -3.64596382938 | -3.64592364498 | +4.018e-5 |
| P22 | KII/KI | -0.002141591607 | -0.002141481647 | +1.100e-7 |
| P23 | theta (deg) | -3.40055813454 | -3.40053055019 | +2.758e-5 |
| P23 | KII/KI | -0.005448833169 | -0.005448897514 | -6.434e-8 |

The M1 local-symmetry interpolation differs from the reference estimate
`84.1463699119 mm` by only about **0.000322 mm**.  This is far below the
4-mm propagation increment and is not interpreted as physical precision.

### Scientific conclusion

This is stronger than the fixed-geometry mesh test.  The alternative
exterior mesh family was allowed to perturb the P1 SIFs, and those
differences were then accumulated recursively through every MTS update.
Nevertheless, the independent M1 calculation reproduces the reference
trajectory, the P17 mode-mixity maximum, the P21--P22 sign-change bracket,
and the onset of positive turning after P22.

Accordingly, the main LEFM trajectory features are robust with respect to
this substantial coarsening of the exterior mesh while the audited local
tip core and EDI extraction region are held fixed.

This result remains a **mesh-family sensitivity result**, not a
crack-increment convergence result.  The next distinct numerical question is
the trajectory obtained with `Delta a = 2 mm`.
