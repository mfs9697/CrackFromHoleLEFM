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
