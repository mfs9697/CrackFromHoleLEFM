# Fixed-geometry coarse-exterior diagnostic

## Purpose

This diagnostic asks a narrow LEFM question before any alternative crack
trajectory is propagated:

> For the **same accepted crack geometry**, how much do the current-tip SIFs
> and MTS turn change when only the exterior mesh-density law is made coarser?

The comparison is performed at:

- **P17**, the accepted positive maximum of `KII/KI`;
- **P22**, the first accepted state with negative `KII/KI`.

This is intentionally not yet a trajectory comparison.  Holding the crack
geometry fixed prevents small SIF changes from recursively changing all
subsequent vertices.

## What is unchanged

The trial mesh retains exactly the same:

- accepted Stage-I initiation point and frozen frame;
- accepted P17/P22 crack polylines;
- 4-mm segment increment;
- 3-mm reflection-paired current-tip core;
- audited core topology and tip size;
- EDI radii `ri=0.4 mm`, `ro=2.6 mm`;
- FE-nodal EDI weight and 16-point quadrature;
- material, plane strain, loading, and minimal anchors;
- PCG/SGS solver tolerances;
- COD windows;
- structural, synthetic-Williams, and physical acceptance gates.

No gate is relaxed for the coarser mesh.

## Trial M1 exterior

The accepted reference exterior uses:

- far cap `0.625*DeltaA = 2.5 mm`;
- far slope `0.10`;
- boundary-metric growth `0.25`;
- transition length `1.0*DeltaA = 4 mm`.

The initial **trial M1** changes only the exterior size law:

- far cap `1.25*DeltaA = 5.0 mm`;
- far slope `0.15`;
- boundary-metric growth `0.35`;
- transition length remains `1.0*DeltaA`.

These values are trial parameters, not accepted settings.  They are retained
only if the unchanged qualification gates pass and the resulting mesh is
usefully less dense.

The production defaults of
`qualify_incremental_crack_candidate.m` remain exactly the old values.  New
name-value arguments merely expose exterior-only sensitivity controls.

## Two-phase protocol

### Phase 1: qualification and mesh inspection only

From the repository root:

```matlab
close all
Q = main_fixed_geometry_coarse_exterior_diagnostic();
```

This performs **no physical FEM solve**.

It builds P17 and P22 on M1, runs the unchanged structural and prescribed
Williams-field qualification, writes a qualification comparison CSV, and
opens/saves mesh previews.

Expected output directory:

```text
verification/crack_path/fixed_geometry_coarse_exterior
```

Inspect:

- T3/T6 element counts versus the reference mesh;
- minimum angle;
- maximum neighboring-size ratio;
- unchanged 11316-element EDI support;
- unchanged physical-boundary clearance;
- pure-I, pure-II, and tiny-mixed synthetic recovery;
- the mesh previews, especially P22.

Do not proceed if either qualification fails.

### Phase 2: two guarded fixed-geometry physical solves

Only after Phase 1 is accepted:

```matlab
D = main_fixed_geometry_coarse_exterior_diagnostic( ...
    'AllowPhysicalSolves',true);
```

This solves exactly the same accepted P17 and P22 geometries on M1 and
compares with the preserved reference physical results.

The main comparison table contains:

- `KI`;
- `KII`;
- `KII/KI`;
- next MTS turn;
- PCG iterations and residuals.

Because the geometry is held fixed, these differences measure
**exterior-mesh sensitivity**, not accumulated path divergence.

## Decision after the diagnostic

If M1:

1. passes every unchanged qualification/physical gate;
2. is materially less dense and visually clearer;
3. preserves the physical interpretation at P17 and P22;

then it becomes a candidate exterior family for an **independent propagated
trajectory**.  That later trajectory must be allowed to diverge naturally
from the reference path; it must not be forced through the accepted
reference vertices.

No claim about trajectory robustness is made from this fixed-geometry
diagnostic alone.
