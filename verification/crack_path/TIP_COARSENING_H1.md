# Fixed-geometry near-tip coarsening: H1 = 2 hTip

This diagnostic introduces the first coarse member of the deterministic
structured crack-tip family while holding the accepted crack geometry fixed.

## H1 definition

Reference H0:

- `CoreScale = 1`;
- `hTip = 0.0270123254 mm`.

Coarse H1:

- `CoreScale = 2`;
- `hTip = 0.0540246508 mm`.

Unchanged:

- crack increment: 4 mm;
- `rInner = 0.4 mm`;
- `rOuter = 2.6 mm`;
- `rCore = 3.0 mm`;
- reference exterior mesh law;
- material, loading and plane-strain setting;
- FE-nodal EDI formulation and quadrature;
- MTS criterion;
- physical solver and acceptance gates.

The default fixed states are P21 and P22 because they bracket the accepted
local-symmetry crossing.

## Deterministic H1 fingerprints

The exact structured-core construction predicts the following H1 fingerprints
before full-domain qualification:

- paired-core T3 elements: 3,318;
- production EDI participating elements: 2,976;
- skip-constant-q mesh-audit subset: 2,700;
- native COD counts in the four established windows: 19/28/23/18.

These are hard qualification expectations, not values learned from a passing
run. Full-domain geometry, topology, quality, support-containment and
prescribed-Williams gates must still pass.

For reference, the historical Step62 H0 distinction is:

- production EDI participation: 11,316 elements;
- skip-constant-q primary mesh-audit subset: 10,278 elements.

## Phase 1: qualification only

Run:

```matlab
addpath(genpath(pwd));

H1 = main_fixed_geometry_tip_coarsening_diagnostic();

disp(H1.qualification);
```

The default is safe: no physical FEM solve is performed.

Expected H1 tip size:

```text
0.0540246508 mm
```

Do not authorize physical solves unless both P21 and P22 pass every structural
and synthetic-Williams gate.

## Phase 2: guarded physical solves

Only after Phase 1 passes:

```matlab
H1 = main_fixed_geometry_tip_coarsening_diagnostic( ...
    'AllowPhysicalSolves', true);

disp(H1.qualification);
disp(H1.physical);
```

Results are written separately under:

```text
verification/crack_path/fixed_geometry_tip_coarsening_H1/
```

so the existing L1 h/2 refinement evidence is not overwritten.

## Planned mesh figure

After H1 is qualified, retain its saved mesh preview/candidate for the planned
common-window comparison:

- H2: `hTip = 0.1080493016 mm` (future, not yet admitted);
- H1: `hTip = 0.0540246508 mm`;
- H0: `hTip = 0.0270123254 mm`.

H2 must receive its own fingerprint and qualification step before use.
