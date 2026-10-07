# Fixed-geometry near-tip h/2 refinement diagnostic

## Purpose

This diagnostic isolates sensitivity to the near-tip finite-element scale after
the exterior-mesh sensitivity study.

The accepted crack geometries at P17 and P22 are held fixed. Their vertices,
reference physical observables, and L0 qualification metrics are read from the
committed audited snapshot `paper/data/evidence_exact.mat`; the diagnostic does
not require the uncommitted `final_clean_run` archive. Only the
reflection-paired structured core is refined from the production L0 family to
the audited L1 family.

## Controlled change

Reference L0:

- core scale = 1;
- hTip/DeltaA = 0.00675308135;
- DeltaA = 4 mm;
- hTip = 0.0270123254 mm;
- paired core T3 = 12,678;
- primary FE-nodal EDI support = 11,316 elements;
- native COD counts = 38/55/44/34.

Refined L1:

- core scale = 0.5;
- hTip = 0.0135061627 mm;
- paired core T3 = 49,518;
- literal radial-support diagnostic = 44,130 elements;
- primary FE-nodal EDI support (skip constant-q) = 40,146 elements;
- native COD counts = 74/108/86/67.

The following remain unchanged:

- accepted crack geometry;
- DeltaA = 4 mm;
- rInner = 0.4 mm;
- rOuter = 2.6 mm;
- rCore = 3.0 mm;
- reference exterior-mesh law;
- material and plane-strain model;
- unit remote-y loading and anchors;
- FE-nodal-q EDI definition and 16-point integration;
- MTS selector;
- SGS-PCG tolerances and residual gates.

## Why P17 and P22

P17 is the accepted state with maximum positive KII/KI.

P22 is the first accepted state with negative KII/KI.

Together they probe the near-tip discretization on both sides of the
late-path sign reversal without changing crack geometry.

## Phase 1: mesh and synthetic qualification only

From the repository root:

    addpath(genpath(pwd));

    Q = main_fixed_geometry_tip_refinement_diagnostic();

    disp(Q.qualification);

No physical FEM solve is performed by default.

Required outcomes include:

- exact L1 core/primary-support/native-sampling fingerprints;
- literal support contains the primary support and remains inside the paired core;
- complete T3/T6 reflection pairing;
- EDI support wholly inside the paired core;
- exterior excluded from the EDI support;
- minimum angle >= 20 deg;
- adjacent-size ratio <= 1.8;
- prescribed Williams recovery gates.

## Phase 2: fixed-geometry physical comparison

Only after Phase 1 is inspected:

    D = main_fixed_geometry_tip_refinement_diagnostic( ...
        'AllowPhysicalSolves', true);

    disp(D.qualification);
    disp(D.physical);

The comparison reports L0 and L1 values of:

- KI;
- KII;
- q = KII/KI;
- MTS turn;
- PCG iterations and residuals.

Because the crack vertices are fixed, any observed difference is a
near-tip discretization effect rather than accumulated path divergence.

## Decision rule

Do not launch a full refined trajectory from the fixed-geometry test alone
until:

1. both L1 meshes pass all structural and synthetic gates;
2. both physical solves pass the unchanged solver/residual gates;
3. the changes in q and MTS turn are interpreted relative to the scale of
   the exterior-mesh M1 sensitivity already measured.

If these conditions are satisfied, the next experiment should be an
independent L1 trajectory. At that stage it is reasonable to combine the L1
tip core with the qualified M1 exterior law to control computational cost,
but that is a separate propagated-path experiment.
