# Stage III-A — two-leg new-tip core qualification

## Purpose

Stage III-A is the first transition from first-segment selection to incremental
crack propagation.

No physical FEM solve is performed in this stage.

The accepted first segment is frozen as

- length: 4 mm;
- direction: theta_1 = 0 deg relative to the frozen Stage-I material normal.

The numerical local-symmetry interpolation from Stage II is approximately
-0.0012 deg, but this is treated as effectively zero at the present spatial
resolution.

A provisional second segment is then appended with

- length: 4 mm;
- relative turn: Delta theta_2 = 0 deg.

Thus the Stage III-A baseline is a two-leg path with path vertices

    P0 = hole mouth,
    P1 = P0 + 4 mm * e1,
    P2 = P1 + 4 mm * e1.

Although the two legs are collinear in this baseline, P1 is retained as an
explicit crack-path vertex. The Stage III-A exterior builder forces this prior
tip to remain an exact duplicated crack-face node after exterior resampling.

## New-tip scaling

The crack-tip/core design is scaled by the current propagation increment
(4 mm), not by the total 8-mm crack length:

- hTip / Delta a = 0.00675308135;
- hTip = 0.0270123254 mm;
- EDI inner radius = 0.10 Delta a = 0.4 mm;
- EDI outer radius = 0.65 Delta a = 2.6 mm;
- paired-core radius = 0.75 Delta a = 3.0 mm;
- exterior transition length = 1.0 Delta a = 4.0 mm;
- far-field cap = 0.625 Delta a = 2.5 mm.

The previous tip P1 is 4 mm behind the new tip P2 and must therefore lie
outside both the 3-mm paired core and the 2.6-mm EDI annulus.

## Qualification gates

The driver requires:

- frozen Stage-I state passed;
- 480-point frozen hole representation;
- exact mouth, prior-tip, and new-tip positions;
- exact 4-mm first and second leg lengths;
- first leg along the frozen material normal;
- straight provisional continuation;
- prior tip represented by distinct upper/lower crack-face node IDs;
- prior tip outside the new paired core and EDI annulus;
- physical boundary outside the new paired core and EDI annulus;
- unchanged 12,678-element paired structured core;
- exact T3/T6 reflection pairing;
- positive T3 areas and T6 Jacobians;
- valid edge incidence and no duplicate triangles;
- preserved physical boundary geometry;
- complete EDI support inside the untouched paired core;
- no exterior element in the EDI support;
- exact native COD sampling fingerprint 38/55/44/34;
- minimum angle >= 20 deg;
- maximum adjacent size ratio <= 1.8.

After structural qualification only, prescribed Williams fields are replayed on
the full assembled T6 mesh:

- KI=1, KII=0;
- KI=0, KII=1;
- KI=1, KII=1e-4.

No physical SIF is calculated.

## Local run

On branch:

    stage3a-two-leg-tip-core-qualification

run:

    addpath(genpath(pwd));

    Q3 = main_stage3a_two_leg_tip_core_qualification( ...
        'FrozenState', R0, ...
        'SecondLegLength', 0.004, ...
        'DeltaTheta2Deg', 0, ...
        'Plot', true);

    disp(Q3.summary);
    disp(Q3.gates);
    disp(Q3.synthetic);
    disp(Q3.syntheticGates);
    disp(Q3.sampleCounts);

If all gates pass, the exact T3 candidate is saved as

    verification/crack_path/stage3a_two_leg_straight_candidate_T3.mat

and is eligible for a separate Stage III-B one-angle physical solve.

## Scope limitation

Stage III-A deliberately does not yet qualify a non-collinear second leg.

The current audited exterior construction splits a straight retained crack in
the current-tip frame. A nonzero MTS turn requires a subsequent polyline
exterior generalization in which the complete kinked crack is preserved as
distinct upper/lower traction-free faces while the new-tip paired core remains
unchanged.
