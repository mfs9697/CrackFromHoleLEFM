# Stage III-C — genuine kinked-polyline qualification

## Purpose

Stage III-C is the final geometry/topology qualification required before the
production incremental crack-path loop.

The production initialization is now:

- theta_1 = 0 deg by prescription at the frozen Stage-I hole-initiation point;
- the SIFs at the end of that 4-mm first segment determine theta_2 by MTS.

Using the accepted Stage-II theta_1=0 physical result,

    KI  = 0.366479612185 MPa sqrt(m)
    KII = 4.19826648316e-06 MPa sqrt(m)

the existing MTS routine gives

    theta_2 = Delta theta_2 = -0.0013127221416 deg.

Stage III-C does not perform a physical solve. It qualifies the mesh machinery
needed to represent that non-collinear second segment exactly.

## Why a visible control is required first

The physical second-leg turn is extremely small. Over a 4-mm increment it
changes the second-tip transverse position by only about 0.092 micrometers,
and the kink vertex lies only about 0.046 micrometers away from the
mouth-to-tip chord.

A geometry implementation that accidentally straightened the crack could
therefore look plausible.

Stage III-C must first be run with a deliberately visible nonphysical control:

    theta_2 = -5 deg.

Only after that case passes every anti-straightening, topology, core, EDI, and
synthetic-field gate should the same code be run with the physical MTS value

    theta_2 = -0.0013127221416 deg.

## True polyline carrier

The historical appended-hole builder was straight-only even when its API
accepted a polyline. Stage III-C introduces a separate audited route:

- build_appended_hole_polyline_loop.m
  - finite-width upper/lower carrier faces follow every crack segment;
  - an interior kink uses a mitered offset;
  - both faces taper to the common sharp current tip.

- build_domain_hole_true_polyline.m
  - merges that genuine polyline pencil with the circular hole.

- identify_polyline_pencil_edge_sets.m
  - recovers every PDE geometry edge belonging to each crack face.

- collapse_polyline_pencil_faces_to_midline.m
  - collapses all upper/lower face edges onto the complete supplied midline;
  - keeps upper/lower node IDs and topology distinct;
  - retains every historical path vertex.

## Polyline exterior

build_stage3c_polyline_exterior.m generalizes the closed-audit C03 exterior
from a negative-x-axis crack to an explicit retained crack polyline.

In the current-tip frame:

- the last crack segment is aligned with +x toward the current tip;
- the previous tip is exactly at [-Delta a,0];
- the paired core intersects the last segment at [-rCore,0];
- the retained exterior crack is
      mouth -> all previous path vertices -> [-rCore,0].

The retained crack itself is a constrained Delaunay boundary. After exterior
triangulation, every retained crack node is duplicated and triangles on the
lower/right side receive the duplicated IDs.

## Last-segment asymptotic frame

native_COD_polyline_audit.m wraps the qualified historical COD extractor.

For a polyline crack it replaces only the extraction descriptor

    Pmid = full crack path

by

    Pmid = [previous tip; current tip]

before calling native_COD_audit.

The displacement field, topology, face IDs, COD formula, and fitting method
are unchanged.

This locks the production rule that SIF/COD extraction at step k uses the
LAST crack segment, never the mouth-to-tip chord.

SIF_LEFM_interaction_EDI already uses the last segment of the supplied
polyline and therefore requires no analogous change.

## New-tip scales

The current increment remains Delta a = 4 mm:

- hTip = 0.0270123254 mm;
- rInner = 0.4 mm;
- rOuter = 2.6 mm;
- rCore = 3.0 mm;
- transition = 4.0 mm;
- far cap = 2.5 mm.

The previous tip is 4 mm from the current tip and must remain outside both the
paired core and the EDI annulus.

## Qualification gates

In addition to all accepted Stage III-A core, exterior, boundary, Jacobian,
sampling, and prescribed-Williams gates, Stage III-C requires:

- genuinely non-collinear second leg;
- exact supplied kink angle;
- multiple PDE geometry edges on both finite-width crack faces;
- exact retained exterior polyline;
- visible kink vertex off the mouth-to-tip chord;
- separate upper/lower carrier nodes at the kink;
- separate upper/lower final crack nodes at the kink;
- exact previous-tip distance of one increment;
- previous tip outside rCore and rOuter;
- last-segment COD frame;
- exact 38/55/44/34 near-tip sampling fingerprint;
- complete EDI support inside the untouched paired core;
- no exterior element in EDI support;
- prescribed pure-I, pure-II, and tiny-mixed Williams recovery.

No physical solve is performed.

## Run 1 — visible control

On branch

    stage3c-kinked-polyline-qualification

run:

    addpath(genpath(pwd));

    Q3c5 = main_stage3c_kinked_two_leg_qualification( ...
        'FrozenState', R0, ...
        'SecondLegLength', 0.004, ...
        'Theta2Deg', -5, ...
        'Plot', true);

    disp(Q3c5.summary);
    disp(Q3c5.gates);
    disp(Q3c5.synthetic);
    disp(Q3c5.syntheticGates);
    disp(Q3c5.sampleCounts);

Do not proceed to the physical small angle unless this visible control passes.

## Run 2 — physical second segment

Derive theta_2 from the accepted theta_1=0 Stage-II physical SIFs rather than
typing a rounded angle:

    [~,theta2Deg] = kink_angle_LEFM_MTS( ...
        Rphys.EDI.KI_unit, Rphys.EDI.KII_unit);

Then run:

    Q3cPhys = main_stage3c_kinked_two_leg_qualification( ...
        'FrozenState', R0, ...
        'SecondLegLength', 0.004, ...
        'Theta2Deg', theta2Deg, ...
        'Plot', true);

If both the visible control and the physical case pass, the mesh machinery is
qualified for the first genuine incremental propagation state.

## Next stage

After Stage III-C passes, Stage III-D will perform one guarded physical solve
of the actual two-leg crack with prescribed theta_2 from MTS.

That solve will determine

    KI(P2), KII(P2)

and then

    Delta theta_3 = MTS(KI(P2),KII(P2)),
    theta_3 = theta_2 + Delta theta_3.

No angle sweep is part of the production incremental path.
