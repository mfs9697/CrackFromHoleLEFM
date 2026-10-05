# Stage III-D — first incremental physical kinked solve

## Purpose

Stage III-D is the first **production incremental crack-path solve**.

The crack geometry is already fixed before the solve:

- first leg: Delta a = 4 mm at theta_1 = 0 deg;
- second leg: Delta a = 4 mm at theta_2 determined from the accepted P1 SIFs.

The accepted P1 physical result is

    KI(P1)  = 0.366479612185 MPa sqrt(m)
    KII(P1) = 4.19826648316e-6 MPa sqrt(m)

and the production MTS helper gives

    theta_2 = -0.00131272214162 deg.

Stage III-D does **not** vary theta_2. It performs one physical LEFM solve of
the already-qualified Stage III-C candidate and predicts the next direction.

The production recurrence is now

    fixed crack through P_k
        -> solve for KI(P_k), KII(P_k)
        -> Delta theta_{k+1} by MTS
        -> theta_{k+1} = theta_k + Delta theta_{k+1}.

For Stage III-D, k=2.

## Qualified input state

Stage III-C has already passed twice:

1. visible nonphysical qualification at theta_2 = -5 deg;
2. physical candidate at theta_2 = -0.00131272214162 deg.

The physical Stage III-C fingerprint is hard-gated:

- T3 nodes/elements: 25,325 / 49,637;
- T6 nodes: 100,287;
- paired-core T3 elements: 12,678;
- literal EDI support: 11,316 elements;
- polyline-safe native COD face sets: 138 / 138;
- native COD windows: 38 / 55 / 44 / 34 points.

All Stage III-C structural and synthetic gates must remain true.

## Physical solver

The solver architecture is inherited unchanged from the accepted Stage III-B
physical solve:

- assemble the unclamped symmetric global stiffness matrix;
- enforce only the frozen minimal anchors through free-DOF elimination;
- exact unit remote-y traction;
- symamd permutation;
- parameter-free symmetric Gauss-Seidel preconditioner;
- PCG tolerance 1e-10;
- maximum 5000 PCG iterations;
- true free-DOF residual gate 5e-10;
- exact constrained DOFs;
- no direct backslash solve.

Exactly one physical linear state is authorized.

## Checkpoint rule

The physical field is checkpointed **before any COD, EDI, or MTS
postprocessing**.

The checkpoint is accepted only when all of the following match the current
qualified candidate:

- Stage III-D checkpoint label;
- T3/T6 fingerprints and connectivity;
- T3 coordinates;
- complete crack polyline P0 -> P1 -> P2;
- current increment;
- prescribed theta_2;
- accepted P1 KI/KII provenance;
- single physical state;
- no angle sweep;
- no third leg generated;
- qualified SGS-PCG solver.

A valid checkpoint is reused without a second physical solve.

Default checkpoint:

    verification/crack_path/
    stage3d_kinked_theta_m0p00131272deg_physical_solved.mat

## Postprocessing

Postprocessing uses the saved physical field only.

### COD

native_COD_polyline_audit.m is mandatory.

It builds the crack-tip extraction frame from the **last segment P1 -> P2**,
never from the mouth-to-tip chord.

The four fixed windows are

    [0.04,0.20] Delta a
    [0.04,0.30] Delta a
    [0.08,0.30] Delta a
    [0.12,0.30] Delta a

with linear and quadratic fits used as diagnostics.

### Interaction EDI

Exactly one physical EDI is run with

    r_inner = 0.10 Delta a = 0.4 mm
    r_outer = 0.65 Delta a = 2.6 mm

using the unchanged 16-point FE-nodal-q production interaction integral.

SIF_LEFM_interaction_EDI.m already constructs its frame from the last crack
segment.

## MTS output

From the measured P2 SIFs,

    KI^(2), KII^(2),

Stage III-D computes

    Delta theta_3 = MTS(KI^(2),KII^(2))

relative to leg 2, followed by

    theta_3 = theta_2 + Delta theta_3

in the frozen local Stage-I frame.

The global direction is also reported.

**No third leg is generated.**

The result is therefore a prediction to be consumed by the next mesh
qualification/propagation stage.

## Run

Checkout

    stage3d-kinked-physical-solve

and run, while R0 and Q3cPhys remain in memory:

    addpath(genpath(pwd));

    R3d = main_stage3d_kinked_two_leg_physical_solve( ...
        'FrozenState', R0, ...
        'Candidate', Q3cPhys.candidate, ...
        'AllowSolve', true);

    disp(R3d.summary);
    disp(R3d.fitTable);
    disp(R3d.EDI);
    disp(R3d.prediction);
    disp(R3d.solverInfo);
    disp(R3d.gates);

If Q3cPhys is no longer in memory, the driver will by default load

    verification/crack_path/
    stage3c_kinked_theta_m0p00131272deg_candidate_T3.mat

provided that file exists locally.

## Acceptance

Stage III-D passes only if all solver, provenance, COD, EDI, and MTS gates
are true.

There is deliberately **no target value** for the P2 KII/KI ratio or for
theta_3. The previously computed straight 8-mm Stage III-B state may be used
afterward as an independent comparison because the physical kink is small,
but it is not a calibration target and is not an acceptance gate.

## Next step

After Stage III-D passes, theta_3 becomes prescribed geometry for the third
4-mm segment. The next stage should construct and qualify

    P0 -> P1 -> P2 -> P3

with the supplied theta_3, before any physical solve at P3.
