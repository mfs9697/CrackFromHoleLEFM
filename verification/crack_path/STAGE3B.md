# Stage III-B — one physical solve at the second tip

## Purpose

Stage III-B performs exactly one physical LEFM solve of the qualified
Stage III-A straight two-leg crack.

The geometry is already prescribed:

- first leg: 4 mm, theta_1 = 0 deg in the frozen Stage-I local frame;
- second leg: 4 mm, theta_2 = 0 deg;
- total crack length: 8 mm.

Stage III-B does **not** search for theta_2. The two-leg crack is fixed
geometry. The purpose of the physical solve is to obtain KI and KII at the
new tip P2 and use them to predict the next turn by the MTS criterion.

No third segment is generated in this stage.

## New-tip extraction scale

The current propagation increment remains Delta a = 4 mm, even though the
total crack length is 8 mm.

Therefore the qualified new-tip extraction scales remain

- r_inner = 0.10 Delta a = 0.4 mm;
- r_outer = 0.65 Delta a = 2.6 mm;
- paired-core radius = 0.75 Delta a = 3.0 mm.

The COD fitting windows are also normalized by Delta a, not by the total
crack length.

## Solver

The solver is unchanged from the qualified Stage-II physical workflow:

- assemble the unclamped symmetric stiffness matrix;
- restrict to free DOFs;
- check matrix symmetry;
- apply symamd permutation;
- use parameter-free symmetric Gauss-Seidel preconditioning;
- solve with PCG, tolerance 1e-10 and maximum 5000 iterations;
- verify the true residual and exact homogeneous constraints;
- checkpoint the physical displacement field before any COD/EDI/MTS
  postprocessing.

A valid exact checkpoint is reused without another physical solve.

## Candidate fingerprint

The Stage III-B driver is intentionally tied to the exact locally qualified
Stage III-A candidate:

- T3: 25,356 nodes / 49,691 elements;
- T6: 100,403 nodes;
- paired core: 12,678 T3 elements;
- literal EDI support: 11,316 elements;
- native crack-face nodes: 220 / 220;
- current-increment COD sampling: 38 / 55 / 44 / 34.

All Stage III-A structural and prescribed-Williams gates must still be true.

## Physical outputs

The primary physical quantities are

    KI(P2), KII(P2), KII(P2)/KI(P2).

The MTS routine is then called exactly once:

    Delta theta_3 = MTS(KI(P2), KII(P2)).

The reported next direction is

    theta_3 = theta_2 + Delta theta_3.

The driver reports both

- theta_3 in the frozen Stage-I local frame;
- theta_3 in global coordinates.

The eight COD fits also receive MTS turns as diagnostics. The primary
propagation prediction is the EDI-based MTS turn.

No expected sign or magnitude of KII or Delta theta_3 is imposed as a gate.

## Local run

On branch

    stage3b-new-tip-physical-solve

use the Stage III-A candidate already in memory:

    addpath(genpath(pwd));

    R3 = main_stage3b_two_leg_physical_solve( ...
        'FrozenState', R0, ...
        'Candidate', Q3.candidate, ...
        'AllowSolve', true);

    disp(R3.summary);
    disp(R3.fitTable);
    disp(R3.EDI);
    disp(R3.prediction);
    disp(R3.solverInfo);
    disp(R3.gates);

The physical checkpoint is written to

    verification/crack_path/stage3b_two_leg_straight_physical_solved.mat

and the compact result to

    verification/crack_path/stage3b_two_leg_straight_physical_small_data.mat

## Stop condition

After a successful Stage III-B run, stop.

Do not generate the third leg automatically.

The next stage must first qualify a genuinely kinked polyline exterior and
new-tip core for the MTS-predicted theta_3 direction.
