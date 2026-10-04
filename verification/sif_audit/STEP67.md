# Step 67: Level-0 iterative-solver qualification

## Purpose

Step67 qualifies a memory-safer linear-algebra formulation on the already established Level-0 C03 physical problem before any Level-1 solve is considered.

It is **not** a mesh-convergence step. It deliberately returns to Level 0, where the direct-solve physical reference is known.

## Authorization and solve count

`AllowSolve=false` by default.

An explicitly authorized run with:

    R67 = main_step67_level0_iterative_solver_qualification( ...
        'AllowSolve', true);

performs exactly one physical PCG solve if no valid Step67 checkpoint exists.

The driver contains:

- zero calls to `solve_cracked_LEFM`;
- one call to `stif_assem`, with `fixvar=[]`;
- one call to `ichol`;
- one call to `pcg`;
- one physical `SIF_LEFM_interaction_EDI` call after checkpointing;
- no Level-1 builder or Level-1 solve path;
- no direct `K\F` solve.

## Exact physical problem

Step67 reconstructs deterministic Level-0 C03 directly from the committed archived Step62 baseline and requires:

- 32,980 T3 elements;
- 66,854 T6 nodes;
- maximum adjacent-size ratio 1.79678451 within tolerance;
- 6-mm paired region;
- complete T3/T6 reflection pairing;
- exact Level-0 native COD sampling 38/55/44/34 in the four established windows.

The physical problem is unchanged from Step63:

- same material;
- unit remote-y traction;
- minimal rigid-body anchoring;
- same crack geometry with a0=8 mm.

## Symmetric free-DOF formulation

The current production solver clamps only constrained rows and therefore stores a nonsymmetric matrix.

Step67 instead calls:

`stif_assem(mesh,mat,quad,[])`

to assemble the original unconstrained symmetric elasticity matrix `K`.

The three homogeneous constrained DOFs are then removed algebraically:

`Kff * uf = Ff`.

This is exactly equivalent to zero prescribed displacements but preserves the SPD structure needed by PCG.

Before solving, Step67 requires the unclamped stiffness symmetry error

`norm(K-K',1)/norm(K,1)`

to be at most `5e-13`.

Only after that gate, machine-level antisymmetry is removed by replacing `Kff` with `(Kff+Kff')/2`.

## Fixed iterative method

There is no solver-parameter sweep.

Predeclared settings:

- method: PCG;
- relative tolerance: 1e-10;
- maximum iterations: 5000;
- symmetric AMD permutation of `Kff`;
- incomplete Cholesky type: `ict`;
- drop tolerance: 1e-3;
- diagonal compensation: 1e-3;
- modified incomplete Cholesky: on.

If this single configuration fails, Step67 stops. It does not silently tune the preconditioner.

The true residual in the original free-DOF ordering must be <=5e-10 and the constrained displacements must remain <=1e-14.

## Checkpoint-first rule

After the one successful PCG solve, the physical field is immediately saved to:

`verification/step67_level0_iterative_physical_solved.mat`

before COD or EDI postprocessing.

This avoids repeating the Step63 lost-checkpoint problem.

## Physical reference fingerprints

Step67 then reloads only the saved iterative field and applies the same physical postprocessing used previously.

### COD

Exactly the Step63 windows and polynomial degrees are used:

- [0.04,0.20] a0, degrees 1 and 2;
- [0.04,0.30] a0, degrees 1 and 2;
- [0.08,0.30] a0, degrees 1 and 2;
- [0.12,0.30] a0, degrees 1 and 2.

The eight established direct-solve COD ratios from the reproducible Step63R field are embedded as the comparison fingerprint.

The maximum relative COD-ratio difference must be <=0.05%.

### EDI

Exactly one matched physical EDI is computed:

- r_inner = 0.8 mm;
- r_outer = 5.2 mm = 0.65 a0;
- 16-point quadrature;
- FE-nodal q weight;
- no radius sweep.

The established Level-0 direct reference is:

- KI = 0.43785;
- KII = 4.6547e-5;
- KII/KI = 1.0631e-4.

Each of KI, KII, and KII/KI must agree within 0.05%.

## Optional direct displacement comparison

If the recovered Step63 physical checkpoint happens to survive locally, Step67 additionally compares the complete iterative displacement vector against that direct-solve field.

This comparison is optional because generated MAT checkpoints have previously disappeared across branch switches.

When available, both relative 2-norm and infinity-norm displacement differences must be <=5e-8.

## Meaning of a pass

A passing Step67 establishes only that the free-DOF PCG/ichol formulation reproduces the established Level-0 direct physical solution to the predeclared tolerances.

It then sets:

`readyForLevel1IterativeProposal=true`.

This does **not** authorize or perform a Level-1 physical solve.

## Local run

On branch:

`audit/step67-level0-iterative-solver-qualification`

run:

    addpath(genpath(pwd));

    R67 = main_step67_level0_iterative_solver_qualification( ...
        'AllowSolve', true);

    disp(R67.solverInfo);
    disp(R67.CODcomparison);
    disp(R67.EDIcomparison);
    disp(R67.gates);

Return the complete console output.
