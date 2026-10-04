# Step 67A: Level-0 SGS-PCG solver qualification

## Why Step67A exists

The first Step67 attempt successfully assembled the symmetric Level-0 free-DOF stiffness, with symmetry error 1.662e-16, but the fixed ICT incomplete-Cholesky preconditioner terminated with a nonpositive pivot before PCG began.

No physical iterative solve was performed in that failed attempt.

Step67A therefore changes **only the preconditioner**, not the mesh, physics, PCG tolerance, COD windows, EDI domain, or qualification tolerances.

## Preconditioner choice

Step67A uses a parameter-free symmetric Gauss-Seidel (SGS) preconditioner for the symamd-permuted free-DOF SPD matrix `A`.

Let

`G = tril(A) = D + L`.

Then the SGS preconditioner is

`M = G * D^(-1) * G'`.

The PCG split is supplied as:

- `M1 = G`;
- `M2 = D^(-1) * G'`.

Thus `M1*M2` is symmetric positive definite whenever the diagonal of `A` is positive.

Unlike incomplete Cholesky, SGS has no drop tolerance, no diagonal compensation, and no incomplete-factor pivot that can break down. There is therefore no post-failure tuning sweep.

## Guard and solve count

`AllowSolve=false` by default.

A run with explicit authorization:

    R67a = main_step67a_level0_sgs_solver_qualification( ...
        'AllowSolve', true);

performs at most one Level-0 PCG physical solve if no valid Step67A checkpoint exists.

The driver contains:

- one unclamped `stif_assem` call;
- one `pcg` call;
- no `ichol` call;
- no direct `K\F` call;
- one matched physical EDI after checkpointing;
- no Level-1 mesh or solve path.

## Fixed numerical settings

Unchanged from Step67:

- PCG relative tolerance: 1e-10;
- PCG maximum iterations: 5000;
- symamd ordering;
- true free-system residual gate: <=5e-10;
- exact homogeneous constraints.

## Physical qualification

After the PCG solution is safely checkpointed, Step67A uses:

- the same eight Level-0 COD fits as Step63;
- one EDI at r_inner=0.8 mm and r_outer=5.2 mm;
- 16-point FE-nodal-q interaction extraction;
- no radius sweep.

Acceptance gates remain unchanged:

- every COD ratio within 0.05% of the reproducible Step63R direct fingerprint;
- EDI KI, KII, and KII/KI each within 0.05% of the Step64 direct reference;
- if the old direct checkpoint still exists, complete displacement-vector agreement within 5e-8 relative norm.

## Meaning of a pass

A pass qualifies the SGS-PCG free-DOF formulation on Level 0 only.

It does not authorize or perform the Level-1 physical solve.

## Local run

On branch:

`audit/step67a-level0-sgs-solver-qualification`

after explicit authorization run:

    addpath(genpath(pwd));

    R67a = main_step67a_level0_sgs_solver_qualification( ...
        'AllowSolve', true);

    disp(R67a.solverInfo);
    disp(R67a.CODcomparison);
    disp(R67a.EDIcomparison);
    disp(R67a.gates);

Return the complete console output.
