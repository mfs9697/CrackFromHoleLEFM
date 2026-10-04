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


## Completed local Step67A result — 2026-10-04

Step67A completed successfully with one explicitly authorized Level-0 physical SGS-PCG solve.

### Linear-algebra result

Problem size:
- T3 elements: 32,980;
- T6 nodes: 66,854;
- total DOFs: 133,708;
- free DOFs: 133,705.

Symmetric free-DOF formulation:
- unclamped stiffness symmetry error: 1.6620e-16;
- SGS preconditioner storage estimate: 0.0494 GiB;
- SGS construction time: 0.0698 s;
- nnz(M1)=1,590,857;
- nnz(M2)=1,590,857.

PCG:
- tolerance: 1e-10;
- maximum iterations: 5000;
- converged flag: 0;
- iterations: 1961;
- reported relative residual: 9.7490e-11;
- true free-system relative residual: 9.7491e-11;
- constrained-displacement infinity norm: 0;
- solve time: 25.3746 s.

The iterative physical field was checkpointed before COD/EDI postprocessing.

### COD comparison against the direct Level-0 reference

All eight COD-ratio fits reproduced the direct Step63R fingerprint extremely closely.

The relative differences ranged from approximately 9.42e-05% to 1.48e-04%.

The largest absolute ratio difference was approximately 1.57e-10.

All COD qualification gates passed.

### EDI comparison against the direct Level-0 reference

Using the single fixed 0.8-5.2 mm interaction domain:

- KI direct reference: 0.43785;
- KI iterative result: 0.43785;
- displayed relative gap: -0.0010305%.

- KII direct reference: 4.6547e-05;
- KII iterative result: 4.6547e-05;
- displayed relative gap: +0.00077941%.

- KII/KI direct reference: 1.0631e-04;
- KII/KI iterative result: 1.0631e-04;
- displayed relative gap: +1.926e-05%.

These direct-reference values were recorded previously at limited displayed precision, so the tiny displayed KI/KII gaps include reference-rounding uncertainty. The ratio agreement is effectively exact on the present reporting scale.

### Qualification gates

All Step67A gates passed:

- PCG converged;
- reported and true residual gates passed;
- exact homogeneous constraints;
- native sampling unchanged;
- all COD ratios matched within 0.05%;
- KI, KII, and KII/KI EDI matched within 0.05%;
- no direct backslash was used;
- no Level-1 physical solve was performed.

The old Step63 direct displacement checkpoint was not available locally, so the optional full-vector displacement comparison could not be performed. This did not affect the required physical-reference gates.

### Interpretation

Step67A establishes that the memory-safer free-DOF SPD + SGS-PCG formulation reproduces the established Level-0 direct physical solution to far tighter accuracy than the predeclared 0.05% tolerance.

The SGS preconditioner is compact and parameter-free, avoiding the incomplete-Cholesky breakdown observed in Step67.

Therefore the SGS-PCG formulation is qualified for a future Level-1 physical convergence solve, subject to separate explicit authorization.
