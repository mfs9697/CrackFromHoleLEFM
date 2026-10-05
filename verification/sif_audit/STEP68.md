# Step 68: Level-1 SGS-PCG physical convergence experiment

## Purpose

Step68 is the first physical Level-0/Level-1 convergence comparison within the same qualified C03 structured mesh family.

It performs exactly one authorized Level-1 physical solve, using the free-DOF SPD + parameter-free SGS-PCG formulation qualified on Level 0 in Step67A.

## Fixed mesh family

Step68 reconstructs Level 1 deterministically from the committed archived Step62 candidate and requires exact agreement with Step65:

- T3 elements: 122,691;
- T6 nodes: 246,701;
- structured scale: 0.5;
- paired radius: 6 mm;
- maximum adjacent-size ratio: 1.79938154 within tolerance;
- complete T3/T6 reflection pairing;
- support fully inside paired region;
- exterior excluded from the primary support;
- six-triangle/seven-edge tip fan;
- native COD sampling: 74/108/86/67 in the four established windows.

No Level-0 physical solve is performed in Step68.

## Physical problem

Unchanged from the Level-0 audit:

- same crack geometry, a0=8 mm;
- same elastic material and plane-strain flag;
- unit remote-y traction;
- minimal rigid-body anchoring;
- same physical plate/hole geometry.

## Solver

The solver is locked to the Step67A-qualified formulation:

- assemble the unclamped symmetric stiffness with `stif_assem(...,[])`;
- restrict to free DOFs;
- symmetrize only machine-level antisymmetry after a strict symmetry gate;
- `symamd` ordering;
- parameter-free SGS preconditioner;
- PCG tolerance 1e-10;
- maximum 5000 iterations;
- no solver-parameter sweep;
- no direct `K\F` solve.

To reduce peak memory, Step68 releases the full stiffness matrix and the unpermuted free matrix immediately after constructing the permuted free system, before SGS is formed.

Step68 records MATLAB-reported available physical memory before assembly, before preconditioner construction, and after the SGS matrices are built when this information is available.

## Checkpoint-first rule

After a successful PCG solve, the Level-1 displacement field is immediately saved to:

`verification/step68_level1_sgs_physical_solved.mat`

before any COD or EDI postprocessing.

## COD convergence

Exactly the same four windows and polynomial degrees as Level 0 are used:

- [0.04,0.20] a0, degrees 1 and 2;
- [0.04,0.30] a0, degrees 1 and 2;
- [0.08,0.30] a0, degrees 1 and 2;
- [0.12,0.30] a0, degrees 1 and 2.

Step68 reports each Level-1 ratio beside its Level-0 counterpart and gives the signed relative Level-1 minus Level-0 change.

## EDI convergence

Exactly one physical interaction integral is evaluated:

- r_inner = 0.8 mm;
- r_outer = 5.2 mm = 0.65 a0;
- 16-point quadrature;
- FE-nodal-q weight;
- no radius sweep.

Step68 reports Level-0 and Level-1 KI, KII, and KII/KI together with their signed relative changes.

## Level-0 reference handling

If the local Step67A small-data result survives branch switching, its exact saved iterative values are used as the Level-0 reference.

If not, Step68 falls back to the recorded audit fingerprints:

- the eight high-precision COD ratios embedded in the Step67A audit;
- displayed/rounded EDI reference KI=0.43785, KII=4.6547e-5, KII/KI=1.0631e-4.

The reported result states which reference source was used.

## Physical-outcome policy

Step68 deliberately does **not** impose a pass/fail threshold on the unknown Level-1 physical convergence result.

The numerical gates check only:

- exact Level-1 mesh identity;
- PCG convergence;
- reported and true residuals;
- exact constraints;
- expected native sampling;
- finite COD and EDI values;
- one fixed EDI domain;
- no radius sweep;
- no solver tuning;
- no direct backslash;
- no Level-0 physical solve.

If those gates pass, the physical Level-0 to Level-1 change is reported and interpreted rather than forced to satisfy a preselected answer.

## Local run

On branch:

`audit/step68-level1-sgs-physical-convergence`

after explicit authorization run:

    addpath(genpath(pwd));

    R68 = main_step68_level1_sgs_physical_convergence( ...
        'AllowSolve', true);

    disp(R68.solverInfo);
    disp(R68.CODconvergence);
    disp(R68.EDIConvergence);
    disp(R68.CrossExtractor);
    disp(R68.gates);

If PCG fails with the fixed qualified settings, do not retune or rerun. Return the complete failure output.


## Completed local Step68 result — 2026-10-04

Step68 completed successfully with one explicitly authorized Level-1 physical SGS-PCG solve and all numerical validity gates passed.

### Level-1 linear solve

Problem size:
- T3 elements: 122,691;
- T6 nodes: 246,701;
- total DOFs: 493,402;
- free DOFs: 493,399.

Solver:
- method: free-DOF SPD SGS-PCG;
- PCG tolerance: 1e-10;
- maximum iterations: 5000;
- converged flag: 0;
- iterations: 3658;
- reported relative residual: 9.9450e-11;
- true relative residual: 9.9450e-11;
- constrained-displacement infinity norm: 0;
- stiffness symmetry error: 1.8420e-16;
- solve time: 160.8521 s.

Memory:
- available physical memory before assembly: 5.8527 GiB;
- available before preconditioner: 5.6755 GiB;
- available after preconditioner: 5.4804 GiB;
- SGS nnz(M1)=nnz(M2)=5,900,546;
- SGS approximate sparse storage: 0.1832 GiB;
- SGS setup time: 0.1012 s.

The Level-1 field was checkpointed before COD/EDI postprocessing.

### Level-0 to Level-1 COD convergence

All eight COD ratios increased slightly under the exact h/2 refinement:

- [0.04,0.20] a0, degree 1: +0.4173%;
- [0.04,0.20] a0, degree 2: +0.52153%;
- [0.04,0.30] a0, degree 1: +0.37750%;
- [0.04,0.30] a0, degree 2: +0.46301%;
- [0.08,0.30] a0, degree 1: +0.32487%;
- [0.08,0.30] a0, degree 2: +0.38465%;
- [0.12,0.30] a0, degree 1: +0.32057%;
- [0.12,0.30] a0, degree 2: +0.35110%.

Thus all COD estimators move consistently in the same direction, with Level-0 to Level-1 changes between approximately 0.32% and 0.52%.

Level-1 COD ratios span approximately 1.0449e-4 to 1.0623e-4, with mean 1.0565e-4 and median 1.0594e-4.

### Level-0 to Level-1 matched EDI convergence

Using the same fixed 0.8-5.2 mm interaction domain:

- KI: 0.43785 -> 0.43786, change +0.0032627%;
- KII: 4.6547e-5 -> 4.6643e-5, change +0.20618%;
- KII/KI: 1.0631e-4 -> 1.0652e-4, change +0.20111%.

The Level-0 EDI reference in this run came from the embedded audit fingerprint because the local Step67A small-data file was not available after branch switching. Its EDI values are therefore displayed/rounded, while the COD reference fingerprint is higher precision.

### Level-1 cross-extractor consistency

At Level 1:

- EDI ratio: 1.0652e-4;
- COD minimum: 1.0449e-4;
- COD maximum: 1.0623e-4;
- COD mean: 1.0565e-4;
- COD median: 1.0594e-4.

The EDI value lies just above the complete COD-fit interval. Relative to EDI:

- COD mean gap: -0.82243%;
- COD median gap: -0.54856%.

This is slightly tighter cross-extractor agreement than at Level 0, where the corresponding mean and median gaps were approximately -1.014% and -0.814%.

### Numerical gates

All Step68 gates passed:

- exact Level-1 mesh identity;
- PCG convergence;
- reported and true residuals;
- exact constraints;
- exact native sampling;
- finite COD and EDI;
- one matched EDI only;
- no radius sweep;
- no direct backslash;
- no solver tuning;
- no Level-0 physical solve.

### Interpretation

The controlled C03 refinement from Level 0 to Level 1 halves the structured scale while changing the matched EDI ratio KII/KI by only +0.201%.

All eight independent COD fits change by only approximately +0.32% to +0.52%, and all move in the same direction.

This provides the first direct mesh-family convergence evidence that the positive Mode-II component of order 1.06e-4 is a genuine feature of the asymmetric physical problem rather than a mesh-pairing or extraction artifact.

With only two physical mesh levels, no observed convergence order or Richardson-extrapolated continuum value is claimed. A third physical level would be required for a formal asymptotic-order estimate.
