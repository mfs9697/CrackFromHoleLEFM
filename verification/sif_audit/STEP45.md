# Step 45: finite-element-induced Mode-II leakage in the symmetric control

**Parent:** completed Step 44 analytical-field extraction replay on `sif-asymmetric-mesh-audit`.

## Scientific question

Step 44 demonstrated that both COD and 16-point FE-nodal interaction EDI recover **prescribed exact leading-order Williams displacements** accurately, including a tiny mixed Mode-II component. This does **not** establish that the actual finite-element displacement solution accurately resolves small Mode II.

Step 45 uses the existing **centered right-half plate with a circular hole**, a horizontal straight crack extending from the exact rightmost point of the hole at `theta = 0`, and remote symmetric vertical tension. The physical problem is symmetric under reflection across `y=0` and should have `KII=0` in the continuum. Numerical FEM discretization, asymmetric element topology, COD fitting and EDI quadrature may introduce residual nonzero estimates.

**Target:** `abs(KII/KI) <= 1e-6` for both extractors. This is a proposed numerical **verification threshold**, not a statistical error bar and not a prerequisite for claiming physical zero kink.

## Discovery: Step 18 output is not stored on the audited GitHub branch

`main_step18_centered_half_stage2_symmetry.m` retains all cases in `O18.Cases` **in MATLAB memory** but does not save a solved checkpoint to GitHub. The audited branch contains no Step 18 solved-field MAT artifact. Do not re-run the full multi-angle Step 18 experiment just to obtain a single symmetric control.

Also, Step 18 specified FE-nodal weight but left the quadrature at the extractor's backward-compatible default **7 points**. For Step 45, compute new 16-point results separately. Historical Step 18 EDI values are retained *only as metadata*, never relabeled as Step 45 16-point EDI.

## Files

- `main_step45_prepare_symmetric_checkpoint.m`: extracts one saved `theta=0` Step 18 solved case; if unavailable, permits **one** new Stage-II symmetric FEM solve only when `AllowSolve=true` is explicitly passed. The already-audited centered-hole geometry builder fixes the exact initiation point; it does not need a preliminary Stage-I solve for a prescribed zero-angle control. Saves only `mesh`, `U`, material, crack and compact provenance (never K/F/R/stress).
- `main_step45_symmetric_field_leakage.m`: validates checkpoint identity; runs native face COD, reports all pointwise distance bands, original fitting windows (if enough nodes), crack-face radial mismatch, and numbers of tip-adjacent T3 triangles above/below the axis. Optional `RunEDI=true` computes the **same actual saved FEM field** using explicit FE-nodal weight / 16-point Dunavant quadrature; integration is checkpointed after each requested radius.
- `test_step45_no_solve_guardrails.m`: cheap check that no absent solved field triggers an unexpected FEM solve and an absent postprocessing checkpoint fails immediately.

No Step 45 files change the production SIF selector, Step 44's test driver, or the existing numerical result files.

## First run (no FEM solve)

```matlab
addpath(genpath(pwd));
test_step45_no_solve_guardrails();

% If the old Step18 output O18 is still in memory or a local MAT file,
% reuse it. With no Npoly specified, selects its most resolved case.
P45 = main_step45_prepare_symmetric_checkpoint(O18);

% Zero-solve and zero-EDI initial stage:
O45 = main_step45_symmetric_field_leakage(P45);
disp(O45.tipTopology);
disp(O45.faceGridMismatch);
disp(O45.rawBands);
disp(O45.fitTable);
disp(O45.CODgates);
```

The original Step 18 output may occupy considerable memory. The checkpoint preparer discards the large original case structures before saving and never copies the original FEM stiffness matrix into its new checkpoint.

**If `O18` is NOT available:** do not rerun all of Step 18. After explicitly choosing to spend **one** FEM solve:

```matlab
P45 = main_step45_prepare_symmetric_checkpoint([], ...
    'AllowSolve',true,'Npoly',240);
O45 = main_step45_symmetric_field_leakage(P45);
```

Alternatively, provide a previously created Step-45 checkpoint directly to the postprocessor:

```matlab
O45 = main_step45_symmetric_field_leakage( ...
    'step45_symmetric_theta0_solved.mat');
```

A default call with neither an existing recognized checkpoint nor `O18` **aborts without computing**. An existing recognized checkpoint is never overwritten unless `Overwrite=true` is explicitly set.

## Matched 16-point interaction-EDI phase

**Only after examining the COD profile and tip topology**, reuse the SAME saved solved field:

```matlab
O45 = main_step45_symmetric_field_leakage(P45, ...
    'RunEDI',true,'ROuterOverA0',0.65);
disp(O45.EDI.table);
disp(O45.EDI.passed);
```

If the first EDI result warrants a domain test, request all three:

```matlab
O45 = main_step45_symmetric_field_leakage(P45, ...
    'RunEDI',true,'ROuterOverA0',[0.50 0.65 0.80]);
disp(O45.EDI.table);
```

The EDI progress cache is keyed by the saved checkpoint and each exact pair `(r_outer,r_inner)` together with the fixed 16-point FE-nodal method. Extending the radius list reuses previously computed matching domains; changing the checkpoint or numerical method is rejected. A changed inner radius for the same outer radius is tracked as a **separate** domain.

The default inner radius matches the Step-18 geometric rule `max(0.1*r_outer, 2*hTip)` to ensure the domain is admissible. Quadrature is explicitly 16-point and `StoreGPDiagnostics=false`. Results are stored after **each** completed radius, allowing interruption and resumption with identical settings.

If too few native nodes lie within a requested COD window, its fit is reported as missing; the driver **does not invent extra interpolated points or silently widen the window**. Inspect `O45.rawBands` and make a scientifically justified decision about additional resolution if necessary.

## Interpretation and stopping rules

- Exact physical symmetry gives `KII=0`, so report **absolute** `abs(KII/KI)`, not percentage error relative to zero.
- Compare COD and EDI on the **same saved actual FEM field**. EDI is not declared the physical truth merely because its domain variation is small.
- Report all tested COD fit windows and polynomial degrees; do not choose the result closest to zero.
- Upper/lower face radial-grid mismatch and asymmetric counts of tip-adjacent T3 elements may explain numerical residuals; neither proves the cause without additional tests.
- Passing `1e-6` on this symmetric geometry does **not** prove the asymmetric geometry's `KII/KI ~ 1.08e-4` is accurate to 1%; geometry-specific errors and higher-order Williams terms remain possible.

**Status:** code prepared; **not yet executed in MATLAB on an actual symmetric Step-18 FEM field**. Do not record the numerical leakage threshold as passed until the investigator provides both matching COD and 16-point EDI outputs.
