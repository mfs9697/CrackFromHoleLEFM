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

## GitHub Desktop + one-command MATLAB workflow (recommended)

The Step-45 code is merged into the **`sif-asymmetric-mesh-audit`** branch, not the production `main` branch. In GitHub Desktop, select the audit branch and **Fetch origin → Pull origin**. The `main_step45_local.m` script locates the repository itself, so the current MATLAB folder does not matter as long as the script is on the MATLAB path.

From a MATLAB session with the repository on its path:

```matlab
main_step45_local
```

The script first runs the cheap no-solve guard test. It then reuses an existing identified checkpoint under `verification/step45_symmetric_theta0_solved.mat`, or the previously solved `O18` variable if still in the MATLAB workspace. It runs **only native COD**, prints the tip-side triangle counts, crack-face grid mismatch, pointwise COD bands, original-window linear/quadratic extrapolations and a proposed leakage gate. Send this complete output to the investigator before starting EDI.

If neither a saved checkpoint nor `O18` exists, it prints a message and **does not run a FEM solve**. If an older Step-18 MAT file is available locally, load its `O18` variable, then rerun the same script. If no solved results can be recovered and the investigator explicitly authorizes **one** new symmetric control solve, run:

```matlab
STEP45_ALLOW_NEW_SOLVE = true;
main_step45_local
```

This creates and immediately saves **exactly one** Npoly=240 horizontal zero-angle Stage-II field and then produces COD results. No Stage-I solve or full Step-18 angle sweep is performed. The existing recognized checkpoint is always reused on subsequent runs. Keep the large `verification/step45_symmetric_theta0_solved.mat` local; only report numerical results or share a compact result file when needed.

**Do not run EDI until we have inspected the COD result.** The EDI phase remains a separate, opt-in call shown below.

## Alternative manual first run (no FEM solve)

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

## Measured results: Step 45 phase 1 and first 16-point EDI (investigator-supplied MATLAB log)

The investigator explicitly authorized **one** centered zero-angle Stage-II FEM solve using `STEP45_ALLOW_NEW_SOLVE=true`; the code immediately checkpointed it locally at `verification/step45_symmetric_theta0_solved.mat`. The original geometry-ID identification fell back to the temporary background-mesh route and then succeeded; this is an identification fallback, **not** a failure of the FEM solve.

- Symmetric control: `Npoly=240`, `a0=0.004 m`, **3,640 T6 nodes**, **1,721 T3 triangles**.
- COD native crack-face samples: **12 upper / 12 lower**, identical radial abscissae (`gridMismatch=0`).
- Original Step-39-style COD windows `[0.04,0.30]`, `[0.08,0.30]`, and `[0.12,0.30]` contain only **4, 3, and 2** native samples, respectively. All requested linear and quadratic fits correctly **skipped**. `CODgates.evaluated=false`; do not interpret `passed=false` here as an accuracy failure.
- Characteristic tip-edge median `hTip=0.00043375122 m`, `hTip/a0=0.108438`, tip-adjacent T3 elements **2 above / 3 below** the crack line.
- Raw pointwise COD ratios (one native sample per listed band): `+2.6801e-3` for `0.04–0.08`, `+3.102e-4` for `0.08–0.12`, `−1.6581e-6` for `0.12–0.20`, and `−5.3765e-5` for `0.20–0.30`. These **are not extrapolated crack-tip SIF ratios**.
- First opt-in EDI on **the identical saved actual FEM field**: FE-nodal weight, 16-point quadrature, `r_outer/a0=0.65`, `r_inner=0.0008675 m` (= `2*hTip`, dominating `0.1*r_outer`). Recovered **`KI=3.613086139e-1`**, **`KII=-3.167574307e-6`**, **`KII/KI=-8.766949318e-6`**.
- Magnitude of this symmetric-control residual corresponds to **8.1206%** of the **different refined asymmetric calculation's** ratio `1.0795940665e-4`. This is a **scale comparison only**, not a transferable error estimate: the symmetric control has a different crack length, material geometry, load distribution and vastly coarser mesh.
- The EDI result **exceeds** the proposed `1e-6` control threshold by ~8.77× for this one domain. The COD gate is **not evaluated**; neither an EDI single-domain residual nor one-point COD bands identify the causal error source. The Step-44 exact-field replay showed only much smaller extractor-only EDI leakage on a different, finer mesh; it does not bound actual FEM errors here.

**Next incremental experiment, no new FEM solve:** extend only the 16-point EDI outer-radius list on the **saved symmetric checkpoint** to `[0.50,0.65,0.80]`. The Step-45 per-domain progress cache should reuse the completed `0.65` domain automatically. On this mesh `2*hTip≈0.00086750244 m` exceeds `0.1*r_outer` at all three radii, so the three tests have the **same inner radius**; this isolates outer-radius dependence more cleanly. Interpret any domain stability as **consistency**, not proof of absence of a common discretization bias.

```matlab
O45 = main_step45_symmetric_field_leakage(P45, ...
    'RunEDI',true, 'ROuterOverA0',[0.50 0.65 0.80], ...
    'SavePrefix',fullfile(fileparts(P45.checkpointPath), ...
                          'step45_symmetric_field_leakage'));
disp(O45.EDI.table);
disp(O45.EDI.domainRatioSpread);
```

Stop and review this table before choosing any mesh refinement. If the EDI residual changes materially with outer radius, investigate quadrature-domain sensitivity on the **same** control field. If it is stable but remains above the target, a controlled symmetric FEM mesh-quality experiment is justified, not an automatic claim that the asymmetric EDI suffers the same relative error.

## Measured Step 45 three-domain EDI (investigator MATLAB, 2026-10-02)

On the **same** saved actual centered-symmetry control (Npoly=240, 3,640 T6 nodes, a0=0.004 m; no additional FEM solution), the investigator evaluated the original FE-nodal **16-point EDI** at three outer radii. The inner radius was fixed at `0.0008675 m` across all domains, set by `2*hTip` (not by `0.1*r_outer`).

| r_outer/a0 | actual KI | actual KII | signed KII/KI | `abs(ratio)/1.0795940665e-4` |
| ---: | ---: | ---: | ---: | ---: |
| 0.50 | 0.3613148899 | -5.649036290e-6 | -1.563466231e-5 | 0.14482 |
| 0.65 | 0.3613086139 | -3.167574307e-6 | -8.766949318e-6 | 0.081206 |
| 0.80 | 0.3613093612 | -2.372589674e-6 | -6.566643239e-6 | 0.060825 |

The 0.65 case **was reused from the progress cache**; only two additional EDI integrations ran. The signed **ratio spread was `9.068019071e-6`**, over nine times the proposed `1e-6` symmetry tolerance. Across these domains, `KI` varies by only ~0.00174% (range divided by ~0.36131). Meanwhile the spurious `KII` changes substantially and retains a negative sign.

**Interpretation:** On this particular coarse control mesh, Mode II is **not domain independent** at the scale required for the tiny asymmetric signal. The decreasing magnitude with increasing outer radius does **not** justify extrapolation to `KII=0`, nor does it prove which component of the numerical calculation is responsible. COD is still **not evaluable** with the original fitting windows (at most four native points). Neither the sign of this spurious symmetric-control residual nor its percentage of the different asymmetric signal is physically transferable.

**Next inexpensive discriminating experiment, Step 45b:** Replay a prescribed **exact pure Mode-I** leading Williams displacement field on **the same coarse Step-45 T6 mesh** and perform **only one** 16-point FE-nodal EDI calculation initially at `r_outer/a0=0.65`, with the **exact same** recorded `r_inner=0.0008675 m`. Verify the exact COD jump on this mesh first. Compare this exact-field EDI pure-I→II leakage directly against the actual-FEM signed residual at the matched annulus. This test introduces **no new FEM solution**, but exact-field replay is not an independent global equilibrium test; generator and EDI auxiliary fields share analytical conventions. If exact-field leakage is negligible relative to actual-FEM leakage, the *computed FEM field* is implicated; if comparable, investigate same-mesh EDI interpolation/quadrature before another solve.

Pull the updated `sif-asymmetric-mesh-audit` branch using GitHub Desktop and run from the repository root with the previously returned `O45` still in MATLAB:

```matlab
addpath(genpath(pwd));
R45b = main_step45_coarse_exact_pureI_replay(O45);
disp(R45b.table);
```

If `O45` is not in the workspace, load it from the **existing compact file** (this does not rerun COD or EDI):

```matlab
load(fullfile('verification', ...
    'step45_symmetric_field_leakage_small_data.mat'),'O45');
R45b = main_step45_coarse_exact_pureI_replay(O45);
disp(R45b.table);
```

The new driver uses the solved checkpoint path stored in `O45`, verifies the exact COD field against its known coefficients, executes the **single** matched exact EDI, and caches its result. If the result warrants exploring domain dependence, call the same driver with `'ROuterOverA0',[0.50 0.65 0.80]`: the 0.65 case will be reused, provided the checkpoint and integration options have not changed. Wait for the 0.65 result before running further domains.

## Interpretation and stopping rules

- Exact physical symmetry gives `KII=0`, so report **absolute** `abs(KII/KI)`, not percentage error relative to zero.
- Compare COD and EDI on the **same saved actual FEM field**. EDI is not declared the physical truth merely because its domain variation is small.
- Report all tested COD fit windows and polynomial degrees; do not choose the result closest to zero.
- Upper/lower face radial-grid mismatch and asymmetric counts of tip-adjacent T3 elements may explain numerical residuals; neither proves the cause without additional tests.
- Passing `1e-6` on this symmetric geometry does **not** prove the asymmetric geometry's `KII/KI ~ 1.08e-4` is accurate to 1%; geometry-specific errors and higher-order Williams terms remain possible.

**Status:** the symmetric FEM solve, COD-only phase and all three matched 16-point EDI domains have now been executed by the investigator. The EDI numerical symmetry target is not met, and Mode II displays considerable domain dependence while Mode I remains stable. COD extrapolation remains unevaluable on this coarse mesh. Step 45b's **same-mesh exact pure-I EDI replay** is prepared but not yet run in MATLAB; wait for that diagnostic before another FEM solve.
