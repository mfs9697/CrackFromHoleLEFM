# Step 47: one gated refined symmetric FEM control (not yet executed)

## Selection from completed Step 46 mesh-only experiments

The investigator selected the final **FaceFactor=0.125 / TipFactor=0.5** symmetric control mesh. This preserves exactly the original Npoly=240 polygonal centered-hole geometry and the prescribed a0=0.004 m horizontal crack, while requesting more crack-face and tip refinement. The saved Step45 baseline was independently regenerated with identical sorted T3 coordinate sets (maximum discrepancy 0 m). The selected trial produced:

| Measured property | Original Step45 mesh | Selected Step46 trial |
| --- | ---: | ---: |
| T3 triangles | 1721 | **2391** |
| T6 nodes | 3640 | **5054** |
| Measured median tip edge, m | 0.00043375122 | **0.00013535** (reported rounding) |
| Measured hTip/a0 | 0.108438 | **0.033837** |
| Tip-adjacent T3 above/below | 2 / 3 | **3 / 3** |
| Native upper/lower crack-face nodes | 12 / 12 | **72 / 72** |
| Maximum crack-face abscissa mismatch, m | 0 | **2.77556e-17** |
| COD native points in window 0.04–0.30 a0 | 4 | **19** |
| COD native points in window 0.08–0.30 a0 | 3 | **16** |
| COD native points in window 0.12–0.30 a0 | 2 | **13** |

**All three windows now exceed both the predeclared eight-native-point minimum for linear COD fits and twelve-native-point minimum for quadratic fits.** A 3/3 tip-triangle count is useful, but it does **not** establish reflection-paired node locations or symmetric numerical tractions. The mesher is also permitted to alter the triangulation away from the crack; this is a *same-geometry, selected-refinement-recipe* comparison, not an identical-exterior-connectivity comparison.

This is sufficient to **stop further Step 46 mesh-parameter trials**, not to declare physical small Mode II validated or numerical convergence established. The control is different from the much finer asymmetric-crack problem.

## Precisely one additional FEM solve (separate, explicit user action)

The new `main_step47_local.m` script reuses the local selected compact Step46 preflight data file `verification/step46_mesh_preflight_face0125_tip05.mat` and the saved *original* symmetric field `verification/step45_symmetric_theta0_solved.mat`. The defaults **cannot run a new solve**. If a recognized Step47 checkpoint exists, the script reuses it; otherwise, it requires explicit `'AllowSolve',true`.

Before solving, it:
1. Checks that the *selected* preflight really used `FaceFactor=0.125`, `TipFactor=0.5`, the correct polygonal hole and straight-crack geometry, all three quadratic COD windows, and the recorded tip/face gates.
2. Reconstructs the original centered right-half **baseline geometry and mesh** and compares sorted original T3 coordinates with the saved Step45 baseline.
3. Reuses the **same geometry-description struct** and reidentified exact sharp-pencil geometry IDs to generate the specified manual edge/vertex-refined candidate. Verifies the observed T3/T6 counts, both crack-face counts, pointwise native COD radial sampling, tip-edge median, tip-triangle counts and all three fitting-window populations against the saved preflight.
4. Rejects degenerate collapsed triangles. The saved preflight contains a compact radial profile and metrics, not the entire T3 connectivity; this reconstruction cannot claim bitwise equality to all unrecorded preflight interior triangles.
5. Calls `solve_cracked_LEFM` **exactly once**, only after explicit approval and all checks, using the unchanged centered-plate tension/constraint/material configuration. Verifies that the solver's T3 coordinates and connectivity match the accepted generated candidate.
6. Immediately saves **only** `mesh,U,mat,crack,a0,meta` to a unique local MAT checkpoint `verification/step47_refined_symmetric_theta0_solved.mat`. It does not overwrite the original Step45 checkpoint and does not save stiffness/stress arrays.
7. Runs the **same native COD postprocessor** used in Step45 and presents all preset linear/quadratic fit windows. **No EDI is run in this first refined-field step.** A rerun reuses the recognized checkpoint rather than solving again.

If the selected preflight file was saved elsewhere, supply its actual path using `'PreflightFile',someAbsolutePath`. The default above matches the investigator's latest successfully reported preflight path.

After pulling `sif-asymmetric-mesh-audit` through GitHub Desktop, the default is a zero-solve safety check:

```matlab
addpath(genpath(pwd));
[P47,O47] = main_step47_local();
```

When the investigator **explicitly approves spending one refined FEM solve**, use:

```matlab
[P47,O47] = main_step47_local('AllowSolve',true);
disp(O47.fitTable);
disp(O47.CODgates);
```

The second command may require MATLAB time for the actual single FEM solution. **Do not run it merely to inspect this report.** The code has been statically reviewed in GitHub, but numerical execution on the selected mesh has not yet been verified in MATLAB.

## How to compare results after that run

First compare the **signed** COD `KII/KI` obtained with the same linear/quadratic degrees across all three windows. A fitted value below `1e-6` at one window is not enough; examine consistency and residuals. Do not choose the fit that happens to be closest to zero.

Only after reviewing the COD results should we choose whether to perform a small 16-point FE-nodal EDI analysis. If we compare refined and original FEM fields, **keep the exact same absolute EDI integration radii**:

- `r_outer/a0 = [0.50 0.65 0.80]`, as in Step45.
- `r_inner` must match the original saved Step45 value `0.0008675...` m, NOT the smaller automatic `2*hTip` calculated from the newly refined mesh.

The original Step45 EDI result at `0.65` was `KII/KI=-8.766949318e-6`; this is an **observed coarse-control residual**, not a transferable error estimate for the refined symmetric solution or asymmetric physical crack. Start with **one** matched EDI domain after COD interpretation, not all three by default. The Step45 postprocessor accepts the Step47 checkpoint's backward-compatible case marker and uses a separate `step47_refined_symmetric_field_leakage` cache prefix.

**Stopping rule:** no further FEM solves or extensive EDI integration until these refined-COD and a deliberately chosen matched-EDI result are independently reviewed.


## Investigator's measured Step 47: refined FEM solved, COD evaluated

The investigator explicitly authorized and executed **exactly one** additional symmetric Stage-II FEM solution using the previously selected `FaceFactor=0.125`, `TipFactor=0.5` recipe. The script verified the achieved mesh before solving: **2,391 T3 triangles, 5,054 T6 nodes, hTip=0.000135346748 m (`hTip/a0=0.0338367`), 72 native nodes on each crack face with mismatch `2.77556e-17 m`, 3/3 tip-adjacent triangles, and 19/16/13 native samples in the three predefined COD windows**. The original centered-hole geometry and horizontal crack were unchanged. The refined solution was saved separately at `verification/step47_refined_symmetric_theta0_solved.mat`; no EDI was run.

All **six** original linear/quadratic COD fits could now be evaluated, but their signed `KII/KI` intercepts differed considerably:

| COD window, r/a0 | Native samples | Linear signed KII/KI | Quadratic signed KII/KI |
| ---: | ---: | ---: | ---: |
| 0.04–0.30 | 19 | -2.9487e-4 | -6.2706e-4 |
| 0.08–0.30 | 16 | -3.7903e-5 | +1.9368e-4 |
| 0.12–0.30 | 13 | -1.2246e-4 | -2.1741e-4 |

The `CODgates` result is `evaluated=true, passed=false, maxAbsRatio=6.27056963e-4` versus the **proposed** `1e-6` verification target. Unlike the coarse Step45 COD test, this is an actually **evaluated failure**. The quadratic fit in the middle window even changes the sign. The smallest Mode-II residuals from the narrowest window are **not** proof of an accurate zero-SIF intercept; fitting the apparent SIFs with ordinary unweighted `polyfit` does not enforce reflection symmetry or produce an uncertainty estimate.

Pointwise median native COD ratios by radial band are `-1.1755e-3` (0–0.04, only 2 samples), `-4.9558e-4` (0.04–0.08, 3 samples), `-7.0878e-5` (0.08–0.12, 3 samples), `-7.9394e-5` (0.12–0.20, 6 samples), `-5.0657e-5` (0.20–0.30, 7 samples). The near-tip bands and resulting extrapolations therefore remain sensitive to nonuniform numerical Mode-II signals despite sufficient sampling counts. The refined mesh's equal 3/3 tip-triangle *counts* do **not** establish exact mesh reflection symmetry.

**Interpretation:** COD now detects a strongly window/degree-dependent apparent tangential opening in the actually solved symmetric FEM field. This does not validate any extrapolated value as a physical `KII` (the continuum control must have `KII=0`), nor does it quantify error in a different, finer asymmetric crack geometry. Potential contributors include non-reflection-paired mesh details, the actual FEM displacement field near a singularity, higher-order fields and polynomial-extrapolation sensitivity.

## Next approved small diagnostic: one matched same-field EDI, zero FEM solves

The new `main_step48_refined_matched_edi.m` driver loads the **already saved Step45 baseline EDI result** and the **already saved Step47 refined FEM checkpoint**. It obtains the exact original **absolute** inner radius for the `r_outer/a0=0.65` domain directly from the saved baseline result, rather than recomputing it from the smaller refined `hTip`. It verifies the selected geometry, mesh provenance, endpoint/checkpoint identity and admissible annulus, then invokes the unchanged **16-point FE-nodal EDI** for **only this one domain**. It prints a matched original/refined table and caches its results separately.

After pulling `sif-asymmetric-mesh-audit` in GitHub Desktop, run from MATLAB:

```matlab
addpath(genpath(pwd));
R48 = main_step48_refined_matched_edi();
disp(R48.comparison);
```

The script generates **no new mesh or FEM solution**, and leaves the original and refined displacement checkpoints intact. It deliberately repeats the inexpensive native COD diagnostic while invoking EDI only at 0.65. Send the full MATLAB output before deciding whether to evaluate further EDI domains or to introduce a reflection-paired mesh test. Although both fields share the same physical problem and exact integration annulus, retriangulation outside the crack-tip neighborhood is possible, so the difference is an observed **mesh/field sensitivity**, not a cleanly isolated crack-tip-only error or a correction transferable to the refined asymmetric case.

**Current status:** the investigator completed the selected refined FEM solution, all six native COD fits and one 16-point FE-nodal EDI integral matched to the original symmetric control's exact absolute annulus. Mode-II EDI changed sign and increased in magnitude with this refinement. The next proposed step is the separate **no-solve/no-EDI [Step 49 reflection-parity audit](STEP49.md)** of both saved FEM meshes and crack-face displacements.

## Investigator's measured Step 48: identical absolute EDI annulus

The investigator ran the prepared Step 48 script on the **already saved** refined FEM solution. The 16-point FE-nodal EDI had precisely the original symmetric control's **absolute** integration annulus, `r_inner=0.00086750243077 m` and `r_outer=0.0026 m` (`r_outer/a0=0.65`). No new FEM solve, mesh or further EDI domains were computed.

| Matched physical control | T6 nodes | KI | KII | Signed KII/KI |
| --- | ---: | ---: | ---: | ---: |
| Original Step45 | 3640 | 0.3613086139 | -3.167574307e-6 | -8.766949318e-6 |
| Refined Step47 | 5054 | 0.3618471493 | +8.120223857e-6 | +2.244103310e-5 |

The reported **signed ratio change was `+3.1207982421e-5`**. Mode I changed by approximately 0.149%, whereas spurious Mode II reversed sign and increased in magnitude by approximately 2.56×. These two meshes **do not establish convergence** of the tiny Mode-II ratio. The difference is approximately 28.9% of the much finer *asymmetric* study's previously reported ratio `1.0795940665e-4`, but is strictly a **scale comparison, not an asymmetric-problem uncertainty bound**. Original and refined meshes share the physical configuration and integration annulus, but PDE remeshing altered connectivity outside the tip neighborhood as well.

The Step48 calculation also reran the saved refined-field COD postprocessor, reproducing the same six discordant linear/quadratic COD fits; no additional COD convergence can be inferred. The remaining causal questions concern actual FEM reflection symmetry and same-mesh integration/interpolation effects. Start with the **zero-solve, zero-EDI Step49 geometric reflection and displacement-parity diagnostic** documented separately at [STEP49.md](STEP49.md), rather than another full FEM solve or broad EDI sweep.
