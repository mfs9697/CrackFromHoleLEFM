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

The new \`main_step47_local.m\` script reuses the local selected compact Step46 preflight data file \`verification/step46_mesh_preflight_face0125_tip05.mat\` and the saved *original* symmetric field \`verification/step45_symmetric_theta0_solved.mat\`. The defaults **cannot run a new solve**. If a recognized Step47 checkpoint exists, the script reuses it; otherwise, it requires explicit \`'AllowSolve',true\`.

Before solving, it:
1. Checks that the *selected* preflight really used \`FaceFactor=0.125\`, \`TipFactor=0.5\`, the correct polygonal hole and straight-crack geometry, all three quadratic COD windows, and the recorded tip/face gates.
2. Reconstructs the original centered right-half **baseline geometry and mesh** and compares sorted original T3 coordinates with the saved Step45 baseline.
3. Reuses the **same geometry-description struct** and reidentified exact sharp-pencil geometry IDs to generate the specified manual edge/vertex-refined candidate. Verifies the observed T3/T6 counts, both crack-face counts, pointwise native COD radial sampling, tip-edge median, tip-triangle counts and all three fitting-window populations against the saved preflight.
4. Rejects degenerate collapsed triangles. The saved preflight contains a compact radial profile and metrics, not the entire T3 connectivity; this reconstruction cannot claim bitwise equality to all unrecorded preflight interior triangles.
5. Calls \`solve_cracked_LEFM\` **exactly once**, only after explicit approval and all checks, using the unchanged centered-plate tension/constraint/material configuration. Verifies that the solver's T3 coordinates and connectivity match the accepted generated candidate.
6. Immediately saves **only** \`mesh,U,mat,crack,a0,meta\` to a unique local MAT checkpoint \`verification/step47_refined_symmetric_theta0_solved.mat\`. It does not overwrite the original Step45 checkpoint and does not save stiffness/stress arrays.
7. Runs the **same native COD postprocessor** used in Step45 and presents all preset linear/quadratic fit windows. **No EDI is run in this first refined-field step.** A rerun reuses the recognized checkpoint rather than solving again.

If the selected preflight file was saved elsewhere, supply its actual path using \`'PreflightFile',someAbsolutePath\`. The default above matches the investigator's latest successfully reported preflight path.

After pulling \`sif-asymmetric-mesh-audit\` through GitHub Desktop, the default is a zero-solve safety check:

\`\`\`matlab
addpath(genpath(pwd));
[P47,O47] = main_step47_local();
\`\`\`

When the investigator **explicitly approves spending one refined FEM solve**, use:

\`\`\`matlab
[P47,O47] = main_step47_local('AllowSolve',true);
disp(O47.fitTable);
disp(O47.CODgates);
\`\`\`

The second command may require MATLAB time for the actual single FEM solution. **Do not run it merely to inspect this report.** The code has been statically reviewed in GitHub, but numerical execution on the selected mesh has not yet been verified in MATLAB.

## How to compare results after that run

First compare the **signed** COD \`KII/KI\` obtained with the same linear/quadratic degrees across all three windows. A fitted value below \`1e-6\` at one window is not enough; examine consistency and residuals. Do not choose the fit that happens to be closest to zero.

Only after reviewing the COD results should we choose whether to perform a small 16-point FE-nodal EDI analysis. If we compare refined and original FEM fields, **keep the exact same absolute EDI integration radii**:

- \`r_outer/a0 = [0.50 0.65 0.80]\`, as in Step45.
- \`r_inner\` must match the original saved Step45 value \`0.0008675...\` m, NOT the smaller automatic \`2*hTip\` calculated from the newly refined mesh.

The original Step45 EDI result at \`0.65\` was \`KII/KI=-8.766949318e-6\`; this is an **observed coarse-control residual**, not a transferable error estimate for the refined symmetric solution or asymmetric physical crack. Start with **one** matched EDI domain after COD interpretation, not all three by default. The Step45 postprocessor accepts the Step47 checkpoint's backward-compatible case marker and uses a separate \`step47_refined_symmetric_field_leakage\` cache prefix.

**Stopping rule:** no further FEM solves or extensive EDI integration until these refined-COD and a deliberately chosen matched-EDI result are independently reviewed.
