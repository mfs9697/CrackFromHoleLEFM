# Step 53: true common-grid reflection-parity comparison (no new solve)

## Motivation

The investigator completed Step52 successfully: all mirrored spatial-point queries were located on both already solved symmetric T6 FEM meshes; affine reflected-field self-check residuals were ~1.4e-17 m, and direct analytical pure Mode I reflection self-checks returned zero.

At *each mesh's own upper-node sampling positions*, actual FEM parity residuals became smaller after refinement:

| Normalized residual | Original Step45 | Refined Step47 |
| --- | ---: | ---: |
| COD tip-disk, even ux | 4.7428e-4 | 1.2113e-4 |
| COD tip-disk, gauge-removed odd uy | 5.4888e-4 | 1.4879e-4 |
| Matched EDI-annulus, even ux | 1.3587e-4 | 7.2582e-5 |
| Matched EDI-annulus, gauge-removed odd uy | 1.0913e-4 | 6.2692e-5 |

The separate exact pure-I **nodal** field also showed generally decreasing parity RMS when processed through the two different T6 meshes. Yet the actual EDI `KII/KI` residual became *larger* in magnitude and reversed sign after refinement. We cannot use parity RMS as a proxy for EDI sign or numerical error.

**Step52 limitation:** its original/refined COD regions used **30 versus 176 upper FEM nodes**, and its EDI regions used **70 versus 342 nodes**, at *different physical positions*. Each FEM field also used its own native face-opening denominator. Therefore Step52 is not an identical-sampling cross-mesh convergence test.

## Precisely one new diagnostic, same saved FEM solutions

`main_step53_common_grid_parity.m` fixes these issues **without any FEM solve, remeshing, EDI integral or SIF fit**. The driver reads the original and refined FEM checkpoints referenced by the previously saved Step48 comparison. It defines two **predeclared identical physical sampling grids**, each with 7 angles `[25 45 65 90 115 135 155]` degrees in the upper half-plane and their exact reflected lower locations:

- COD-tip region: 6 radii `r/a0=[0.08 0.12 0.16 0.20 0.24 0.28]` (42 planned mirror pairs).
- Step48-matched EDI annulus: 7 radii `r/a0=[0.24 0.30 0.36 0.42 0.48 0.54 0.60]` (49 planned mirror pairs). The original absolute inner and outer radii are verified from Step48 before evaluating.

Each mesh's saved T3 geometry is used **only** for `triangulation/pointLocation`. Evaluations use the saved matching T6 connectivity and quadratic barycentric shape functions. The script retains **only** those physical sample pairs successfully located in BOTH meshes; it refuses a cross-mesh comparison if fewer than 80% of planned pairs are shared. Thus reported two-mesh RMS values are evaluated at **exactly the same physical points**, even if generated mesh node sets differ.

For each mesh and region, the script reports:

1. Actual FEM RMS even `u_x` and gauge-removed odd `u_y` parity residuals in **physical displacement units**, and normalized by **one shared normal-opening denominator**. This denominator uses interpolated actual normal openings from BOTH solved meshes at the **same four predetermined native-face physical radii** `r/a0=[0.14 0.18 0.22 0.26]`.
2. Synthetic prescribed **pure-I nodal** Williams-field parity RMS on the same points, normalized by **one shared exact pure-I opening** at those same physical reference radii. These synthetic values are diagnostics of T6 interpolation for that field, not corrections to the solved FEM error.
3. Exact affine even/odd reflected-field self-check. Since an affine displacement is represented exactly by T6 interpolation, parity must hold to tight absolute tolerance.
4. Directly evaluated pure-I Williams field at the same upper and lower physical points; exact physical reflection must hold independently of interpolation.

**Interpretation:** this is a controlled same-physical-sampling *displacement-parity comparison*. Even if parity decreases with refinement under this test, the FEM field and EDI extraction are distinct and we still cannot attribute the small actual EDI sign reversal to any single cause, establish physical `KII` accuracy, or transfer a numerical correction to the separate asymmetric plate. The original and refined meshes also differ in connectivity away from the tip.

## One local run

After pulling the `sif-asymmetric-mesh-audit` branch from GitHub Desktop, execute:

```matlab
addpath(genpath(pwd));

R53 = main_step53_common_grid_parity();

disp(R53.summary);
```

The code only reads existing local checkpoints plus `verification/step48_refined_matched_edi_comparison_small_data.mat`. It creates one compact diagnostic `verification/step53_common_grid_parity_small_data.mat`, containing summary statistics and the fixed small common query-point grids. The two saved full FEM checkpoints remain unchanged.

**Stop after this one run** and inspect query coverage, the affine and analytical self-checks, raw physical displacement errors, and shared-denominator actual/synthetic parity statistics. Depending on these results, consider whether designing a fully reflection-paired mesh is worth an additional separately authorized FEM solution. Do not treat the static source audit as a substitute for MATLAB execution.


## Investigator's measured Step 53: identical-point, common-normalization parity

The investigator completed the script on **both previously solved symmetric control checkpoints**. All planned upper/lower mirror queries were found on both meshes: **42/42 COD-region pairs** and **49/49 original-EDI-annulus pairs**. This time, unlike Step52, each reported coarse/refined pair is evaluated at **identical predetermined physical coordinates**, and the two actual FEM fields use **the same reference opening**. The affine T6 interpolation self-check passed with maximum parity errors ~2.14–2.82e-17 m, and direct analytical exact pure-I parity errors were zero.

The normal crack openings on the predetermined common physical reference radii were `1.341309942e-7 m` (original) and `1.371813053e-7 m` (refined); the **single shared** normal-opening denominator was `1.356561497e-7 m`. The shared exact pure-I synthetic opening denominator was `3.886608669e-7 m`.

| Common physical region and parity measure | Original Step45 | Refined Step47 |
| --- | ---: | ---: |
| COD disk: actual RMS even ux / shared actual opening | 3.8058e-4 | 9.8011e-5 |
| COD disk: actual RMS gauge-corrected odd uy / shared actual opening | 4.7716e-4 | 2.0162e-4 |
| Matched EDI annulus: actual RMS even ux / shared actual opening | 1.6691e-4 | 9.5775e-5 |
| Matched EDI annulus: actual RMS gauge-corrected odd uy / shared actual opening | 1.2908e-4 | 1.1304e-4 |
| COD disk: nodal exact pure-I RMS even ux / exact opening | 3.1129e-4 | 9.833e-5 |
| COD disk: nodal exact pure-I RMS gauge-corrected odd uy / exact opening | 3.933e-4 | 1.9506e-4 |
| Matched EDI annulus: nodal exact pure-I RMS even ux / exact opening | 1.1564e-4 | 8.0147e-5 |
| Matched EDI annulus: nodal exact pure-I RMS gauge-corrected odd uy / exact opening | 1.2846e-4 | 1.1539e-4 |

At the **same spatial locations and with shared normalization**, the actual FEM parity residuals decline between the saved meshes: approximately 74%/58% for even/odd COD-region components, and 43%/12% for even/odd EDI-annulus components. These are the actual displacement RMS changes under this experiment, **not** changes in physical SIF or EDI recovery. The prescribed exact-nodal pure-I interpolation parity similarly declines under this test, reinforcing the hypothesis that spatial representation contributes to the observed numerical asymmetries.

**Crucial contrast:** the original/refined actual FEM 16-point EDI ratio at the exact matched annulus changed from `-8.766949318e-6` to `+2.244103310e-5` (larger absolute value despite smaller parity RMS). Meanwhile prescribed exact pure-I **nodal** EDI changed `-3.88687e-5 → +3.37937e-5`, whereas **direct Gauss analytical** EDI yielded approximately `-1.147788e-7` and `+6.95969e-10`. Thus decreasing spatial displacement-parity RMS **does not imply** converged EDI Mode II or quantitatively explain the actual-FEM sign reversal.

### Next boundary-only gate before designing a fully mirror-paired mesh

The two solved meshes are not reflection-paired off the crack faces (Step49). Before building a new mesh, verify that the **prescribed polygonal geometry** reconstructed using the exact Step45/47 Npoly=240, a0=0.004 m configuration is itself reflection-symmetric at its **complete vertex and boundary-edge sets**. If it is not, constructing a perfectly mirrored volume mesh could silently change the benchmark rather than merely improve its discretization.

The separate [Step 54 geometry symmetry preflight](STEP54.md) checks that prerequisite **without running any mesh generator, FEM solve or EDI**. No additional FEM solve has been authorized.
