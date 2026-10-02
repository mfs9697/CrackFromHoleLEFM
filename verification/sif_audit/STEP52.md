# Step 52: reflected-point T6 displacement-parity diagnostic (no FEM solve)

## Scientific motivation

**Completed Step51:** on the selected 5,054-node refined T6 symmetric control mesh at the same Step48 absolute annulus (r_inner=0.00086750243077 m, r_outer=0.0026 m), prescribed pure-I **exact nodal** Williams field gave signed EDI q≈+3.37937e-5; evaluating the **same** prescribed analytical field directly at Gauss points gave signed q≈+6.95969e-10 (approximately **48,600-fold** reduction). On the coarse control mesh, the analogous nodal/Gauss values were −3.88687e-5/−1.14779e-7. Exact-field interpolation contaminates EDI significantly and with mesh-dependent sign. This **does not** quantitatively explain the actual-FEM sign reversal (−8.76695e-6 → +2.24410e-5) or validate the much-finer asymmetric physical q≈1.08e-4.

**Completed Step49:** neither saved symmetric mesh is reflection-paired off the crack faces. For the selected refined mesh, there were 0 exact reflected partners among 494 upper off-face nodes in the tip disk, and the nominally equal 3/3 tip triangle fan was **not** a complete mirror assignment. Yet both meshes' *crack-face* radial grids match to coordinate roundoff, and the refined actual FEM field has a persistent raw negative tangential opening. We cannot evaluate actual off-face displacement parity using **coincident nodes**, but the existing T6 finite-element field can be evaluated at any **off-face point inside an element**.

## One deliberately limited next experiment

The new driver `main_step52_interpolated_reflection.m` evaluates **both saved FEM solutions** at a collection of upper off-face T6 nodes and at their geometrically reflected points on the lower side. It uses MATLAB `triangulation/pointLocation` on the saved T3 geometry to find the triangle containing each reflected point and then evaluates the **correct six T6 interpolation shape functions** in the corresponding saved T6 element. There are **no new FEM equations or remeshing**.

It performs the same spatial sampling with **four deliberately distinct inputs/checks**:

1. **Actual computed FEM displacement** on each mesh. Measure the RMS even-parity residual `u_x(x,+y)-u_x(x,-y)` and the gauge-removed RMS odd-parity residual `u_y(x,+y)+u_y(x,-y)`. Normalize by measured RMS native mode-I crack opening in the narrowest original window `0.12–0.30 a0`. Do NOT interpret these RMS values as fitted `KII/KI`.
2. **Prescribed pure-I leading Williams displacements sampled at the exact same T6 nodes**. Repeat the identical reflected-point T6 interpolation on each mesh. Report residuals with **their own** pure-I reference opening. This controls *interpolation* of a known symmetric singular field, not the actual FEM displacement error.
3. **An affine, reflection-compatible synthetic field** `u_x=x-x_tip,\ u_y=y-y_tip`. T6 polynomial interpolation must reproduce this exactly at reflected spatial points. The driver aborts if the maximum parity discrepancy exceeds its tight numerical tolerance, catching element connectivity/shape-function errors.
4. **Direct analytical pure-I Williams values evaluated at the reflected point**. Confirm that the analytic upper sample and directly evaluated lower sample satisfy even/odd reflection to numerical tolerance. This checks angular branch conventions independently of T6 interpolation.

The driver checks the saved Step48 annulus and the recognized, unchanged physical Step45 and Step47 checkpoints. It uses exactly two regions: the **fixed physical COD tip disk** `r/a0 in [0.04,0.30]` and the **previously matched Step48 EDI annulus** `r in [0.00086750243077,0.0026] m`. It excludes points with `y-y_tip<=0.00008 m` so interpolation does not sample locations arbitrarily close to the collapsed crack faces. Unlike Step49, unmatched nodal layouts are **not** grounds for refusing physical-point interpolation: `pointLocation` finds an arbitrary lower element where the mirrored spatial point lies. It reports the count and fraction of actual located mirror queries and returns `NaN` for parity if fewer than two succeed.

**Interpretation limits:** T6-interpolated reflected displacements measure the already-solved FE fields at the same physical locations and resolve the Step49 *missing exact nodal partner* obstacle. An interpolation residual for the nodally prescribed Williams field is a *calibration for that synthetic field*, not an additive correction for actual FEM. Results may vary with the sampled spatial region and near-slit exclusion, and a reflection-symmetric continuum control does not prove physical `KII` recovery on an asymmetric geometry. Do not claim that a smaller parity RMS is equivalent to smaller EDI Mode-II leakage.

## One local MATLAB run

After updating branch `sif-asymmetric-mesh-audit` through GitHub Desktop:

```matlab
addpath(genpath(pwd));
R52 = main_step52_interpolated_reflection();
disp(R52.summary);
```

The driver reads the existing `verification/step48_refined_matched_edi_comparison_small_data.mat` plus the saved original/refined FEM checkpoints referenced by that file. It uses the existing checked `native_COD_audit.m` and `exact_williams_displacement_audit.m`. It saves only compact summaries in `verification/step52_interpolated_reflection_small_data.mat`. **No new FEM solve, geometry creation, mesh generation or EDI integration** is performed.

Send the complete MATLAB output, especially the affine and directly evaluated exact-field self-checks, mirror-query coverage and normalized *actual versus prescribed-field* parity values. Only after interpreting this genuinely same-physical-point test should we commit to the much larger development and separately approved FEM solution on an exactly reflection-paired mesh.


## Measured Step 52: all reflected queries located; reduced mesh-node-sampled parity on refinement

The investigator completed the Step52 two-checkpoint MATLAB run **without any new mesh, FEM solve, EDI or polynomial fit**. Both affine and directly analytical controls passed, and **100%** of reflected lower spatial queries were found inside the saved T6 meshes. This resolves Step49's inability to compare noncoincident opposite-side nodal positions, but each mesh was still sampled at **its own upper T6 node coordinates**.

| Field / region / parity | Original Step45 | Refined Step47 |
| --- | ---: | ---: |
| Number of off-face upper query points, COD tip disk | 30 (30 located) | 176 (176 located) |
| Number of off-face upper query points, matched EDI annulus | 70 (70 located) | 342 (342 located) |
| Max affine reflected-point T6 parity residual [m] | ~1.4e-17 | ~1.4e-17 |
| Direct analytic pure-I reflected-point parity residual | 0 | 0 |
| Actual FEM COD disk: RMS even ux / native opening | 4.7428e-4 | 1.2113e-4 |
| Actual FEM COD disk: RMS gauge-corrected odd uy / native opening | 5.4888e-4 | 1.4879e-4 |
| Actual FEM EDI annulus: RMS even ux / native opening | 1.3587e-4 | 7.2582e-5 |
| Actual FEM EDI annulus: RMS gauge-corrected odd uy / native opening | 1.0913e-4 | 6.2692e-5 |
| Nodal exact pure-I COD disk: RMS even ux / pure-I opening | 6.5483e-4 | 1.1051e-4 |
| Nodal exact pure-I COD disk: RMS gauge-corrected odd uy / pure-I opening | 5.8234e-4 | 1.3444e-4 |
| Nodal exact pure-I EDI annulus: RMS even ux / pure-I opening | 7.5261e-5 | 4.4385e-5 |
| Nodal exact pure-I EDI annulus: RMS gauge-corrected odd uy / pure-I opening | 6.8307e-5 | 5.347e-5 |

The actual FEM displacement-parity residual decreases in both sampled regions as the mesh is refined, **under this particular mesh-dependent sampling and each field's own opening normalization**. The nodally prescribed pure-I synthetic field exhibits the same broad decrease in parity RMS. Thus Step52 makes it plausible that near-tip spatial interpolation and discrete field asymmetry are changing together with mesh refinement. However, **the actual EDI q changed in the opposite direction in magnitude** and reversed sign (Step48, `-8.76695e-6 → +2.24410e-5`). A smaller parity RMS is **not** a validated predictor of smaller EDI Mode-II contamination.

**Important sampling limitation:** 30 vs 176 COD-disk locations and 70 vs 342 EDI-annulus locations are **different spatial point sets**. The prior Step52 reference opening denominator is also computed separately from each mesh's actual native face nodes. These results support a within-mesh physical-point parity test and an exploratory cross-mesh trend, **not** a formal at-identical-points convergence measurement or direct decomposition of actual FEM error.

### Next smallest zero-solve control: two meshes at one shared physical sampling grid

[Step 53](STEP53.md) replaces mesh-dependent upper-node locations with the **same predeclared polar-grid physical points on both meshes**. Every evaluation uses each saved mesh's own T6 interpolation at the exact same upper and reflected lower coordinates; parity RMS comparisons retain only spatial points successfully located **in both meshes**. It also uses the **same physical crack-face reference radii** on both solved fields and one shared FEM normal-opening normalization, while preserving separate pure-I analytical normalization. This reduces two important sampling confounders before we decide whether to construct and separately authorize a new *reflection-paired* FEM mesh.

No further FEM solve, mesh generation, EDI integration, or wide radius sweep is proposed at this stage.
