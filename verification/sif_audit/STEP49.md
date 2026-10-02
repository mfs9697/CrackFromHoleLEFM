# Step 49: reflection geometry and displacement-parity audit

## Why Step 48 requires a different diagnostic

The investigator completed the matched 16-point FE-nodal-q EDI at exactly the **same absolute annulus** on the original Step45 and selected refined Step47 saved FEM fields: `r_inner=0.00086750243077 m`, `r_outer=0.0026 m` (`r_outer/a0=0.65`). Both are the same physical centered-hole, zero-angle horizontal crack problem, for which the **continuum** `KII=0`.

| Quantity | Original Step45 | Refined Step47 |
| --- | ---: | ---: |
| T6 nodes | 3640 | 5054 |
| Measured hTip/a0 | 0.108438 | 0.0338367 |
| KI at identical EDI annulus | 0.3613086139 | 0.3618471493 |
| KII at identical EDI annulus | -3.167574307e-6 | +8.120223857e-6 |
| Signed KII/KI | -8.766949318e-6 | +2.244103310e-5 |

The investigator's Step48 code reported `delta(q_refined-q_original)=+3.1207982421e-5`, approximately **28.9%** of the separate fine asymmetric signal ratio `1.0795940665e-4` **as a magnitude comparison only**. The measured KI changes by only **~0.149%**, while the Mode-II spurious ratio changes **sign** and its magnitude grows **~2.56×**. These two meshes therefore **do not establish convergence for Mode II** at this scale. A matched EDI annulus rules out changes in the integration radii, but not different T6 interpolation properties, retriangulation outside the local refinement, or differences in actual FEM displacement fields.

The refined field's fully evaluable six COD fits also fail the proposed symmetry target and strongly depend on window and polynomial degree, including a sign reversal of the fitted quadratic ratio in the middle window. More COD nodes alone did not guarantee a stable near-tip Mode-II extrapolation.

## Next experiment: no-solve/no-EDI reflection and parity test

Rather than paying for another FEM solve or sweeping more EDI radii, compare *the two already solved fields directly*. The new `main_step49_reflection_parity.m` driver loads only the already saved Step48 compact matched comparison and the **original/refined** checkpoints it references. It refuses different physical crack geometry or an unrecognized refined checkpoint. There is **no** EDI call, FEM solve, mesh generation, or polynomial SIF fit.

Three separate diagnostics are reported:

1. **Native crack-face displacement parity.** The code reuses the audited `native_COD_audit` face labels, requires both actual native radial grids to match without interpolating the opposite face, and computes `(u_x^upper-u_x^lower)/(u_y^upper-u_y^lower)` at each paired native radius. It reports signed/absolute median ratios and normalized RMS tangential jumps separately in all three original fit windows (0.04–0.30, 0.08–0.30 and 0.12–0.30). These are **pointwise displacement-jump diagnostics, not extrapolated SIF ratios**, and cannot be equated directly to physical `KII/KI`.
2. **Off-face geometric reflection pairing.** Within the **same physical tip disk** `0<r<=r_outer` and **the same EDI annulus** `r_inner<=r<=r_outer`, compare each upper T6 coordinate `(x,+y)` with the nearest lower T6 coordinate `(x,-y)`. Report the number/fraction of upper nodes with an exact reflection match (within `max(1e-12,1e-8*a0)`) and the median nearest reflected-node distance normalized by the actual mesh tip-edge median. Only on **genuinely paired nodes**, report the even `u_x` and gauge-corrected odd `u_y` displacement-parity residuals normalized by representative native Mode-I opening. If there are too few exact geometric pairs, the parity result is **unevaluable** and appears as `NaN`, not zero. Nearest-neighbor pairing statistics alone do not establish element-to-element reflection.
3. **Entire tip-adjacent T3 triangle pairing.** Identify the complete upper and lower T3 tip fans. Reflect **all three vertices** of each upper triangle across the crack line and seek a distinct lower triangle, allowing vertex permutations. Report the best worst reflected-vertex mismatch, normalized by `hTip`, and whether **all** upper and lower tip triangles admit a complete reflection assignment. Thus equal upper/lower 3/3 fan counts can still fail this stronger test.

A small geometric test cannot prove the FEM solution physically converges: even a reflection-paired local tip fan does not establish a reflection-paired **full domain mesh**, nor does it control all T6 singular-field interpolation or EDI discretization errors. The parity audit is intended to decide whether deliberately generating a reflection-paired mesh is warranted, or whether the next isolated check should instead be an exact-nodal/Gauss replay **on the selected refined mesh** at the same original annulus.

### One local MATLAB run

Pull `sif-asymmetric-mesh-audit` using GitHub Desktop, then:

```matlab
addpath(genpath(pwd));
R49 = main_step49_reflection_parity();
disp(R49.faceWindows);
disp(R49.offFaceMirror);
disp(R49.tipFanMirror);
```

The driver locates `verification/step48_refined_matched_edi_comparison_small_data.mat` and reads its original and refined checkpoint paths. It writes **only compact tables** to `verification/step49_reflection_parity_small_data.mat`; the saved FEM meshes and solutions remain unchanged. There are no new EDI domains or FEM solves.

**Stopping rule:** interpret the face tangential jumps, percentage of geometrically pairable nodes, and *complete-triangle* tip-fan symmetry for **both meshes**. Do not attribute the Step48 EDI sign reversal to a particular source solely because one of these tests reports a mismatch; it establishes a potential mechanism to test, not numerical error causality.


## Investigator's measured Step 49 results: missing full mesh reflection

The investigator executed `main_step49_reflection_parity` on the two **already solved** Step45/47 control meshes using the exact same saved Step48 annulus. No new FEM solve, mesh, EDI or COD polynomial extrapolation occurred.

**Geometric mirror test:**

| Property | Original Step45 | Refined Step47 |
| --- | ---: | ---: |
| T6 nodes | 3640 | 5054 |
| Native upper/lower crack-face nodes | 12 / 12 | 72 / 72 |
| Native upper/lower coordinate mismatch [m] | 0 | 2.77556e-17 |
| Tip-adjacent T3 upper/lower | 2 / 3 | 3 / 3 |
| Best worst *complete upper-to-lower assignment* reflected-tip-triangle vertex mismatch [m] | 2.775110e-5 | 1.340087e-5 |
| Mismatch normalized by actual hTip | 0.063979 | 0.099011 |
| Full reflected tip-fan match? | false (also unequal counts) | false |
| Upper off-face T6 nodes inside tip disk / exactly mirrored | 89 / 0 | 494 / 0 |
| Upper off-face T6 nodes inside matched EDI annulus / exactly mirrored | 71 / 0 | 360 / 0 |
| Median nearest off-face mirror-coordinate distance / hTip, full disk | 0.042535 | 0.21773 |
| Median nearest off-face mirror-coordinate distance / hTip, annulus | 0.04418 | 0.24128 |

The strict off-face exact-match tolerance is `max(1e-12,1e-8*a0)` m; the only exact opposite-face pairing is on the coincident crack faces, not in the sampled interior. Therefore **off-face displacement parity was NOT evaluated** (appropriately `NaN` in both meshes). Zero exact off-face node matches is geometric evidence of a **non-reflection-paired mesh**, not evidence of failed actual displacement parity at nonmatching spatial points. Even the refined 3/3 tip-fan triangles are not fully mirrored: its absolute best assignment mismatch is smaller but **relative to hTip is larger** than in the coarse fan. For the coarse 2/3 fan, the assignment mismatch characterizes its best possible upper-to-subset-of-lower assignment; a complete bijective reflection is impossible because the counts differ.

**Raw actual-FEM crack-face parity, NOT a fitted SIF ratio:**

| Window r/a0 | Original native points | Original median signed jumpX/openY | Refined native points | Refined median signed jumpX/openY | Refined RMS jumpX / RMS openY |
| --- | ---: | ---: | ---: | ---: | ---: |
| 0.04–0.30 | 4 | +1.5427e-4 | 19 | -6.2179e-5 | 1.2754e-4 |
| 0.08–0.30 | 3 | -1.6581e-6 | 16 | -5.6655e-5 | 7.0047e-5 |
| 0.12–0.30 | 2 | -2.7711e-5 | 13 | -5.5841e-5 | 6.1084e-5 |

In the narrowest refined window `0.12–0.30 a0`, the RMS residual of the **gauge-removed upper/lower u_y parity** divided by representative opening is `4.0463e-6`, while the RMS tangential jump/opening is `6.1084e-5`. The finite native tangential jump is therefore a real property of the *numerically solved field* (not opposite-face radial interpolation error on these geometrically matched native faces), but it is not a validated Mode-II SIF. The sign of the refined native tangential ratio in these windows is **negative**, while Step48 refined 16-point EDI at the 0.65 domain is **positive** `+2.244103310e-5`: these are *different estimands* and the sign difference reinforces the need to scrutinize extractor/field errors rather than declare physical Mode II.

**What we know and do not know:** The physical continuum control has exact `KII=0`. Both FEM meshes demonstrably lack full nodal/element reflection pairing, and the refined actual FEM field exhibits a persistent native tangential crack-face jump. But this test alone **cannot prove that mesh asymmetry caused the EDI sign change**, because off-face physical displacement parity cannot be evaluated at exactly corresponding nodes on these nonmatching meshes.

### Next single, no-solve discriminating test: prescribed pure-I replay on the refined mesh

Before building/solving a deliberately reflection-paired mesh, determine how the **existing selected refined T6 mesh** treats a *known prescribed* pure Mode-I Williams displacement field. This mirrors Step45b on the coarse mesh and uses the **same absolute Step48 EDI annulus** (`r_outer/a0=0.65`, original `r_inner=0.00086750243077 m`). The existing `main_step45_coarse_exact_pureI_replay` accepts any verified zero-angle Step45/Step47 checkpoint and imports the actual EDI domain from the `R48.refined` output; its function name and print header remain historical. **Use a distinct SavePrefix** to avoid colliding with the coarse-mesh exact-EDI cache, which correctly rejects a different checkpoint.

The investigator already has `R48` in MATLAB from Step48. After pulling the updated audit branch (documentation only), run:

```matlab
dataDir = fileparts(R48.refined.checkpointPath);
R50 = main_step45_coarse_exact_pureI_replay( ...
    R48.refined, ...
    'ROuterOverA0',0.65, ...
    'SavePrefix',fullfile(dataDir,'step50_refined_exact_pureI'));
disp(R50.table);
```

If `R48` was cleared, restore it **without recomputing Step48**:

```matlab
load(fullfile(fileparts(P47.checkpointPath), ...
    'step48_refined_matched_edi_comparison_small_data.mat'),'R48');
```

This performs **one exact-nodal, 16-point FE-nodal-q EDI** on the already saved refined mesh after verifying its exact COD jumps. It generates **no new mesh and performs no FEM solve**. Compare its **signed** pure-I→II leakage and KI error to the previously measured coarse exact-nodal result (`KII/KI ≈ -3.88687e-5` at 0.65) and the actual-field refined EDI residual `+2.24410e-5`. Neither the exact-nodal difference nor the actual-FEM residual is a validated additive physical correction. After seeing this result, decide whether to perform a similarly matched **single exact Gauss-point EDI** on the refined mesh (existing Step45c driver, another distinct cache prefix) or move directly to a symmetry-paired meshing strategy with separate explicit authorization for any additional FEM solve.


## Measured Step 50: exact pure-I nodal replay on the selected refined mesh

The investigator completed the single prescribed-field Step50 control on the **already solved** Step47 refined T6 mesh with no new mesh/FEM solve. It used precisely the saved Step48 annulus: `r_outer/a0=0.65`, `r_inner=0.00086750243077 m`, unchanged FE-nodal q and 16-point quadrature. The replay has its own `step50_refined_exact_pureI` cache, separate from the coarse-mesh Step45b results.

The exact prescribed pure-I COD self-check on **72 native crack-face nodes** recovered unit `KI` with maximum error `2.220e-16` and maximum spurious `KII=7.125e-33`. Yet the interpolated nodal displacement field gave **`KI_EDI=1.000061439`**, **`KII_EDI=+3.379579752e-5`** and signed **`KII/KI≈+3.37937e-5`**.

### Same-annulus signed comparison across fields and meshes

| Source of actual field supplied to EDI | Original coarse mesh | Selected refined mesh |
| --- | ---: | ---: |
| Actual solved FEM signed `KII/KI` | -8.766949318e-6 | +2.244103310e-5 |
| Prescribed exact pure-I **nodal** Williams field signed `KII/KI` | ≈-3.88687e-5 | ≈+3.37937e-5 |
| Prescribed exact pure-I Williams field evaluated directly **at Gauss points** | ≈-1.147788e-7 | **+6.9597e-10 (measured Step51)** |

**Finding:** the exact-nodal field is physically pure Mode I, so its nonzero extracted Mode II is numerical. It changes sign between exactly the same two meshes as the actual-FEM EDI result does. This is strong evidence of mesh-dependent **exact-nodal interpolation/EDI contamination** with the same directional change as the FEM signal, especially since both exact nodal fields passed analytical COD self-checks. However, **sign co-variation does not establish what fraction of the actual FEM residual was caused by this mechanism**. The exact nodal test isolates a *prescribed leading Williams field*, not the numerically computed equilibrium displacement (including mesh-dependent higher-order behavior). Its leakage is ≈4.43× the coarse actual FEM residual and ≈1.51× the refined actual FEM residual in magnitude. These are descriptive ratios, **not correction coefficients**. The refined native actual-FEM tangential crack-face jump remains negative over the main fit windows while its EDI ratio is positive: direct opening ratios and EDI-extracted crack-tip SIFs are different estimands.

### Completed single matched no-solve test: exact Williams **Gauss-point** replay on the refined mesh

The **existing** `main_step45c_exact_gauss_isolation` driver accepts the Step50 `R50` output and its one measured `r_outer/a0=0.65` domain. It evaluates the exact pure-I Williams actual field *directly at Gauss points* (diagnostic `AnalyticActualK=[1,0]`), while keeping the same refined T6 mesh, same physical absolute annulus, FE-nodal-q weighting, 16-point integration, elasticity, and auxiliary mode convention. This leaves only the nodal-versus-direct-Gauss representation changed within the prescribed-field test. It performs **ONE EDI integral**, no FEM solve or mesh generation. The historical coarse-mesh Step45c at the same annulus returned signed leakage ≈`-1.147788e-7`.

Use a **different prefix** to protect all coarse and refined nodal replay caches. With `R50` in MATLAB's workspace:

```matlab
dataDir = fileparts(R50.checkpointPath);
R51 = main_step45c_exact_gauss_isolation( ...
    R50, ...
    'ROuterOverA0',0.65, ...
    'SavePrefix',fullfile(dataDir,'step51_refined_exact_gauss'));
disp(R51.table);
```

If `R50` was cleared, restore its measured compact output first; do not repeat Step50:

```matlab
load(fullfile(fileparts(R48.refined.checkpointPath), ...
    'step50_refined_exact_pureI_small_data.mat'),'Out');
R50 = Out;
```

**Interpretation rule:** If refined exact-Gauss leakage is much smaller than refined exact-nodal leakage, the mesh-dependent exact-nodal interpolation mechanism is demonstrated on **both** symmetric meshes. If it remains comparable on the refined mesh, investigate refined FE-nodal q/quadrature or analytical/convention differences before inferring interpolation dominance. Neither result validates an error correction or uncertainty bound for the separate much finer asymmetric tiny Mode-II calculation.


## Investigator's measured Step 51: exact Gauss-point control on refined T6 mesh

The investigator completed the prepared **no-solve** Step51 replay on the same saved 5,054-node selected refined symmetric T6 mesh, at the **exact** previously matched absolute annulus (`r_inner=0.00086750243077 m`, `r_outer/a0=0.65`), with the same FE-nodal q and 16-point integration.

The analytical actual field was evaluated **directly at Gauss points**, bypassing T6 interpolation of prescribed nodal exact Williams displacements. It recovered **`KI=1` at displayed precision**, **`KII=+6.959687749e-10`** and therefore signed `KII/KI≈+6.9597e-10`, relative to the same-mesh prescribed **exact-nodal** result `KI=1.000061439`, `KII=+3.379579752e-5`, signed `KII/KI≈+3.37937e-5`.

The exact-nodal → exact-Gauss reduction in apparent spurious Mode II is approximately **48,600×** on the refined mesh. On the original coarse mesh, the analogous measured reduction at the same annulus was **~339×** (`-3.88687e-5 → -1.147788e-7`). This establishes that **ordinary T6 interpolation of the prescribed singular leading Williams field** is sufficient to generate mesh-dependent, sign-changing artificial Mode II in the exact-nodal EDI replay, whereas direct analytical Gauss-point evaluation of that same field produces extremely small leakage on both meshes.

**Important distinction:** the original and refined actual *computed FEM* EDI ratios (`-8.76695e-6` and `+2.24410e-5`) also change sign, but this co-variation is not a numerical decomposition of the computed-field error. Prescribed analytical leading Williams fields and actual solved elastic fields are different inputs, and the exact-Gauss diagnostic shares analytical conventions with the auxiliary fields; it is **not** an independent validation of physical SIF accuracy or a transferable correction to the asymmetrically cracked plate's very small `KII`.

The matched prescribed-field EDI testing is complete for both symmetric meshes. Additional exact-field radius sweeps are not justified at present. Proceed instead to [Step 52: interpolated reflected-point parity](STEP52.md) on the two ALREADY SAVED actual FEM fields, with an affine-field interpolation self-check and synthetic exact pure Mode-I interpolation control. Do **not** build or solve a new FEM mesh yet.
