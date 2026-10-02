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
