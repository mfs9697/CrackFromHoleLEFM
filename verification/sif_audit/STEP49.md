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
