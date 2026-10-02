# Step 44: exact Williams-field SIF extraction replay

**Parent:** `sif-asymmetric-mesh-audit`, Steps 34–43. This change is a **verification experiment**, not a replacement of the production SIF selector.

## Motivation

Step 39–43:
- Halving tip scale from ~0.108 mm to ~0.054 mm barely changed the three EDI ratios (~0.39%) but increased the native-face COD estimates.
- Matching all 133 old nodes showed the COD shift is predominantly a change in the computed FE displacement field, **not** the extra face samples.
- Fitted COD intercepts still change with lower/upper cutoff and polynomial order. Fitting-window drift is **not** a physical uncertainty estimate.
- A single near-tip midside Mode II COD point is negative, although more distant values are positive. Do not interpret it as a physical kink sign.
- Published DEM/IIM comparisons often use **quarter-point** tip elements; our current T6 tip is **conventional**, so published benchmark accuracies cannot simply be transferred.

## Files

- `native_COD_audit.m`: the original Step-38 local COD implementation, moved unchanged numerically into a shared function. Optional sixth argument `includeFace=true` exposes crack-face classification; the default keeps Step-38 outputs compact.
- `main_step38_postprocess_tip_checkpoint.m`: now calls this shared extractor instead of a private local copy.
- `test_step44_cod_micro_mesh.m`: quick, **checkpoint-free** exact pure-I / pure-II / tiny-mixed COD recovery on a manufactured, topologically separated T6 crack-face micro mesh.
- `main_step44_exact_williams_replay.m`: loads the **already solved** Step-38 checkpoint and verifies shared-real-COD reproducibility before any analytical replay. By default, replays three exact cases in COD only. Expensive interaction-EDI work is strictly opt-in.

## Run in MATLAB

```matlab
addpath(genpath(pwd));
test_step44_cod_micro_mesh();

% Load the compact Step-38 MATLAB result (NOT Step-39 small_data).
load('step38_tip_refined_solved_results.mat','O38');

% Fast phase, no EDI integration and no FEM solve:
O44 = main_step44_exact_williams_replay(O38);
O44.CODgates
O44.CODmatrix
O44.fitTable

% Only if COD gates and regression pass:
O44 = main_step44_exact_williams_replay(O38, 'RunEDI',true);
O44.EDI
```

The driver uses `O38.checkpointPath` to find the **original solved-field** `step38_tip_refined_solved.mat`, which contains the full `mesh`, `U`, material, and crack. Step-38 **small_data** and Step-39 `compact` data alone are insufficient for exact-field replay.

The default synthetic mixed-mode amplitudes are the saved refined-mesh EDI `KI,KII` from the existing annulus closest to `r_outer/a0=0.65`. These only define a **known input signal**; they are not treated as verified physical SIFs.

### COD acceptance / reporting

- The first gate re-evaluates the previously solved **actual** field with the shared routine and compares it with `O38.nativeR` / `O38.nativeApparent`; abort if they differ.
- Pure unit I and unit II identify leakage / diagonal normalization. The mixed case tests whether that leakage is negligible relative to `KII ≈ 4.7×10⁻⁵`.
- Evaluate actual native-face points without creating artificial regression points. Retain original Step-39 fitting windows and degrees.
- Recovering an exact analytical COD field says nothing by itself about accuracy of the **computed FE** displacement field.
- Do not treat a small regression RMSE or matching EDI estimate as independent proof of physical accuracy.

### Optional EDI phase

`'RunEDI',true` computes **two** expensive 16-point FE-nodal-q interaction integrals at the same existing annulus: one exact unit-I field and one exact unit-II field, using the **existing** Step-38 T6 mesh. The recovery matrix `M_EDI` should approach identity, and its pure-I→II leakage must be compared with the ~10⁻⁴ mixed-mode ratio of interest.

The tiny-mixed EDI result is **inferred from linearity** (`M_EDI * [KI;KII]`), **not** a third independently integrated field. Each completed unit integration is saved immediately to `step44_exact_williams_edi_progress.mat` so MATLAB interruptions do not lose finished cases; the cache rejects a different checkpoint file size/timestamp or EDI annulus.

**No FEM solver, stiffness assembly, global remeshing, or geometry modification is called.** Full MATLAB runs must be completed on the machine holding the saved checkpoint; source inspection alone is not a numerical verification.

## MATLAB results supplied 2026-10-02 — Step 44 recovery passed

The investigator ran both phases in MATLAB using the existing solved Step-38 checkpoint. **These are measured outputs supplied by the investigator, not results of a GitHub Actions run.** They must not be confused with direct verification of the unknown physical SIF.

- Micro-mesh COD test: `passed=1` (16 upper and 16 lower native face samples).
- Real-field refactoring regression: **exact agreement**, `n=152`, `max dr=0`, `max dK=0` on the existing 79,769-node T6 mesh.
- Max crack-face transverse coordinate roundoff: `3.849e-18 m`; upper/lower abscissa mismatch: `1.11022e-16 m`.
- COD unit I/II/mixed replay: `||M_COD-I||_F=3.140e-16`, maximum cross-mode leakage `6.029e-18`, relative tiny-mixed Mode-II raw deviation `7.526e-14`, maximum fitted deviation `1.300e-13`. All COD gates passed.
- Exact unit-field EDI, same original `r_inner=0.0008 m`, `r_outer=0.0052 m`, 16-point FE-nodal-q:
  ```text
  M_EDI = [ +1.0000000174e+00, +6.7874710151e-09;
            -2.5354410276e-08, +1.0000000376e+00 ]
  ||M_EDI-I||_F = 4.90843e-08
  ```
- For imposed `[KI;KII]=[0.43783617;4.7268533e-5]`, the **linearity-inferred, not independently integrated** EDI mixed Mode-II relative deviation was `-0.000234814` (= `-0.0234814%`). The contribution from pure-I→II leakage is `M(2,1)*KI ≈ -1.11e-8`, or ~0.0235% of the tiny mixed Mode-II signal.
- The refined physical-field linear COD fit `0.12 <= r/a0 <= 0.30` had an ~8.57% discrepancy from the corresponding actual-field EDI ratio. The exact-field EDI cross-mode contamination at these settings is roughly **365 times smaller** in relative scale. This comparison excludes a basic exact-leading-field numerical leakage of the measured size as the main explanation; it is **not** an estimate of physical EDI or COD accuracy.

**Scientific limits:** Exact nodal field replay tests the two extraction pipelines on *prescribed leading Williams fields*, not the computed FEM equilibrium error or the higher-order terms in its actual crack-tip solution. The exact-field generator and EDI auxiliary fields share analytical conventions, which limits independence of a sign/normalization check. Only one original EDI domain was replayed with exact fields; the earlier three-domain stability is for the actual FE solution.

The huge `ratio_if_defined` for the exact pure-II case is **undefined mathematically** because the prescribed `KI=0`. Ignore that ratio and assess the individual recovered SIFs / off-diagonal leakage instead.

**Next gate:** prefer a saved, independent symmetric pure Mode-I *FEM solution* (if available) on comparable conventional-T6 topology to quantify field-induced spurious `KII` using both extractors. If none is saved, approve and checkpoint exactly one new symmetric control solve; do not blindly increase crack-tip refinement. Examine an analytically defined next Williams term separately to characterize the finite-distance COD extrapolation effect.

## Remaining distinct tests

1. If COD exact-field replay passes but COD/EDI disagree on actual `U`, consider finite-element representation and higher-order terms, not just extractor implementation.
2. If exact EDI has nontrivial unit-I→II leakage at the relevant scale, audit EDI auxiliary derivatives, q construction, and integration on the same mesh.
3. A separate **symmetric zero-`KII` FE solution** is ultimately required to measure method-wide leakage from numerical FEM approximation. That would require a separately approved FEM solve.
