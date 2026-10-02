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

## Remaining distinct tests

1. If COD exact-field replay passes but COD/EDI disagree on actual `U`, consider finite-element representation and higher-order terms, not just extractor implementation.
2. If exact EDI has nontrivial unit-I→II leakage at the relevant scale, audit EDI auxiliary derivatives, q construction, and integration on the same mesh.
3. A separate **symmetric zero-`KII` FE solution** is ultimately required to measure method-wide leakage from numerical FEM approximation. That would require a separately approved FEM solve.
