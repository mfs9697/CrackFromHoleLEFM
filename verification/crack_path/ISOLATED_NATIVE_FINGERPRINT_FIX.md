# Isolated P1 native-face fingerprint and study continuation

The interrupted study is `isolated_tip_resolution_study_20261009T163513088`. Its eight fixed reference-state physical calculations passed. The interruption occurred during P1 **before stiffness assembly or any P1 physical solve**, at the pre-solve native sampling check in `main_stage2_theta0_physical_solve`.

## Creation and verification trace

1. `build_stage2_scaled_audited_core.m:80` constructs distinct upper/lower T3 crack-axis chains. At line 90 it collects their T6 midsides, orders nodes by distance from the tip, and exposes the canonical face nodes in `Core.crack.upperT6/lowerT6` at line 125. This topology depends on `CoreScale`.
2. `main_stage2_embed_scaled_core_full_domain_theta0.m:276` combines the core chains and exterior crack-face nodes. Its full crack metadata therefore also depends on the exterior subdivision. `local_topology_audit`, called at line 342, uses `native_COD_audit` to label T6 face midsides and measure the paired native grid.
3. Qualification counts samples in the four unchanged COD windows at line 370. The scale-2 canonical counts are `[19;28;23;18]` (line 175), tested by `gates.nativeSamplingExact` at line 410 and saved in `candidate.coreMeshControls.expectedNativeSamples` at line 557. These windows lie entirely inside the core. The old qualification did **not** promise a 70-node complete crack face.
4. `verification/sif_audit/native_COD_audit.m:17` classifies faces by topology; at line 37 it selects the non-tip straight-crack nodes and at lines 42–44 sorts/uniquifies the abscissae. Its `diag.nUpper/nLower` count the **complete tip-to-mouth T6 faces**, including the exterior portion. The COD formulas and this extractor are unchanged.
5. Before the fix, `main_stage2_theta0_physical_solve.m` at commit `a89be2e` rebuilt T6, called that same extractor in the global frame (lines 196–209), and compared its complete-face counts with `coreFP.expectedFaceNodes`. The private scale-2 branch hard-coded `faceNodes=70` (line 707). That was a fingerprint of the historical *coupled* core-2/exterior-2 full face, not the core alone. Its scale-1 branch likewise hard-coded the historical full-face count 138.

## First differing quantity

Read-only extraction from the exact saved qualified P1 candidate gave:

| Quantity | Qualified isolated candidate / solver reconstruction | Historical solver expectation | Interpretation |
| --- | --- | --- | --- |
| Full upper/lower native face nodes | 80 / 80 | 70 / 70 | **First failing quantity**; full-face subdivision includes the independently finer exterior |
| Core native face nodes, excluding tip | 60 / 60 | Canonical scale-2 core: 60 / 60 | Unchanged |
| COD window counts | `[19,28,23,18]` | `[19,28,23,18]` | Unchanged; column-vector representation also matches |
| Upper/lower grid mismatch | 0 m | At most `1e-12` m | Unchanged pairing |
| Core abscissae versus canonical core | Maximum difference `2.298508605669269e-17` m | Same core grid | Floating-point local/global frame representation only |
| Full upper/lower T3 chain nodes, including tip | 41 / 41 | Implies exactly `2*(41-1)=80` non-tip T6 nodes per face | Consistent T3/T6 representation |

The full face contains 60 core nodes plus 20 exterior nodes. The coupled historical double-size family contains the same 60 core nodes plus 10 exterior nodes. `ExteriorScale=1` deliberately changes the exterior's mesh representation of the same physical crack; it does not change the frozen material, loading, crack length, mouth/tip, or the core/EDI geometry. No ordering mismatch caused this failure. Qualification and the physical entry point recover the same 80-node grid from the saved mesh.

## Correction

`stage2_native_sampling_fingerprint.m` separates the two identities:

- The expected **full-face grid** is derived independently from each declared T3 tip-to-mouth boundary chain and its exact endpoint midsides. Every chain edge must be a singly incident domain boundary. T6 full-face counts remain exact, and its abscissae must agree with that T3-derived grid.
- The **core-face grid** is checked against a separately generated canonical structured core at the recorded `CoreScale`; full-face count changes cannot conceal a core change.
- The original four native sampling counts, face pairing tolerance and all physical/material, EDI, residual, solver and MTS gates remain in force. No acceptance tolerance is relaxed.

New qualification records this separate fingerprint in the candidate and compact qualification. The physical solver verifies it before assembly and records it in the field checkpoint and compact result. Existing qualified candidates lacking the new field are checked directly against their T3 boundary chains and the independent canonical core. They need not be remeshed, requalified, or assigned a manually patched expected count.

The post-solve native pairing gate now uses the independently verified full-face counts, while its exact COD fit sample counts remain the canonical core-window counts. The erroneous hard-coded full-face counts are removed from the core fingerprint.

## Resume without repeating accepted work

The study driver now accepts:

```matlab
Riso = main_isolated_tip_resolution_study( ...
    'AllowPhysicalSolves',true, ...
    'ResumeStudyDir',fullfile(pwd,'verification','crack_path', ...
        'isolated_tip_resolution_study_20261009T163513088'));
```

An explicit resume validates the study controls, exact saved candidate paths, material, passed qualification/synthetic gates and mesh-family identity. For accepted fixed-state compact results, `validate_isolated_fixed_result_cache` verifies the original T3/T6 field geometry, physics, solver provenance and EDI/core metadata without loading U, repeating PCG, or repeating EDI/COD postprocessing. Missing compact results can still be recovered through the existing physical-field checkpoint architecture. P1 uses its saved qualified candidate; accepted P1 compact/field data are reused when available. Subsequent trajectory continuation uses the existing atomic path-state resume checks.

Fresh invocations still require a new output directory. Explicit continuation updates only the selected study's progress and new continuation results. Original qualified candidate files and accepted fixed-state field/compact files are retained.

## Verification evidence

`test_stage2_native_fingerprint` checks the actual saved isolated P1, reversed face-node ordering, complete T3 node renumbering/T6 ID changes, paired T6 midside corruption, paired T3/core corruption, stored fingerprint corruption, the public solver's explicit no-solve guard, the historical 70-node coupled candidate and mismatched fixed-result cache controls. All eight actual accepted fixed-state caches are independently validated. `test_isolated_tip_resolution` is rerun to preserve scale/default and checkpoint incompatibility regressions.

The saved P1 passes the pre-solve fingerprint and reaches the explicit physical solve guard with `AllowSolve=false`. Nine fingerprint/ordering/corruption/public-entry regressions and eight accepted fixed-cache checks passed. The existing scale-separation suite also passed all 32 portable and two archive checks. A fresh prescribed-field P1 qualification passed and emitted the new 80-node full-face / 60-node core fingerprint. MATLAB code analysis found no syntax or uninitialized-variable issues after the correction; existing unused-variable advisories remain.

## Completed physical continuation

The selected study completed through **P23**, with all eight fixed-state results and all 23 independent physical states passing their existing gates. The driver reused nine saved qualified candidates and retained all eight accepted fixed physical results without repeating either PCG or their postprocessing.

P1 was solved once during the authorized continuation and its field was checkpointed. A missing COD-window variable introduced during this refactor stopped the first postprocessing attempt. Restoring `windows=nativeFP.windows` preserved the exact original windows. The subsequent attempt reused that field (`R.newSolve=false`) and passed every P1 physical gate; no P1 physical solve was repeated. P2–P23 each required one new physical solve and passed all qualification, solver, residual, EDI, COD and MTS checks.

Read-only compact-result/state verification confirmed:

| Completion evidence | Result |
| --- | --- |
| Atomic path state | `completedPhysicalSegments=23`; core scale 2 / exterior scale 1 |
| Stop reason | `max_segments_reached` |
| Total physical crack length | `92.00000000000003 mm` (floating-point representation of 92 mm) |
| P1 native fingerprint | Full faces 80/80; core 60; windows `[19,28,23,18]` |
| P2–P23 COD fit sample counts | Exact `[19,19,28,28,23,23,18,18]` for every state |
| P2–P23 EDI support | 2,976 elements for every state |
| Largest P2–P23 true relative residual | `9.987669875849209e-11` |
| Original saved MAT files | SHA-256 unchanged for all 34: eight fixed candidates, eight qualifications, eight physical fields, eight compact physical results, plus the original P1 candidate and qualification |

The completed study remains in `verification/crack_path/isolated_tip_resolution_study_20261009T163513088`. Its updated study MAT contains both fixed-state comparisons and the independent trajectory; the trajectory directory contains the accepted compact results, physical checkpoints and final atomic resume state. `independent_validation_summary.csv` records the read-only completion checks. No branch was merged, and no physical formula or acceptance tolerance was changed.
