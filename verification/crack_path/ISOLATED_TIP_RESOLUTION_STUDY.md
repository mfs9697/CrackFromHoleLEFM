# Independent core and exterior resolution

Implemented on `tip-refinement-fixed-geometry`, 2026-10-09. No physical FEM solve was performed during implementation or verification.

`ExteriorScale` is optional in `main_stage2_embed_scaled_core_full_domain_theta0`, `qualify_incremental_crack_candidate`, and `run_incremental_crack_path`. Omission, or an empty value, resolves to `CoreScale`. Existing drivers therefore retain their historical coupled mesh family, including the historical half-size and double-size studies.

The structured patch still uses `CoreScale`. Both exterior builders now use the resolved `ExteriorScale` for the requested radial target-size law, retained-crack subdivisions, ring seeds, boundary-layer spacing and boundary metric. The original C1 target-size formula is shared without changing its arithmetic. The far-field cap, transition length, calibration defaults, physical geometry, EDI radii and all numerical gates retain their existing definitions.

For `CoreScale=2, ExteriorScale=1`, the core has 3,318 T3 elements instead of the reference 12,678, while the **requested exterior size law** is exactly the reference law. Its boundary target is `hBase + slope*rCore`, with the reference default transition and far-field cap. The assembled exterior triangles can differ because they conform to a different core-interface subdivision; this does not imply a change to the requested exterior law.

## Provenance and reuse

New candidates record `structuredDesign.exteriorScale`, `exteriorDesign.exteriorScale` and `exteriorMeshControls.exteriorScale`. Qualification and physical summaries append `exterior_scale`; P1 summaries also record `core_scale`. Path summaries include both scales, and path resume states retain both in their mesh controls. New physical checkpoints record a canonical exterior identity alongside the core scale.

Candidate reuse, P1 compact-result reuse, physical checkpoint reuse and path resume compare the resolved exterior scale as well as the existing family controls. A legacy record without an exterior scale means its **own recorded core scale**, not an unconditional exterior scale of one. Legacy physical fields still require the existing exact mesh and physics checks. They cannot authorize reuse by an isolated family whose exterior scale differs from its core scale. Labels and historical reference flags do not replace numeric family identity.

The existing `isReferenceProductionExterior` flag retains its historical meaning: reference transition, cap and calibration parameters, regardless of the coupled scale. This preserves the old double-size trajectory driver. New `isReferenceRequestedExteriorLaw` metadata additionally requires an exterior scale of one. Compatibility uses numeric controls and the parent core provenance, not either descriptive flag.

## Study driver

From the repository in MATLAB:

```matlab
addpath(genpath(pwd));
D = main_isolated_tip_resolution_study(); % qualification only

% Execute the study's physical solves explicitly:
D = main_isolated_tip_resolution_study('AllowPhysicalSolves',true);
```

Each fresh invocation creates a timestamped directory under `verification/crack_path/isolated_tip_resolution_study_*`. An explicit `OutputDir` must not already exist. To continue a particular interrupted study, pass its directory as `ResumeStudyDir` instead of `OutputDir`. That explicit continuation validates and reuses its saved qualified candidates and accepted compact/field results. See [the native-face fingerprint fix and resume trace](ISOLATED_NATIVE_FINGERPRINT_FIX.md). Historical output directories and archives are never selected automatically.

The driver reads the exact committed Stage-I state and accepted reference vertices. It verifies that the requested P17/P21/P22/P23 states are accepted and have segment lengths compatible with the frozen increment. It then:

1. Qualifies each fixed reference path at `CoreScale=2` and `0.5`, always with `ExteriorScale=1`; when physical solves are enabled, solves those eight prescribed states and saves their EDI/COD/MTS results and reference comparisons in separate family directories.
2. Qualifies and, when enabled, solves a fresh P1 at `CoreScale=2, ExteriorScale=1` in its own seed directory.
3. Seeds an independent trajectory from that P1's own KI/KII and propagates through P23 with the same scale pair. Fixed reference paths and their results do not seed this propagation.

Synthetic replay is always enabled. All qualification, solver, residual, EDI, COD and MTS validity gates and tolerances remain active in the existing routines. As in the existing alternative-family trajectory driver, `RegressionGates=false` disables only the hard-coded reference-family P1/P2 numeric anchors; an independent alternative mesh family cannot be required to reproduce those reference numbers. No gate formula or tolerance is edited. If an ordinary gate or the physical-clearance stop prevents completion through P23, the driver preserves results and reports the stop rather than relaxing the gate.

The top-level study MAT and CSV summaries are updated after each accepted fixed-state result. Every full physical field and compact result has a separate file. The independent trajectory also writes the existing atomic `path_run_state.mat`; it can subsequently be continued using the production driver's resume options with `CoreScale=2, ExteriorScale=1` and its own P1 seed.

## Verification

MATLAB R2023a checks performed without physical solves:

| Check | Result |
| --- | --- |
| `test_isolated_tip_resolution` | 32 portable checks passed; 2 investigator-local archive checks passed |
| Existing `test_crack_reproducibility` | 15 portable and 22 archive checks passed |
| Original-source exterior comparison | Both builders matched original `dbd66d9` source bitwise for omitted exterior scale at core scales 0.5, 1 and 2; coordinates, connectivity, boundary IDs/edges and prior metadata matched |
| Full P1 qualification, core 2 / exterior 1 | Structural and prescribed-field gates passed |
| Full fixed P17 qualification, core 2 / exterior 1 | Structural and prescribed-field gates passed |
| Full fixed P23 qualification, core 0.5 / exterior 1 | Structural and prescribed-field gates passed |
| Existing double-size independent-trajectory preflight | Passed with omitted exterior scale resolving to 2; historical reference flag preserved and new effective-law flag false |
| MATLAB code analysis and Git whitespace check | No syntax errors; existing performance/unused-variable advisories remain |

The portable suite verifies the historical size-law arithmetic, exact omitted-versus-explicit mesh equality, changed core fingerprints with unchanged reference requested exterior law, candidate-cache mismatch rejection, new and legacy checkpoint mismatch rejection, new and legacy resume mismatch rejection, P1 compact-cache mismatch rejection and the study's output-overwrite guard. The archive checks exercise scale rejection through both public physical solver entry points with solves disabled.

These verification results describe the original scale-separation implementation. Representative prescribed-field qualifications establish mesh and synthetic eligibility, not physical convergence results. The later physical study, interrupted P1 fingerprint diagnosis, and continuation evidence are recorded in [ISOLATED_NATIVE_FINGERPRINT_FIX.md](ISOLATED_NATIVE_FINGERPRINT_FIX.md).
