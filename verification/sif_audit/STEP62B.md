# Step 62B: calibrate Step62 exterior grading without physical SIFs

## Purpose

Step62 established a strong structured-mesh baseline, but its realized maximum adjacent characteristic-size ratio was about **2.311**, above the preferred prospective design range of roughly **1.5–1.8**.

Step62B performs a deliberately small **mesh-only calibration** of the exterior grading. It does **not** tune against the physical asymmetric KII and performs **no physical FEM solve**.

The untracked Step38 solved checkpoint is no longer available locally. Step62B therefore uses the **committed exact Step62 candidate** as its immutable source geometry. This is sufficient for exterior mesh calibration because the candidate already stores the qualified physical boundary, crack topology, material, structured core, EDI radii, support IDs, source provenance, and prescribed-field qualification. No physical displacement field is needed.

The Step62 structured core is held fixed: paired radius 6 mm; explicit three-sector ring topology; six tip triangles / seven incident tip edges; exact T3/T6 reflection pairing; crack topology; radial core law; and physical geometry.

Only the exterior transition law is calibrated.

## Scientific rule

The allowed calibration variables are numerical mesh-design quantities only. Step62B must never inspect or optimize against the Step38 physical KI, physical KII, historical asymmetric KII/KI, or any anticipated physical SIF result.

Selection is lexicographic:

1. topology/geometry gates;
2. maximum adjacent size ratio <= 1.8;
3. smallest achieved ratio;
4. lower T3 count;
5. prescribed-field qualification.

The value 1.8 is a prospective engineering mesh-quality target, **not** an uncertainty estimate.

## Archived-source recovery and baseline diagnostic

The driver first opens the exact committed Step62 candidate `verification/step62_structured_graded_mesh_candidate_T3.mat`. Before regenerating any mesh it requires the archived candidate to report structural qualification, prescribed-field qualification, scientific readiness for an FEM *proposal*, and the seven-edge tip topology. The missing Step38 checkpoint is not required.

It then locates the actual adjacent element pairs responsible for the largest size ratio.

For each worst pair it records element IDs, shared edge node IDs, longest-edge sizes, ratio, crack-tip-local centroid radii, shared-edge midpoint radius, and whether each element is in the paired core or exterior.

## Parameterization added to Step62

build_step62_graded_exterior.m now accepts an optional design.exteriorCalibration structure. If absent, the historical Step62 defaults are preserved: transition length 8 mm; far-field slope 0.15; boundary-metric growth 0.30; six smoothing steps; 100 refinement passes; 25-degree minimum angle; longest-edge factor 1.65; neighbor-ratio target 2.5.

main_step62_structured_graded_mesh.m also receives `SourceCandidateFile`, `ExteriorCalibration`, `WriteArtifacts`, and `Verbose`. `SourceCandidateFile` permits checkpoint-independent rebuilding from the archived Step62 candidate; `WriteArtifacts=false` permits in-memory screening without producing a large collection of MAT/PNG files.

## Predeclared calibration grid

Exactly eight candidates are evaluated:

- transition length L = 6, 8, 10, 12 mm;
- far-field slope = 0.10 or 0.15.

All use fixed boundary-metric growth 0.25, six smoothing steps, 140 maximum refinement passes, 25-degree minimum angle target, longest-edge factor 1.65, and adjacent-size target 1.8.

Do not enlarge the sweep merely because more parameters exist.

## Two-stage qualification

### Stage A: mesh-only screening

All eight candidates are generated with RunSynthetic=false and WriteArtifacts=false. The driver records structural pass, maximum adjacent-size ratio, T3/T6 counts, paired/exterior minimum angles, q-support containment, transition exclusion from q-support, tip-fan gate, and native crack-face sampling.

Candidates failing structural gates or the 1.8 ratio target are rejected.

### Stage B: prescribed fields

Only the best **two** surviving candidates are regenerated in memory with the existing prescribed-field controls: pure mode I, pure mode II, and KI=1 with KII=1e-4.

The best synthetic-pass candidate is then regenerated once with full artifacts and saved under prefix verification/step62b_calibrated_mesh_selected with fresh source hashes and calibration provenance.

## Outputs

The calibration driver is verification/sif_audit/main_step62b_mesh_calibration.m.

It writes baseline worst-neighbor CSV, calibration CSV, synthetic CSV if reached, selected worst-neighbor CSV if a winner exists, a compact MAT report, and the normal selected Step62 PNG/CSV/MAT artifacts.

If no candidate satisfies the predeclared grading gate, the driver stops without forcing a winner and without running prescribed fields.

## Safety

Step62B contains no call to `solve_cracked_LEFM`. In archived-candidate source mode the Step62 builder loads no Step38 checkpoint at all and no physical displacement vector U. The only SIF extraction permitted is prescribed-field EDI after structural qualification.

No physical asymmetric FEM solve is authorized by this step.

## Local run

On branch audit/step62b-mesh-calibration run:

    addpath(genpath(pwd));
    R62b = main_step62b_mesh_calibration();
    disp(R62b.baselineWorstNeighbors(1:min(10,height(R62b.baselineWorstNeighbors)),:));
    disp(R62b.calibrationTable);
    if isfield(R62b,'syntheticTable'), disp(R62b.syntheticTable); end

Return the complete console output. If a candidate is selected, also return the overview and tip PNGs with prefix step62b_calibrated_mesh_selected.

Do **not** run a physical FEM solve afterward. The exact selected mesh must be reviewed first.


## Completed local calibration result — 2026-10-04

The investigator completed the full Step62B run locally using the archived Step62 candidate as the checkpoint-independent source.

### Baseline grading defect localized

The archived Step62 baseline maximum adjacent characteristic-size ratio was

[
2.3110209728.
]

The worst pair lies entirely in the far-field exterior, at crack-tip-local radii approximately (70.6) and (72.6) mm. The next-worst pairs are likewise exterior/exterior pairs, typically at radii (18)–(67) mm. Therefore the problematic grading is **not** in the structured paired crack-tip / EDI region.

### Eight-candidate mesh-only screen

All eight predeclared candidates passed every structural gate and the prospective ratio target (<=1.8).

| candidate | L (mm) | far slope | max neighbor ratio | T3 | T6 nodes | exterior min angle (deg) |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| C03 | 8 | 0.10 | **1.7968** | 32980 | 66854 | 25.072 |
| C01 | 6 | 0.10 | 1.7969 | 32338 | 65570 | 25.043 |
| C07 | 12 | 0.10 | 1.7969 | 34126 | 69152 | 25.037 |
| C04 | 8 | 0.15 | 1.7972 | 30246 | 61386 | 25.036 |
| C08 | 12 | 0.15 | 1.7972 | 31411 | 63719 | 25.052 |
| C02 | 6 | 0.15 | 1.7994 | 29593 | 60081 | 25.038 |
| C06 | 10 | 0.15 | 1.7996 | 30843 | 62579 | 25.017 |
| C05 | 10 | 0.10 | 1.7998 | 33514 | 67922 | 25.009 |

The structured paired patch minimum angle stayed fixed at (40.654^circ) for every candidate. All candidates retained complete q-support containment, excluded the exterior from the primary support, preserved the six-triangle/seven-edge tip fan, and retained adequate native crack-face sampling.

The predeclared lexicographic rule selected **C03** because it had the smallest achieved neighbor ratio. The differences among the leading ratios are very small; this selection should be interpreted as protocol-following rather than as evidence that (L=8) mm and far slope (0.10) are physically superior parameters.

### Prescribed-field qualification

Only C03 and C01 were subjected to prescribed-field qualification, as predeclared. Their results were identical to displayed precision:

- pure-I: (K_I=1.00000006621), (K_{II}=9.05395	imes10^{-15});
- pure-II: (K_I=2.10975	imes10^{-15}), (K_{II}=1.00000003838);
- tiny mixed input (K_I=1, K_{II}=10^{-4}): recovered (K_{II}=1.00000003847	imes10^{-4});
- affine interpolation error: (8.7311	imes10^{-11});
- recovery-matrix error: (7.6533	imes10^{-8});
- relative tiny-(K_{II}) error: (3.847	imes10^{-8}).

Both synthetic qualifications passed.

### Selected C03 candidate

The exact selected candidate was regenerated with full artifacts and provenance:

`verification/step62b_calibrated_mesh_selected_candidate_T3.mat`

Key properties:

- T3 triangles: **32980**;
- T6 nodes: **66854**;
- paired radius: (6) mm;
- tip median edge: (0.054025) mm;
- patch minimum angle: (40.654^circ);
- exterior minimum angle: (25.072^circ);
- maximum adjacent-size ratio: **1.79678451**;
- all topology, geometry, support-containment, Jacobian, grading, and sampling gates passed;
- prescribed-field qualification passed;
- **zero physical FEM solves** were performed.

Compared with the archived Step62 baseline, the ratio was reduced from (2.3110) to (1.7968), at the cost of increasing the mesh from 29777 to 32980 T3 triangles. This remains below the historical Step38 count of 39441 T3 triangles.

### Reporting caveat in checkpoint-independent mode

In the current console/figure helper text, some comparison labels still say `Step38 affected`, `Saved Step38`, or `Step38 original fan`. During Step62B checkpoint-independent regeneration, those labels refer to the **archived Step62 source candidate**, not to a recomputed Step38 mesh. The numerical candidate is unaffected; the labels should be corrected before final publication-quality artifact generation.

The exact selected mesh should be visually reviewed before any physical asymmetric solve is proposed.
