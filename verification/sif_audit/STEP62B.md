# Step 62B: calibrate Step62 exterior grading without physical SIFs

## Purpose

Step62 established a strong structured-mesh baseline, but its realized maximum adjacent characteristic-size ratio was about **2.311**, above the preferred prospective design range of roughly **1.5–1.8**.

Step62B performs a deliberately small **mesh-only calibration** of the exterior grading. It does **not** tune against the physical asymmetric KII and performs **no physical FEM solve**.

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

## Baseline diagnostic

The driver first opens the exact archived Step62 candidate verification/step62_structured_graded_mesh_candidate_T3.mat and locates the actual adjacent element pairs responsible for the largest size ratio.

For each worst pair it records element IDs, shared edge node IDs, longest-edge sizes, ratio, crack-tip-local centroid radii, shared-edge midpoint radius, and whether each element is in the paired core or exterior.

## Parameterization added to Step62

build_step62_graded_exterior.m now accepts an optional design.exteriorCalibration structure. If absent, the historical Step62 defaults are preserved: transition length 8 mm; far-field slope 0.15; boundary-metric growth 0.30; six smoothing steps; 100 refinement passes; 25-degree minimum angle; longest-edge factor 1.65; neighbor-ratio target 2.5.

main_step62_structured_graded_mesh.m also receives ExteriorCalibration and WriteArtifacts. WriteArtifacts=false allows candidate screening in memory without producing a large collection of MAT/PNG files.

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

Step62B contains no call to solve_cracked_LEFM. The Step62 builder still loads no physical displacement vector U. The only SIF extraction permitted is prescribed-field EDI after structural qualification.

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
