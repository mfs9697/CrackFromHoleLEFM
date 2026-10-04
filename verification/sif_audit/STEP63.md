# Step 63: one physical solve on the calibrated Step62B C03 mesh

## Purpose

Step63 is the first physical asymmetric FEM calculation on the deliberately structured, reflection-balanced and calibrated mesh family.

The investigator has explicitly authorized exactly one solve on the selected Step62B C03 mesh.

Step63 deliberately does **not** run a physical interaction EDI. It checkpoints the solved field first and then evaluates native crack-face COD only. Physical EDI is reserved for a separate later step after the displacement/COD field is reviewed.

## Exact mesh source and recovery

The preferred local Step62B outputs are:

- `verification/step62b_calibrated_mesh_selected_candidate_T3.mat`
- `verification/step62b_mesh_calibration_small_data.mat`

These files are generated artifacts and can disappear when switching branches/clones. Step63 therefore has a deterministic recovery path.

If the selected C03 candidate is missing, Step63 regenerates it from the committed archived Step62 baseline:

- `verification/step62_structured_graded_mesh_candidate_T3.mat`;
- fixed C03 parameters (L=8) mm, far slope (0.10), boundary-metric growth (0.25), neighbor target (1.8);
- unchanged Step62 structured core.

The recovery performs **no physical FEM solve**. It reruns only mesh construction and the already accepted prescribed-field controls, and requires reproduction of the recorded C03 counts, grading result, minimum angles and synthetic pass before proceeding.

If the compact Step62B calibration MAT is also missing, Step63 reconstructs only the minimal calibration provenance from the exact qualified C03 candidate and saves it locally.

It then requires the selected candidate to be C03 with:

- transition length 8 mm;
- far-field slope 0.10;
- neighbor-ratio target 1.8;
- 32,980 T3 triangles;
- 66,854 T6 nodes after deterministic upgrade;
- maximum adjacent-size ratio consistent with the completed Step62B result (~1.79678451);
- six tip triangles / seven incident tip edges;
- complete structural pass;
- prescribed-field qualification pass;
- scientific-ready flag from Step62B.

The final candidate and compact calibration provenance files—whether pre-existing or deterministically recovered—are SHA-256 hashed. The hashes are stored in the solved checkpoint and checked before any existing Step63 result is reused.

## Physical setup

No nominal hole geometry or mesh generator is invoked.

The full physical domain is the exact saved candidate geometry.

The non-geometric physics is cross-checked against `cfg_hole_initiation.m`:

- material `E`, `nu`, plane-strain flag;
- remote tension in y;
- unit nominal traction;
- minimal rigid-body anchoring.

The outer plate extents are measured directly from the selected mesh and required to remain x=[0,0.30] m and y=[-0.10,0.10] m.

## Solve policy

`AllowSolve=false` by default.

If a valid Step63 checkpoint already exists for the exact candidate/calibration hashes, it is reused and no solve is repeated.

Only when the investigator explicitly runs with `AllowSolve=true` and no valid checkpoint exists does the driver call `solve_cracked_LEFM` exactly once.

The solver output is required to reproduce the exact candidate T3 and deterministic T6 coordinates/connectivity. Any mesh change aborts the step.

The solved checkpoint is written first to an `.incomplete.mat` temporary file and atomically moved into place only after the solve has completed.

Default checkpoint:

`verification/step63_calibrated_asymmetric_physical_solved.mat`

The checkpoint stores mesh, U, material, crack, a0, compact physical configuration and provenance metadata. Large transient stiffness/stress objects are not retained after checkpointing.

## COD-only postprocessing

After checkpointing, Step63 evaluates the already audited `native_COD_audit` on the physical displacement field.

Four windows are reported:

- 0.04–0.20 a0
- 0.04–0.30 a0
- 0.08–0.30 a0
- 0.12–0.30 a0

with linear and quadratic intercept fits.

The output also reports raw pointwise apparent KII/KI ratios in five radial bands.

These COD fits are diagnostic. Previous audit steps showed that intercept sensitivity can be substantial, so Step63 does not treat a particular COD intercept as a validated physical KII.

## Safety / scope

Static review of the updated Step63 driver confirms:

- exactly one `solve_cracked_LEFM` call site;
- one mesh-only `main_step62_structured_graded_mesh` recovery call site;
- deterministic C03 recovery must reproduce 32,980 T3 / 66,854 T6, max ratio 1.79678451 within tolerance, minimum angles, and prescribed-field pass;
- zero `generateMesh` calls;
- zero `SIF_LEFM_interaction_EDI` calls;
- exact-mesh checks after the solver;
- explicit candidate/calibration SHA-256 provenance;
- existing checkpoint reuse instead of repeated solving;
- physical EDI explicitly marked not performed.

## Authorized local run

On branch:

`audit/step63-calibrated-physical-solve`

run:

    addpath(genpath(pwd));

    [P63,O63] = main_step63_calibrated_asymmetric_physical_solve( ...
        'AllowSolve', true);

    disp(O63.rawTable);
    disp(O63.fitTable);

Return the complete console output.

Do not run a separate physical EDI command yet. The Step63 COD result should be interpreted first.
