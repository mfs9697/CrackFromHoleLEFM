# Step 64: one matched physical EDI on the saved Step63 field

## Purpose

Step64 performs exactly one physical interaction-integral extraction on the already solved Step63 C03 displacement field.

It is deliberately a **postprocessing-only** step:

- no FEM solve;
- no mesh generation;
- no remeshing;
- no radius sweep;
- one fixed physical annulus;
- same 16-point FE-nodal-q interaction EDI used throughout the audit.

## Source field

Required local files:

- `verification/step63_calibrated_asymmetric_physical_solved.mat`
- `verification/step63_calibrated_asymmetric_cod_small_data.mat`

The Step63 checkpoint must identify:

- stage `step63_calibrated_asymmetric_physical`;
- selected mesh `C03`;
- 32,980 T3 triangles;
- 66,854 T6 nodes;
- crack length 8 mm;
- unit remote-y traction;
- minimal rigid-body anchoring;
- max calibrated neighbor ratio consistent with 1.79678451;
- no physical EDI previously performed by Step63;
- no mesh regeneration.

The Step63 COD result must carry the same candidate and calibration SHA-256 provenance as the solved checkpoint.

## Fixed EDI domain

Exactly one domain is used:

- `r_inner = 0.0008 m = 0.8 mm`;
- `r_outer = 0.65*a0 = 0.0052 m = 5.2 mm`.

The driver independently measures the Step63 tip-edge median and requires the fixed inner radius to remain outside `2*hTip`.

## Interaction integral

The unchanged production extractor is called exactly once with:

- `WeightFunction='fe_nodal'`;
- `QuadratureRule=16`;
- plane-strain flag inherited from the solved material;
- no stored GP diagnostics.

## COD comparison

Step64 does not refit COD. It reads the eight Step63 COD intercept fits from the same physical field and reports for every fit:

- COD KI;
- COD KII;
- COD ratio;
- single matched EDI ratio;
- signed COD-minus-EDI ratio difference;
- relative percentage gap.

It also reports the min, max, mean, and median COD fit ratios and whether the EDI ratio falls inside the complete COD-fit range.

This is cross-extractor evidence on **one physical mesh**. It is not yet a mesh-convergence proof.

## Idempotence

The solved Step63 checkpoint and COD file are SHA-256 hashed.

If a completed Step64 result already exists with matching hashes and the exact same fixed annulus, Step64 reuses it and does not repeat the interaction integration.

## Output

`verification/step64_matched_physical_edi_small_data.mat`

contains:

- the single EDI result;
- the COD-vs-EDI comparison table;
- source hashes;
- exact annulus;
- no-solve/no-remesh/no-sweep flags.

## Local run

On branch:

`audit/step64-matched-physical-edi`

run:

    addpath(genpath(pwd));

    R64 = main_step64_matched_physical_edi();

    disp(R64.EDI);
    disp(R64.Summary);
    disp(R64.CODcomparison);

Return the complete console output.
