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


## Completed local Step64 result — 2026-10-04

Step64 ran on the recovered Step63R physical field with zero additional FEM solves.

Fixed matched annulus:
- r_inner = 0.8 mm;
- r_outer = 5.2 mm = 0.65 a0;
- tip median edge = 0.0540246508 mm.

Single 16-point FE-nodal-q interaction EDI:

[
K_I = 0.43785,qquad
K_{II} = 4.6547	imes10^{-5},qquad
K_{II}/K_I = 1.0631	imes10^{-4}.
]

The eight Step63 COD-fit ratios span 1.0416e-4 to 1.0583e-4, with mean 1.0523e-4 and median 1.0544e-4. The EDI value lies 0.456% above the closest/highest COD fit, 0.814% above the COD median, and 1.014% above the COD mean. Individual COD-vs-EDI gaps range from 0.456% to 2.027%.

The four quadratic COD fits are the closest subset: their relative gaps to EDI are approximately 0.645%, 0.566%, 0.462%, and 0.456%.

Interpretation: COD and interaction EDI independently identify the same positive Mode-II signal of order 1.05e-4 on the same deliberately paired physical mesh. This is strong cross-extractor evidence, but it remains a one-mesh result rather than a mesh-convergence proof.
