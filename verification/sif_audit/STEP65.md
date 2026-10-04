# Step 65: Level-1 qualification of the calibrated C03 mesh family

## Purpose

Step65 constructs and qualifies the next spatial level of the deliberately structured asymmetric mesh family.

It is **mesh-only**. No physical displacement field is loaded, no stiffness matrix is assembled, and no physical FEM solve is permitted.

The scientific question is whether the accepted Level-0 C03 discretization admits a clean Level-1 successor with the same physical geometry and extraction support but one-half the structured tip/core length scale.

## Source and reproducibility

Step65 depends only on the committed archived Step62 candidate:

`verification/step62_structured_graded_mesh_candidate_T3.mat`

Generated Level-0 C03 or Step63 MAT files are not required.

The driver first reconstructs Level 0 in memory using the fixed C03 calibration and requires exact reproduction of the accepted baseline:

- 32,980 T3 triangles;
- 66,854 T6 nodes;
- maximum adjacent-size ratio 1.79678451 within tolerance;
- paired radius 6 mm;
- six tip triangles and seven topological tip edges.

If this baseline reproduction fails, Step65 stops before Level 1.

## Family definition

Both levels use:

- exact same physical plate/hole/crack geometry inherited from the archived source;
- paired structured radius 6 mm, fully containing the primary 0.8–5.2 mm EDI support;
- C03 exterior transition length 8 mm;
- far-field slope 0.10;
- boundary-metric growth 0.25;
- maximum adjacent-size target 1.8;
- same deterministic three-sector radial-zipper topology.

The structured size law is

`h(r) = 2^(-level) * (hBase + 0.028*r)`.

Therefore Level 1 has exactly one-half the Level-0 structured target scale.

## Level-1 structural gates

Step65 requires:

- unchanged physical geometry and physical boundary nodes;
- unchanged material;
- complete T3 and T6 reflection pairing;
- crack faces distinct and intact ligament shared;
- positive T3 areas and T6 Jacobians;
- no hanging gaps/overlaps;
- complete production and literal q-support containment inside the paired region;
- exterior excluded from the primary support;
- six-triangle / seven-edge tip fan;
- target tip scale exactly halved;
- measured tip scale approximately halved;
- unchanged 6-mm paired radius;
- maximum adjacent-size ratio <=1.8;
- minimum angle >=20 deg;
- Level-1 native crack-face sampling not lower than Level 0.

These are engineering/qualification gates, not estimates of physical KII error.

## Prescribed-field controls

Only after all structural gates pass, the existing Step62 synthetic qualification is run on Level 1:

- pure Mode I;
- pure Mode II;
- mixed input KI=1, KII=1e-4.

These controls qualify topology/interpolation/extraction behavior. They do **not** establish physical SIF convergence.

## Outputs

The Level-1 builder appends `_L1` to the save prefix. Expected files include:

- `verification/step65_c03_family_L1_candidate_T3.mat`;
- `verification/step65_c03_family_L1_small_data.mat`;
- `verification/step65_c03_family_L1_overview.png`;
- `verification/step65_c03_family_L1_tip.png`.

Step65 also saves:

- `verification/step65_c03_family_qualification_small_data.mat`

containing the Level-0/Level-1 comparison, sampling table, qualification gates, and Level-1 proposal status.

## Physical-solve policy

Even if Step65 passes every gate, it performs **zero physical FEM solves**.

A successful result sets only `readyForOneLevel1FEMProposal=true`. A Level-1 physical solve would require a separate explicit investigator authorization and a new guarded step.

## Local run

On branch:

`audit/step65-level1-mesh-qualification`

run:

    addpath(genpath(pwd));

    R65 = main_step65_level1_mesh_qualification();

    disp(R65.FamilyComparison);
    disp(R65.SamplingComparison);
    disp(R65.gates);

    if R65.level1.synthetic.performed
        disp(R65.level1.synthetic.table);
    end

Return the complete console output and visually inspect the generated Level-1 overview/tip figures before any physical Level-1 solve is considered.


## Completed local Step65 result — 2026-10-04

The Level-1 C03 mesh-family member passed every predeclared qualification gate with zero physical FEM solves.

### Level-0 reproduction

The family generator first reproduced the accepted Level-0 C03 baseline:

- T3 triangles: 32,980;
- T6 nodes: 66,854;
- paired radius: 6 mm;
- target/actual tip median: 0.054025 mm;
- maximum adjacent-size ratio: 1.7968.

### Level-1 geometry and resolution

Level 1 was generated with scale 0.5:

- T3 triangles: 122,690;
- T6 nodes: 246,700;
- paired radius: 6 mm;
- target tip median: 0.027012 mm;
- actual tip median: 0.027012 mm;
- paired-patch minimum angle: 40.815 deg;
- exterior minimum angle: 25.340 deg;
- maximum adjacent-size ratio: 1.7994;
- primary q-support elements: 40,146;
- literal q-support elements: 44,130.

Thus the structured tip scale was halved exactly while the same 6-mm paired region and physical geometry were preserved.

### Native crack-face sampling

The native points increased approximately by a factor of two:

- [0.04,0.20] a0: 38 -> 74;
- [0.04,0.30] a0: 55 -> 108;
- [0.08,0.30] a0: 44 -> 86;
- [0.12,0.30] a0: 34 -> 67.

### Structural gates

All Step65 family gates passed:

- Level-0 baseline reproduced;
- Level-1 structural pass;
- physical geometry/boundary/material unchanged;
- complete T3/T6 reflection pairing;
- q-support entirely inside paired region;
- exterior excluded from q-support;
- six-triangle/seven-edge tip fan;
- target and measured tip scales halved;
- same 6-mm paired radius;
- maximum adjacent-size ratio <=1.8;
- minimum angle >=20 deg;
- native sampling not reduced;
- synthetic qualification performed and passed;
- zero physical FEM solves.

### Prescribed-field qualification

Level-1 prescribed controls:

- pure I: recovered KI=1.00000000325, KII=1.5167e-14;
- pure II: recovered KI=7.7047e-16, KII=1.0000000024;
- mixed KI=1, KII=1e-4: recovered KII=1.00000000255e-4.

The relative tiny-KII recovery error is approximately 2.55e-9, substantially smaller than on Level 0.

### Interpretation

Step65 establishes a clean, deterministic Level-1 member of the same structured C03 mesh family. It preserves the physical problem and extraction support while halving the local structured scale. This creates the first meaningful direct mesh-family convergence test for the tiny physical Mode-II signal.

The Level-1 candidate is qualified for a future one-solve physical proposal, but no physical Level-1 FEM solve has yet been authorized or performed.
