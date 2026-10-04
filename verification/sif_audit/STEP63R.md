# Step 63R: recover the lost Step63 physical checkpoint

## Why this step exists

The first Step63 physical solve completed successfully and saved a checkpoint, but that local untracked MAT file was later lost. The first-run console log preserved the solved-field COD fingerprint.

Step63R is therefore a **replacement solve**, not a new mesh experiment and not a new physical configuration.

## Guard

`AllowRecoverySolve=false` by default.

Only an explicitly authorized run with:

    R63r = main_step63r_recover_lost_physical_field( ...
        'AllowRecoverySolve', true);

may recreate the lost field.

## Exact problem reused

Step63R calls the already audited Step63 driver. That driver:

- deterministically reconstructs C03 if its generated MAT is absent;
- requires 32,980 T3 triangles and 66,854 T6 nodes;
- requires the calibrated ratio 1.79678451 within tolerance;
- requires the six-triangle / seven-edge tip fan;
- reruns and passes the prescribed pure-I, pure-II and KII=1e-4 controls;
- uses unit remote-y traction and minimal anchoring;
- verifies the solver used the exact candidate T3/T6 mesh;
- saves the replacement physical field before COD postprocessing.

Step63R itself contains no direct solver call. Across Step63R plus the invoked Step63 implementation there is exactly one `solve_cracked_LEFM` call site.

## First-run COD fingerprint

The recorded first successful Step63 run reported the eight COD-fit ratios:

- 1.0526e-4
- 1.0562e-4
- 1.0489e-4
- 1.0571e-4
- 1.0456e-4
- 1.0582e-4
- 1.0416e-4
- 1.0583e-4

and raw-band median ratios:

- 1.0365e-4
- 1.0178e-4
- 9.9274e-5
- 9.5709e-5
- 9.0841e-5.

These values came from displayed console output and are rounded. Step63R therefore uses an absolute comparison tolerance of 5e-8 on the dimensionless ratio, approximately 0.05% of the physical signal.

The replacement field must also preserve:

- the same COD windows and native sample counts;
- positive signed Mode-II signal in every reported fit/band;
- complete fit-ratio range inside 1.03e-4 to 1.07e-4.

If any fingerprint gate fails, Step63R stops and explicitly forbids proceeding to Step64.

## No EDI

Step63R performs no physical EDI.

After a successful recovery, **do not switch branches**. The Step64 driver is permitted on the same Step63R branch so the newly restored checkpoint can be used immediately without another artifact-loss opportunity.

## Run

On branch:

`audit/step63r-recover-lost-physical-field`

run only after explicit authorization:

    addpath(genpath(pwd));

    R63r = main_step63r_recover_lost_physical_field( ...
        'AllowRecoverySolve', true);

    disp(R63r.fitFingerprintTable);
    disp(R63r.rawFingerprintTable);
    disp(R63r.gates);

If and only if Step63R prints `STEP63R PASS`, the already-authorized Step64 EDI can then be run on the **same checked-out branch**:

    R64 = main_step64_matched_physical_edi();

Do not switch branches between recovery and Step64.
