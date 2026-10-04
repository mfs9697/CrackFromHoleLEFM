# Step 69: additional s=1/sqrt(2) physical convergence point

## Purpose

Step69 inserts one additional mesh scale between the already solved physical scales:

- coarse: s0 = 1;
- midpoint: sm = 1/sqrt(2) = 0.7071067811865476;
- fine: s1 = 1/2.

The refinement ratio is constant:

`r = s0/sm = sm/s1 = sqrt(2)`.

This is the most efficient additional physical mesh for estimating an observed convergence order without changing the mesh design.

## Continuous-scale extension

The structured patch builder previously exposed only integer levels with `scale=2^(-level)`.

Step69 adds an optional explicit continuous `Scale` to the shared mesh-only builder. Existing Level-based calls are unchanged when `Scale` is omitted.

The size law remains:

`h(r) = scale * (hTip + 0.028*r)`.

Thus the midpoint target tip size is approximately:

`0.0540246508 mm / sqrt(2) = 0.0382011969 mm`.

All other design choices remain fixed:

- paired radius 6 mm;
- exterior transition length 8 mm;
- far slope 0.10;
- boundary-metric growth 0.25;
- neighbor ratio target 1.8;
- same deterministic radial-zipper topology;
- same physical boundary geometry.

## Mesh qualification before physical solve

Before any physical solve, the midpoint mesh must pass:

- all structural gates;
- complete T3/T6 reflection pairing;
- q-support entirely inside the paired region;
- exterior excluded from q-support;
- six-triangle/seven-edge tip fan;
- maximum adjacent-size ratio <=1.8;
- same prescribed pure-I, pure-II, and KI=1/KII=1e-4 synthetic controls.

If any mesh or synthetic gate fails, the physical solve is not reached.

## Physical solve

`AllowSolve=false` by default.

With explicit authorization Step69 performs exactly one midpoint physical solve using the already qualified solver:

- free-DOF SPD formulation;
- symamd ordering;
- parameter-free SGS preconditioner;
- PCG tolerance 1e-10;
- maximum 5000 iterations;
- no direct backslash;
- no solver tuning.

The midpoint field is checkpointed before COD or EDI.

## Postprocessing

Step69 uses exactly:

- the same four COD windows and degrees 1/2;
- one 16-point FE-nodal-q interaction EDI;
- r_inner=0.8 mm;
- r_outer=5.2 mm;
- no radius sweep.

## Three-level convergence

For each COD ratio and for the matched EDI ratio, the three values are arranged as:

`q0 = q(s=1), qm = q(s=1/sqrt(2)), q1 = q(s=1/2)`.

When consecutive changes have the same sign and decrease consistently, the observed order is:

`p = log((q0-qm)/(qm-q1)) / log(sqrt(2))`.

The Richardson estimate is then:

`q_inf = q1 + (q1-qm)/(sqrt(2)^p - 1)`.

No order is reported when the three points do not satisfy the monotone/decreasing-difference condition.

## Endpoint references

Step69 prefers exact saved Step67A and Step68 small-data files when they survive locally.

If they are unavailable after branch switching:

- COD uses the high-precision Level-0 fingerprint and the Level-1 values reconstructed from the recorded Step68 Level-1-minus-Level-0 deltas;
- EDI uses the displayed/rounded audit fingerprints.

Therefore the driver explicitly marks EDI `p` and Richardson values as approximate unless both exact endpoint files are available.

## Local run

On branch:

`audit/step69-sqrt2-scale-convergence`

run:

    addpath(genpath(pwd));

    R69 = main_step69_sqrt2_scale_convergence( ...
        'AllowSolve', true);

    disp(R69.midMeshSummary);
    disp(R69.solverInfo);
    disp(R69.ThreeScaleCOD);
    disp(R69.ThreeScaleEDI);
    disp(R69.gates);

If the midpoint mesh fails qualification or PCG fails with the fixed settings, do not alter parameters and rerun. Return the complete output.
