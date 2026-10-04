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


## Completed local Step69 result — 2026-10-04

Step69 completed successfully with one explicitly authorized physical solve at

`s = 1/sqrt(2) = 0.707106781187`.

All midpoint mesh, synthetic-field, solver, and postprocessing gates passed.

### Midpoint mesh

The continuous-scale mesh has:

- target tip size: 0.0382011969 mm;
- actual tip median: 0.0382011969 mm;
- T3 elements: 62,570;
- T6 nodes: 126,112;
- paired radius: 6 mm;
- paired-patch minimum angle: 40.7561 deg;
- exterior minimum angle: 25.1372 deg;
- maximum adjacent-size ratio: 1.79593993;
- primary q-support elements: 20,058;
- literal q-support elements: 22,134.

Native COD sampling:

- [0.04,0.20] a0: 53 points;
- [0.04,0.30] a0: 77 points;
- [0.08,0.30] a0: 61 points;
- [0.12,0.30] a0: 47 points.

### Prescribed-field qualification

The midpoint mesh passed all three prescribed controls:

- pure I: recovered KI=0.999999984677, KII=8.9265e-15;
- pure II: recovered KI=-1.3404e-14, KII=0.999999994039;
- mixed KI=1, KII=1e-4: recovered KII=9.99999994128e-5.

### Midpoint SGS-PCG solve

- PCG tolerance: 1e-10;
- maximum iterations: 5000;
- converged flag: 0;
- iterations: 2626;
- reported relative residual: 9.7535e-11;
- true relative residual: 9.7535e-11;
- exact homogeneous constraints;
- stiffness symmetry error: 1.6366e-16;
- SGS approximate sparse storage: 0.0935 GiB;
- solve time: 58.2396 s;
- available physical memory before assembly: 5.9739 GiB;
- available after preconditioner construction: 5.9143 GiB.

The physical midpoint field was checkpointed before COD/EDI postprocessing.

### Three-scale COD convergence

The three scales are:

- s=1;
- s=1/sqrt(2);
- s=1/2.

All eight COD estimators are monotone and have decreasing successive changes.

Observed orders and Richardson estimates reported by Step69:

- [0.04,0.20] a0, degree 1: p=1.8017, q_inf=1.0588e-4;
- [0.04,0.20] a0, degree 2: p=1.3677, q_inf=1.0652e-4;
- [0.04,0.30] a0, degree 1: p=2.1021, q_inf=1.0541e-4;
- [0.04,0.30] a0, degree 2: p=1.4616, q_inf=1.0648e-4;
- [0.08,0.30] a0, degree 1: p=2.0102, q_inf=1.0502e-4;
- [0.08,0.30] a0, degree 2: p=1.6616, q_inf=1.0641e-4;
- [0.12,0.30] a0, degree 1: p=1.7224, q_inf=1.0463e-4;
- [0.12,0.30] a0, degree 2: p=1.7087, q_inf=1.0636e-4.

Thus every COD extraction window/degree combination is admissible under the predeclared monotonicity/decreasing-difference rule.

The quadratic COD extrapolations cluster especially tightly around approximately 1.064e-4 to 1.065e-4.

### Three-scale matched EDI convergence

Using the same fixed 0.8-5.2 mm interaction domain:

- coarse s=1: KII/KI = 1.0631e-4;
- midpoint s=1/sqrt(2): KII/KI = 1.0646e-4;
- fine s=1/2: KII/KI = 1.0652e-4.

Successive relative changes:

- coarse -> midpoint: +0.13783%;
- midpoint -> fine: +0.059626%.

The sequence is monotone and the increment decreases by more than a factor of two.

Step69 reports:

- approximate observed order p = 2.4137;
- approximate Richardson limit KII/KI = 1.0657e-4.

The endpoint-reference flag is false because the exact saved Step67A/Step68 compact files were not available after branch switching. Therefore these EDI order/extrapolation values use displayed/rounded endpoint audit fingerprints and must be described as approximate.

### Interpretation

The additional geometrically spaced midpoint confirms that the tiny positive Mode-II signal is not only stable under one h/2 refinement, but follows a smooth monotone three-scale trend.

For EDI, the refinement increments decrease from approximately +0.138% to +0.060%. For all eight COD fits, both refinement increments have the same sign and the second is smaller than the first.

This is strong evidence that the sequence has entered a convergent regime.

Because the endpoint audit fingerprints are rounded/reconstructed rather than exact surviving machine-precision small-data files, the reported observed orders and Richardson values should be treated as approximate numerical indicators rather than formal high-precision estimates.

A defensible summary value from the matched EDI sequence is:

`KII/KI approximately 1.066e-4`

with an approximate three-level extrapolated value near

`1.0657e-4`.
