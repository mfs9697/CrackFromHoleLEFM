# Incremental crack-path optimization report

Measured 2026-10-05, MATLAB R2023a on Ryzen 7 4800H/Windows.
Scientific base: `10c39ac1e38a76520817d44f4a2d31ef28cb25e2`.
Implementation: `91ea3001dead93a94c87fc57b49614e8178fd4e5`, on
`optimize-incremental-crack-path`. The validated branch and main were not
modified. Independent Stage III-C/D driver files were not changed.

## Outcome

The explicit, default-off **FastEDI** mode improves measured runtime while
preserving all scientific results. All baseline and optimized fresh three-
and five-segment runs pass; both resumed five-segment runs perform zero new
solves. **Every P2–P5 mesh, material, crack topology, solved displacement
vector, structural/synthetic result, KI/KII and trajectory angle is bitwise
identical between the baseline and optimized runs.**

Phase A found the expected PCG bottleneck, 56.37% of baseline fresh-five wall
time. No unexpected dominant cost was found, so the user's conditional stop
did not apply. See [the Phase A report](INCREMENTAL_PROFILE_PHASE_A.md).

| Run | Strict baseline (s) | FastEDI (s) | Observed speedup | Time reduction |
|---|---:|---:|---:|---:|
| Fresh three segments | 169.347027 | 156.755652 | 1.080× | 7.44% |
| Fresh five segments | 337.828553 | 306.202659 | 1.103× | 9.36% |
| Resumed five segments | 20.592993 | 15.002831 | 1.373× | 27.15% |

These are single sequential end-to-end runs, with the same phase clocks and
strict gates. PCG timing varies despite unchanged assembly inputs, identical U and iteration
counts: 190.445 s before versus 182.095 s after across P2–P5. That variation
is not a solver optimization. The EDI-specific improvement is independently
supported by three warmed, alternating-order paired measurements: median
4.521341 s strict versus 3.256949 s fast, **1.388214×**. Do not interpret
the full end-to-end difference as a guaranteed speedup on another machine.

## One conservative optimization

Every auxiliary-field evaluation previously computed:

- u0 at the unshifted point, whose value is never used;
- Dmat\\sig, a strain used only for the auxiliary mismatch diagnostic;
- mismatch denominator and mismatch norms, even with StoreGPDiagnostics=false.

FastEDI omits precisely these discarded calculations. Stress, finite-
difference displacement derivatives in both directions, strains, q gradients,
the interaction formula, all 16 quadrature points, element/GP traversal and
summation order are unchanged. All three synthetic EDI controls still run.
No verification gate or tolerance is weakened.

The extractor option `SkipUnusedAuxWork` defaults to false. Even if requested,
it takes effect only when `StoreGPDiagnostics=false`; diagnostics-enabled
calls retain the original work and outputs. `FastEDI=false` is the default
in the driver, qualifier, physical solver, regression and benchmark harness.
R, F and Path record which mode was requested. Candidate and field reuse
remain the existing validated behavior. No quantity is cached across changing
paths or meshes, so no new cache validity assumption is introduced.

## Measured exclusive phases: fresh five segments

| Phase across P2–P5 | Baseline (s) | FastEDI (s) |
|---|---:|---:|
| PCG | 190.445 | 182.095 |
| Synthetic qualification | 56.158 | 40.323 |
| Exterior | 21.522 | 20.957 |
| Physical EDI | 18.519 | 13.003 |
| Assembly | 13.674 | 13.697 |
| Geometry mapping | 9.838 | 9.647 |
| Carrier | 8.976 | 8.851 |

Whole qualification decreases from 107.790 to 90.491 s. Its child phases are
already included in that total; do not add them again. About 21.350 s is
removed from synthetic qualification plus physical EDI. Other phase timing
differences are observational variation. PCG still consumes about 59.5% of
optimized fresh-five runtime and remains the primary bottleneck.

## Per-tip timing tables

P1 is the accepted seed, not a newly solved or benchmarked state. Each fresh
three-segment run solves P2/P3; each fresh-five run solves P2–P5 once.
The tables below use the five-segment runs. Separate three-segment and resumed
tables, including all requested fields, are supplied in the CSV evidence.

Strict baseline:

| Tip | T3 elements | T6 nodes | Qualification s | Assembly s | PCG iterations | PCG s | Postprocessing s | Total step s |
|---|---:|---:|---:|---:|---:|---:|---:|---:|
| P2 | 49637 | 100287 | 27.386 | 3.537 | 2442 | 46.524 | 4.700 | 83.859 |
| P3 | 50518 | 102122 | 26.577 | 3.353 | 2457 | 47.274 | 4.743 | 83.668 |
| P4 | 51162 | 103478 | 26.892 | 3.353 | 2482 | 48.384 | 4.781 | 85.135 |
| P5 | 51512 | 104242 | 26.936 | 3.431 | 2509 | 48.263 | 4.654 | 84.964 |

FastEDI:

| Tip | T3 elements | T6 nodes | Qualification s | Assembly s | PCG iterations | PCG s | Postprocessing s | Total step s |
|---|---:|---:|---:|---:|---:|---:|---:|---:|
| P2 | 49637 | 100287 | 23.298 | 3.702 | 2442 | 45.066 | 3.330 | 77.132 |
| P3 | 50518 | 102122 | 21.971 | 3.268 | 2457 | 44.694 | 3.311 | 74.793 |
| P4 | 51162 | 103478 | 22.682 | 3.339 | 2482 | 45.916 | 3.332 | 76.851 |
| P5 | 51512 | 104242 | 22.540 | 3.388 | 2509 | 46.418 | 3.333 | 77.241 |

Shared scientific results (SIF units MPa sqrt(m), angles in degrees):

| Tip | KI | KII | KII/KI | Current angle | Incremental next turn | Absolute next angle |
|---|---:|---:|---:|---:|---:|---:|
| P2 | 0.43786411992502 | 4.1052401599805e-05 | 9.37560300826539e-05 | -0.00131272214159878 | -0.0107436494349182 | -0.012056371576517 |
| P3 | 0.474334420575756 | 0.000109830836465289 | 0.000231547262228987 | -0.0120563715765322 | -0.0265333584477464 | -0.0385897300242786 |
| P4 | 0.499809523097675 | 0.000188172855946581 | 0.00037648913686226 | -0.0385897300242805 | -0.0431424628806929 | -0.0817321929049734 |
| P5 | 0.521402328367184 | 0.000269154507476197 | 0.000516212707218008 | -0.0817321929049561 | -0.0591535821289529 | -0.140885775033909 |

CSV decimal serialization is for inspection. Exact binary doubles and
full-precision comparisons are retained in the compact MAT files.

## Validation

- The existing main_incremental_path_regression passes in the baseline and
  optimized fresh runs, and again in post-change default strict mode using
  its matching saved fields.
- All structural, pure-I, pure-II, tiny-mixed, solver, COD and EDI gates pass.
  The current-tip core remains 12678 triangles; primary EDI support remains
  11316 elements and native sampling remains 38/55/44/34.
- All five tracked numeric comparison maxima are exactly zero relative to
  the accepted stored five-step reference: KI, KII, current angle,
  incremental next turn and absolute next angle. Existing coded SIF/angle
  tolerances remain 5e-10 and 5e-9 degrees.
- All four P2–P5 meshes, materials, crack structures and full U vectors
  compare bitwise equal, as do structural gates, synthetic tables, synthetic
  gates and sample-count tables. True residuals remain below 1e-10 and
  below the unchanged 5e-10 acceptance gate.
- 24 standalone EDI cases compare all returned SIFs and Aux diagnostics
  bitwise, covering 7/12/16 quadrature, analytic/FE-nodal q, diagnostic
  enablement and exact/FE actual-field paths. Stored P5 was also compared
  bitwise in warmup and every paired timing sample.
- Optional clock inactive/nested-scope tests pass. MATLAB Code Analyzer
  reports only existing style/performance advisories, no syntax failures.

The 4-mm increment, theta1=0, MTS signs and absolute-angle recurrence, true
polyline topology and last-segment frame, mesh scales, COD extractor, EDI
formula/quadrature, free-DOF SPD assembly, symamd, SGS, PCG 1e-10/5000,
exact constrained DOFs and checkpoint-before-postprocessing order are retained.
No angle sweep, parallelization, new Stage-I solve or global direct solve was
introduced.

## Memory, remaining costs and long-path recommendation

Phase-boundary MATLAB memory samples peak at 2.304 GB before and 2.310 GB
after; this optimization does not demonstrate a memory reduction. Assembly
triplets and sparse solver storage are unchanged. Windows observed peak
working sets were about 1.449 and 1.451 GB, respectively. See Phase A for
sampling limitations, SGS storage and the existing external BN_local helper
dependency and its SHA256. That helper was preserved throughout the tests.

For a later long-path run, preserve and explicitly load the accepted R0,
keep strict gates enabled, and choose a separate durable output directory.
FastEDI=true is now qualified for this five-step case; record that choice.
Retain the existing physical checkpoints and use validated resume behavior
after interruption. Monitor per-tip wall times, PCG iterations/true residual,
all support/sampling fingerprints and the existing core-clearance termination.

A future solver experiment needs independent numerical qualification: do not
infer that warm starts, permutation reuse or preconditioner reuse are safe
from this EDI-only result. Exterior mapping/full-mesh scans are the next
measured mesh candidates. T3-to-T6, native COD and SGS construction are too
small to prioritize. Qualification frequency and mesh density stay unchanged.

## Reproduction and evidence

```matlab
% Default strict benchmark; use a new empty output directory.
Baseline = main_incremental_path_profile('FrozenState',R0, ...
    'AllowPhysicalSolves',true,'OutputDir','<new strict benchmark directory>');

% Explicit optional mode; use another new empty output directory.
Fast = main_incremental_path_profile('FrozenState',R0, ...
    'AllowPhysicalSolves',true,'FastEDI',true, ...
    'OutputDir','<new fast benchmark directory>');

test_incremental_fast_edi('CheckpointFile','<accepted P5 solved MAT>');
```

[Baseline evidence](profiling/baseline/summary.json) and
[optimized evidence](profiling/fast/summary.json) provide wall times and
comparison maxima. Each folder contains fresh3/fresh5/resumed5 phase and
per-tip CSVs and profile_small.mat. Optimized evidence also includes the
paired EDI probe and final bitwise-verification JSON/MAT. Solved field
checkpoints and R0 remain investigator-local and are not published with
these compact artifacts. Report.baselineCommit records the checkout HEAD
at measurement; scientificBaseCommit identifies the validated source.
The optimization was tested before committing it as 91ea300.
