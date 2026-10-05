# Incremental crack-path profiling: Phase A

Measured on 2026-10-05 with MATLAB R2023a (9.14.0.2206163), Ryzen 7 4800H,
Windows, from scientific base `10c39ac1e38a76520817d44f4a2d31ef28cb25e2`.
The exact accepted, investigator-saved R0 was supplied in memory. No Stage-I
state was reconstructed. Runs used separate empty output directories and the
original accepted five-step results as full-precision reference. The existing
interactive MATLAB session and original checkpoints were preserved.

## Result and stop-condition decision

**PCG is the dominant cost, as anticipated by the supplied 44.6–48.0 s
baseline.** It consumes 190.445 s, or 56.37% of the fresh five-segment run.
No unexpected dominant cost was discovered. The user-requested conditional
stop therefore does not apply. Phase B may investigate measured secondary
costs while retaining the solver and all scientific gates.

| Strict baseline | Wall time (s) | New physical solves | Result |
|---|---:|---:|---|
| Fresh three segments | 169.347027 | 2 | PASS |
| Fresh five segments | 337.828553 | 4 | PASS |
| Resumed five segments | 20.592993 | 0 | PASS |

The maximum full-precision differences in KI, KII, incremental next angle,
absolute next angle and current angle are **exactly zero** for all three runs
relative to the accepted stored reference. Existing coded P2 tolerances remain
5e-10 for SIFs and 5e-9 degrees for angles; the harness applies those same
limits to P3–P5 as additional comparisons without changing production gates.
All structural, synthetic, solver, COD and EDI gates pass. Native samples
remain 38/55/44/34, EDI support 11316, and core triangles 12678.
PCG iterations remain 2442/2457/2482/2509 and true residuals remain below 1e-10.

## Dominant exclusive phases: fresh five segments

| Rank | Phase | Total across P2–P5 (s) |
|---:|---|---:|
| 1 | PCG | 190.445 |
| 2 | Synthetic qualification, principally three EDI replays per tip | 56.161 |
| 3 | Exterior triangulation/refinement | 21.517 |
| 4 | Physical EDI | 18.520 |
| 5 | Unclamped stiffness assembly | 13.674 |
| 6 | Geometry mapping/bookkeeping | 9.841 |
| 7 | Source carrier, including its internal geometry-ID fallback | 8.979 |
| 8 | Structured core | 5.526 |
| 9 | Structural qualification | 3.614 |
| 10 | symamd and sparse permutation | 2.050 |

Whole qualification totals 107.790 s. These inclusive totals must not be
added to their child phases. The CSVs separately record Williams replay and
synthetic EDI; their times are already included in synthetic qualification.
Carrier geometry-only ID identification falls back to a temporary PDE mesh;
that internal work is included in source_carrier. The separately timed
geometry_id_recovery is the subsequent polyline edge-set identification.
T3-to-T6 totals about 0.155 s over all three conversions per tip, SGS 0.188 s,
native COD 0.048 s, and physical checkpoint writes 0.807 s. These are poor
first optimization targets. Resumed runtime is predominantly physical EDI;
qualification, assembly and PCG are skipped by existing validated reuse.

## Memory and environment

MATLAB phase-boundary memory samples peak at 2.304 GB (2.146 GiB) in the
fresh five-segment run; these samples are not continuous peak measurements.
The separately observed Windows process peak working set was 1.449 GB.
MATLAB's memory accounting differs from resident working set. The P5 SGS
stored-factor estimate is 0.0771 GiB. Assembly triplet arrays alone reserve
approximately 170 MiB at P5 (144 entries/element, three double arrays), before
sparse finalization. These figures suggest useful future memory work but do
not establish allocation as the leading runtime cost.

The pre-existing assembly helper resolves outside the repository:
`C:\Users\Mikhailo\Google Drive\Projects\8 - CohZonFEM\czm_project\BN_local.m`,
SHA256 `A52287AB4F3394C01551B044E13F4CADCD7AD51789813774BDF16B3C5DE587CA`.
The repository stif_assem is the active assembler. This dependency was not
modified; another machine must provide the same helper to reproduce assembly.

## Measured opportunities

A separate MATLAB profiler replayed P5 EDI from its saved field without a
physical solve. Its instrumented 13.572 s is **not** a production benchmark
(uninstrumented EDI is about 4.6 s). It identifies 359100 auxiliary-field
calls and 1795500 auxiliary-displacement evaluations. In particular, every
auxiliary evaluation calculates a displacement u0 that is never used, and
computes a 3-by-3 constitutive solve plus mismatch norms even when
StoreGPDiagnostics=false discards all mismatch diagnostics.

A conservative first experiment is an explicit, default-off EDI mode that
omits only this unused work when GP diagnostics are already disabled. It
must preserve every value entering the integral, quadrature, accumulation
order, support discovery and reported diagnostics. Protect it with bitwise
stored-field comparisons and fresh three/five-segment regression runs.
No evolving-geometry cache is needed for this experiment.

PCG remains the longer-term limit. Warm starts, reused permutations or
preconditioners would require separate numerical qualification and are not
justified as the first optimization. Exterior bookkeeping is the next
substantial mesh target; measure individual scans before changing it.

## Reproduction and evidence

```matlab
Report = main_incremental_path_profile('FrozenState',R0, ...
    'AllowPhysicalSolves',true,'OutputDir','<new empty benchmark directory>');
```

See `profiling/baseline/*_tips.csv` for every requested per-tip value,
`*_phases.csv` for all measured phases, `summary.json` for wall times and
full-precision differences, and `profile_small.mat` for exact stored timing
and comparison data. `edi_function_profile.csv` contains the separate
profiler's inclusive function times; never sum inclusive parent/child rows.
Timing instrumentation is disabled by default and enabled only for these
explicit sessions. A nested-clock test passed. Phase-transition memory
sampling introduces observer overhead; parent totals include that overhead.
These are single sequential runs, not a statistical estimate of speedup.
No optimization or before/after speedup is claimed in Phase A.
