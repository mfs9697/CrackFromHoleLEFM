# Evidence and provenance

Prepared 2026-10-07. Scientific code/plot snapshot:
`cfdf0f110010d688a4e3c48f6d88a00fd17dc698`. Authoritative numerical input:
`verification/crack_path/final_clean_run`. These investigator-local MAT
files are untracked in the starting checkout; the Git commit alone does not
identify the data. `data/source_manifest.csv` records source bytes and SHA256.
`data/evidence_exact.mat` retains the exact extracted double-precision data;
CSV is a convenient readable export.

## Audit performed for the manuscript

`export_pilot_evidence.m` loads the atomic State and compact records, without
calling any mesh generator, EDI replay, or physical solver. It checks all
22 physical files P2–P23, all 23 qualifications P2–P24, all physical and
synthetic gates, paths, solver tolerances, native sampling, and support.
It also recomputes MTS from the stored SIFs as a convention check.
SIFs, ratios, turns, and next angles agree exactly with the State rows;
current angles agree within 4.352074256530614e-14 degrees. This is semantic
agreement, not a claim that every stored angle is bitwise equal.

The original Stage-I file at the standard verification path is absent.
The exact R0 saved by the investigator for the preceding profiling task was
available and is now copied byte-for-byte to `data/accepted_stage1_source.mat`.
R0.C equals the clean CP2 C exactly; its mouth and normal agree with the
first stored segment to 2e-12. No Stage-I solve or nominal reconstruction was
performed. The original P1 field is not present in this clean archive;
P1 is its accepted seed record, not an additional newly solved clean-run file.

## Important statement-to-source map

Paths in the table are repository-relative. `R` and `Small` are the variables
inside their respective compact MAT files. Pk uses `step_kkk_*`.

| Manuscript statement | Exact source/field or derivation |
|---|---|
| A=0.300 m, B=0.100 m; hole center (0.170,-0.020) m, radius 0.030 m, 480 segments | Clean CP2 `C.A`, `C.B`, `C.hole`; restored `R0.summary.{A_m,B_m,hole_x_m,hole_y_m,hole_R_m,hole_npoly}` |
| Plane strain, E=210000 MPa, nu=0.30; unit y tension; minimal corner anchors | CP2 `mat.{E,nu,ps,D}` and `C.{load,bc}`; `solve_incremental_crack_tip.m` loading/anchor construction |
| Continuum/plane-strain equations | Algebraic form of the D matrix in `solve_incremental_crack_tip.m` and equilibrium/strain statements in `Documents/CrPathLEFM.tex`; no literature attribution invented |
| Strength 300 MPa and fitted peak 3.714070740812109 MPa; lambda_ini=80.77390575883399 | `R0.C.sig_c`, `R0.summary.{sigma_tt_peak_unit,lambda_ini}`; freeze implementation in `verification/crack_path/main_stage1_freeze_starting_state.m` |
| phi=-1.56061277813902 deg, exact mouth and frame | `R0.summary.{phi_star_deg,x_star_m,y_star_m,nmat_x,nmat_y,that_x,that_y}`; `State.vertices(1,:)` cross-check |
| First segment theta1=0 and 4-mm length; no production sweep | `run_incremental_crack_path.m` initializer/recurrence; `State.thetaDeg(1)`, first two vertices; uniform lengths checked over all 24 geometrical segments |
| P1 KI=0.366479612185, KII=4.19826648316e-6, theta2 about -0.001312722142 deg | `State.rowsThroughCompleted(1,:)`, using `State.rowVariableNames`; accepted seed constants and MTS call in the driver |
| Exact MTS S/T formulas, interval +/-80 deg, +/-10 deg bracket, -2KII/KI guess, zero threshold 1e-12, nearly-pure-II guess, compressive safeguard | `kink_angle_LEFM_MTS.m`; the manuscript reproduces the implemented implicit root rule, not an alternative closed-form expression |
| Last-segment frame and COD apparent-SIF formula | `native_COD_polyline_audit.m`, `verification/sif_audit/native_COD_audit.m`; current-tip frame in `SIF_LEFM_interaction_EDI.m` |
| EDI density, E'/(2Kaux) normalization, FE-nodal q, auxiliary finite-difference h=max(1e-7 m,1e-5 r), 16 points | `SIF_LEFM_interaction_EDI.m`; qualifier/solver call options; CP2 `mat.ps=1` |
| Mesh scales, hTip=0.0270123254 mm, ri=0.4, ro=2.6, core=3, transition=4, cap=2.5 mm | All `Small.summary.{increment_m,hTip_m,rInner_m,rOuter_m,rCore_m,transition_m,farCap_m}`; fixed-ratio assertions in the qualifier |
| Core 12678 T3, EDI support 11316 | Every `Small.summary.{core_T3_elements,EDI_elements}` and physical `R.EDI.EDI_elements` |
| 20-deg minimum-angle/1.8 grading gates and positive elements/topology | `qualify_incremental_crack_candidate.m` gates and every `Small.gates`; actual summaries are exported in `qualification.csv` |
| Native sampling 38/55/44/34 and eight fits | `Small.sampleCounts.nativePoints`; `R.fitTable.n_native`, bounds/degree columns; zero-field/source topology gates |
| Pure-I/II/tiny-mixed controls, gates 2e-4 and 1e-10; maximum tiny-mixed relative error 3.861737685184608e-8 | Every `Small.synthetic`, `Small.syntheticGates`; qualifier gate definitions; aggregated `E.audit.syntheticErrors` |
| PCG/SGS/symamd/free-DOF architecture, 1e-10/5000, true gate 5e-10, exact constraints, checkpoint before postprocessing | Solver implementation and every `R.solverInfo`/`R.gates` |
| Accepted iterations 2442–3011; max true residual 9.996283318623952e-11 | State accepted rows and all `R.solverInfo`; audited aggregates `E.audit` |
| Previously accepted P2 regression | Stored `State.regression` (copied to `E.regression`) and coded bridge in `run_incremental_crack_path.m`; this is reproducibility, not discretization-error calibration |
| Historical straight 8-mm verification is separate | `verification/sif_audit/AUDIT_CLOSURE.md`; no historical SIF/mesh-family values are inserted into the trajectory or used to calibrate it |
| 22 accepted compact files, 23 qualifications; P24 unsolved | Directory inventory, `State.completedPhysicalSegments=23`, 25 vertices; `Small` P24 passes; both P24 physical small and solved MAT files absent |
| KI grows monotonically from P1 to P23 | Positive successive differences in State `KI_unit`, cross-checked against physical compact records |
| Positive KII and ratio maxima at P17, a=68 mm; KII=0.002548368739746823, ratio=0.003238930802982266 | Row17 and physical `R.summary`/`R.EDI`; independent max searches of the two stored sequences |
| P21 ratio +8.134267373465366e-5; P22 ratio -0.002141591607115312 | Rows21/22 and corresponding physical compact files |
| Linear q=0 at 84.14636991194091 mm, displayed 84.15 | `-q21/(q22-q21)` between a21=84 and a22=88 mm; only visualization. Position marker interpolates the two tip coordinates, never assigns an elastic solution |
| P23 KI=1.244360443656587, KII=-0.006780312460038702, ratio=-0.005448833169362544 | Row23 and `step_023_physical_small.mat` `R.summary`/`R.EDI` |
| P23 next turn +0.6243470380022705 deg, theta24=-2.776211096541598 deg | Row23 `delta_theta_next_deg`, `theta_next_deg`; corresponding R prediction. This predicts P24 geometry only |
| Absolute direction minimum theta22=-3.645963829383184 deg | Min search of the accepted State `theta_deg`; P22 R current angle agrees within roundoff |
| P23 (291.8404840685052,-25.57894113432961) mm, right-side horizontal ligament 8.159515931494777 mm | Row23 coordinates; `(C.A-tip_x_m)*1000` |
| Characteristic-state table | Rows1,17,21,22,23; rounded only for presentation. Interpolated row has no KI/KII/direction values |
| 176 COD definitions all agree on positive P21/negative P22 | All eight stored `R.fitTable.ratio_COD` values at each of those two states; `derived_metrics.json` sign checks |
| Largest COD-EDI turn difference 0.01023605925082449 deg | Max absolute `R.fitTable.delta_theta_next_MTS_deg - State.delta_theta_next_deg`; P23, degree1, window [0.12,0.30] |
| P23 quadratic turns 0.6204355782965233–0.6224779975149831 deg | All four degree2 rows of P23 fit table |
| Near-zero P21 rear linear ratio 26.31962452791384% above EDI, absolute turn difference 0.002453300431192729 deg | P21 degree1 [0.12,0.30] row: `100*(ratio_COD/ratio_EDI-1)` and absolute turning difference; this is not a relative physical error estimate |
| Secondary Stage-I tensile maximum about 1.34% lower (remaining-work report) | `R0.summary.primary_secondary_gap_rel=0.01336323222072076` |

## Centered-hole symmetry and local-stability benchmark

A separate physical control was run on 2026-10-07 from
`symmetric-path-stability` commit
`6945371dfdfac441b90970690b7e1728c7f06259`. Its complete MATLAB console
log has SHA256
`f4dd0397a86fb1ad395270a496b1eac209effe450dddc8dc3603ac1852695301`.
The readable numerical rows are preserved in
`data/symmetric_stability_benchmark.csv`; detailed source/run provenance is
in `data/SYMMETRIC_STABILITY_PROVENANCE.md`.

The benchmark uses a centered circular hole, an exact horizontal first
4-mm segment, and prescribed second-segment perturbations of
-0.10, -0.05, 0, +0.05, and +0.10 degrees. Each full-domain geometry was
independently qualified and solved once. All qualification and physical
gates passed. The straight control gives
`KII/KI=6.22351746186701e-09` and
`theta3=-7.1316256897125e-07 deg`. The pairwise one-step multipliers are
`-0.114158317407462` and `-0.114150920881333`, and the four nonzero
probes give the origin-constrained fit
`F'(0)=-0.114152400186559`.

This control is evidence for symmetry preservation, sign conventions, SIF
extraction, and the local EDI/MTS directional update. It is deliberately
kept separate from the asymmetric production trajectory. It is not a
4/2/1-mm trajectory-convergence study and it does not establish that the
free boundary is the sole cause of the late mode-mixity reversal.

## Claim limits and missing provenance

The data support the sign change and rising KI. The nearby-boundary
interpretation is not isolated by a separate control. There is no
increment-size convergence study for this trajectory, and the historical
8-mm mesh study cannot supply one. COD shares the FEM field with EDI and
therefore is an independent extractor, not an independent physical solution.

P24 stagnation just outside the strict tolerance is supplied by the task
description. The archive confirms nonacceptance but has no retained P24
failure trace, final residual, iteration count, or stagnation diagnosis.
The manuscript labels the attempt as reported and supplies no invented
failure metrics. There is no claim of arrest, branching, dynamics, or
boundary intersection.

No complete bibliography entry was found in the inspected source material.
Literature TODOs are explicit. The source/manuscript layers remain separate;
the pilot adds no new numerical validation experiment.
