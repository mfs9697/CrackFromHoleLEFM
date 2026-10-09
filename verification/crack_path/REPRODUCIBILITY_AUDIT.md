# Reproducibility forensic audit

Date: 2026-10-09. Branch: `tip-refinement-fixed-geometry`.
Starting scientific/code snapshot: `79c2a9edb3542d5bbe22eba6cb82e704e75aa22a`.

**The audit found a critical unresolved experimental-isolation issue:
CoreScale changes the exterior size function as well as the tip core.**
Four other critical reuse/publication defects were mechanically corrected.
The mesh laws, EDI/COD/MTS formulas, scientific interpretation, numerical
acceptance criteria, and their tolerance values were not changed. No physical
solve, PDE mesh generation, or physical-field EDI replay was run.

The existing dirty manuscript PDFs/PNGs and investigator numerical archives
were preserved and are not part of the audit commits.

## Scope and evidence

An exhaustive pattern search covered **297 working-tree text files / 72,211
lines**, including 219 MATLAB files, 3 Python files, 10 TeX files, 55 Markdown
files, a bibliography, and 9 JSON snapshots. No shell, batch, PowerShell, or
JavaScript files of the searched types existed in this checkout. The search
included untracked/ignored supported source files and excluded Git internals,
build directories, and Python bytecode. The immutable/generated JSON records
are classified as data snapshots rather than parameter logic.

The scanner recorded **2,262 matched lines** across eight assumption classes.
This is an exhaustive textual search, not a claim that every match is a bug
or that static inspection proves all dynamic MATLAB/path-dependent behavior.
Active parameter/cache/publication paths and the relevant historical/snapshot
callers were manually traced. The raw hits give exact locations for the
remaining literal occurrences, including intentional reference assertions.

Evidence is in [reproducibility_audit](reproducibility_audit/):

- `source_inventory_before.csv`, `assumption_hits_before.csv`,
  `scan_summary_before.json`: original inventory, SHA256, and every hit.
- `stage1_reuse_sites.csv`: all source-code R0/frozen-state reuse occurrences.
- `findings.json` / `findings.csv`: machine-readable finding table.
- `baseline_checks.log`: reproduced failures without a solve.
- `checks.json` / `checks.mat`: lightweight regression results and scale table.
- `exterior_scale_probe.json`: analytical evaluation of the unchanged mesh law.
- `scan_repository.py`: rerunnable scanner; output location is explicit.

## Concise finding table

Severity is one of **critical**, **important**, **cosmetic**, **confirmed safe**.
“Fixed” means a mechanical identity/dispatch/input bug was repaired; “open”
means it was documented and left unchanged because a policy, mesh-law,
interpretation, or archive-provenance decision is needed. References below
use the final audited source locations; the before-search preserves original
line references as well.

| ID | Severity | Status | Finding and disposition | Exact file/line evidence |
|---|---|---|---|---|
| F01 | critical | open | CoreScale is also an exterior-sizing parameter; tip-only attribution is not isolated. No mesh-law or interpretation change permitted; see the nonphysical scale probe. | `build_stage3c_polyline_exterior.m:77`; `build_stage3c_polyline_exterior.m:308`; `build_stage2_scaled_audited_exterior.m:39`; `build_stage2_scaled_audited_exterior.m:224` |
| F02 | critical | fixed | P1 compact reuse accepted a 2-mm seed for a 4-mm request because it checked only pass and scale=2. Validate absolute core/radius fingerprints, exterior controls, source mesh, increment, and frozen physical data before reuse. | `verification/crack_path/main_increment_sensitivity_coarse.m:206`; `verification/crack_path/validate_increment_study_p1_cache.m:8` |
| F03 | critical | fixed | Full field checkpoints could reuse another material/load/physical configuration or changed T6 ordering/coordinates while T3 still matched. Compare saved material/C and full T6 geometry with the same existing coordinate tolerance. | `solve_incremental_crack_tip.m:823`; `main_stage2_theta0_physical_solve.m:766`; `assert_crack_checkpoint_physics.m:1` |
| F04 | critical | fixed | Compact accepted-result promotion checked path and 11316 elements, but not its exterior/physical source identity. New result metadata records exterior/physics; legacy compacts require source candidate/field evidence. | `run_incremental_crack_path.m:805`; `run_incremental_crack_path.m:847`; `solve_incremental_crack_tip.m:660` |
| F05 | important | fixed | CoreScale=2 compact promotion/hydration was rejected by the unconditional reference EDI count 11316. Dispatch the existing fingerprint: 11316 for scale 1, 2976 for scale 2; do not permit arbitrary counts. | `run_incremental_crack_path.m:818`; `run_incremental_crack_path.m:827` |
| F06 | important | fixed | State row names and declared schema were stored but not validated before positional history parsing. Reject unsupported declared schemas and reordered/renamed columns before geometry/solve work. | `run_incremental_crack_path.m:560`; `run_incremental_crack_path.m:648` |
| F07 | important | fixed | NaN family numbers could compare as equal because abs(NaN-x)>tol is false. Reject nonfinite scalar family/core identity values; tolerance values unchanged. | `run_incremental_crack_path.m:1018`; `run_incremental_crack_path.m:1030` |
| F08 | critical | fixed | The 4/2/1 publication loader trusted CSV names/lengths without checking the advertised M1+2h0 family or its authoritative summary. Require study/family provenance, unambiguous grids, vertex agreement and summary concordance. | `verification/crack_path/plot_increment_sensitivity_coarse_publication.m:431`; `verification/crack_path/load_increment_comparison_run.m:6` |
| F09 | important | fixed | M1/reference and 2h0/reference comparison paired equal segment indices without confirming equal increments. Reject mismatched/nonuniform increments using the existing driver geometry tolerance. | `verification/crack_path/plot_m1_vs_reference_publication.m:108`; `verification/crack_path/plot_tip2h0_vs_reference_publication.m:98` |
| F10 | important | fixed | Missing explicitly requested M1 history could silently return canonical committed panels instead of the requested experiment. Canonical fallback retained only for canonical default requests. | `verification/crack_path/plot_m1_vs_reference_publication.m:74`; `verification/crack_path/plot_m1_vs_reference_publication.m:370` |
| F11 | important | fixed | An integer target shorter than two increments could reach P1 work before MaxSegments rejected it. Reject the unsupported target before reading/building/saving/solving. | `verification/crack_path/main_increment_sensitivity_coarse.m:68`; `verification/crack_path/main_increment_sensitivity_coarse.m:74` |
| F12 | important | open | Older fully embedded State histories lack a frozen physical signature; their rows remain a trusted atomic archive, not independently requalified evidence. New states are keyed; do not silently manufacture missing historical provenance. | `run_incremental_crack_path.m:575`; `run_incremental_crack_path.m:648`; `run_incremental_crack_path.m:908` |
| F13 | important | open | Metadata-free historical candidates are inferred to be reference exterior/core=1 when reference controls are requested. Legacy compatibility assumption remains explicit. Known source archives are required for legacy compact promotion. | `run_incremental_crack_path.m:959`; `run_incremental_crack_path.m:982` |
| F14 | important | open | Shared output names can overwrite candidate/qualification files before an incompatible field checkpoint is rejected. Use one isolated directory per increment/core/exterior/frozen source. This audit does not redesign archive transactions. | `run_incremental_crack_path.m:99`; `run_incremental_crack_path.m:238`; `verification/crack_path/main_increment_sensitivity_coarse.m:139` |
| F15 | important | open | Direct generic runs keep the accepted reference 4-mm P1 SIF defaults even if increment/core controls are changed. Use an independently qualified/solved seed. Named sensitivity drivers supply their own P1; the generic API cannot infer one. | `run_incremental_crack_path.m:65`; `run_incremental_crack_path.m:114`; `verification/crack_path/main_increment_sensitivity_coarse.m:227` |
| F16 | important | open | An absolute calibration transitionLength_m overrides the dimensionless transition option; summary transition fields describe the nominal design. Named M1 protocol overrides slopes only. Custom callers must inspect effective calibration, not only requested ratios. | `build_stage3c_polyline_exterior.m:287`; `build_stage3c_polyline_exterior.m:297`; `qualify_incremental_crack_candidate.m:198`; `main_stage2_embed_scaled_core_full_domain_theta0.m:147` |
| F17 | confirmed safe | verified | 92 mm is a default target, not a hidden segment count: N=target/increment with explicit divisibility. 4/2/1 mm gives 23/46/92 segments. Nondivisible targets fail; no shortening or floor conversion. | `verification/crack_path/main_increment_sensitivity_coarse.m:51`; `verification/crack_path/main_increment_sensitivity_coarse.m:68` |
| F18 | confirmed safe | verified | Accepted Stage-I geometry/material/load/frame stay unchanged; only the in-memory reserved Stage-II increment changes. Full R0 restores bitwise after resetting that one summary cell in the tests. | `verification/crack_path/main_increment_sensitivity_coarse.m:78`; `verification/crack_path/main_increment_sensitivity_coarse.m:93`; `crack_physics_signature.m:1` |
| F19 | confirmed safe | verified | Historical C.a0=4 mm does not leak into the active P1 carrier: A0Override explicitly forwards the selected reservation. Do not mutate accepted C merely to remove a seemingly stale literal. | `main_stage2_embed_scaled_core_full_domain_theta0.m:87`; `main_stage2_embed_scaled_core_full_domain_theta0.m:118`; `build_stage2_cracked_mesh_for_theta.m:34`; `geom_hole_shortcrack.m:52` |
| F20 | confirmed safe | verified | Tip size, annulus and core radii use the selected increment; CoreScale changes tip size, not physical core/EDI radii. Coarse 1/2/4-mm core element and COD sampling invariance verified without a PDE/physical solve. | `build_stage2_scaled_audited_core.m:65`; `build_stage2_scaled_audited_core.m:118`; `qualify_incremental_crack_candidate.m:192`; `qualify_incremental_crack_candidate.m:210` |
| F21 | confirmed safe | verified | 4/2/1 runs intentionally couple increment and absolute core/exterior mesh scales; they are not fixed-absolute-mesh experiments. Keep this existing interpretation; compare the same dimensionless coarse family. | `verification/crack_path/main_increment_sensitivity_coarse.m:5`; `verification/crack_path/main_increment_sensitivity_coarse.m:25`; `verification/crack_path/main_increment_sensitivity_coarse.m:96` |
| F22 | important | open | Source-carrier Hmin/Hmax/Hgrad, mouth width, and physical hole polygon are inherited absolute quantities, not globally rescaled with Delta a. Core/exterior laws scale; the full plate/carrier/boundary seed problem is not a homothetic mesh family. | `qualify_incremental_crack_candidate.m:133`; `qualify_incremental_crack_candidate.m:138`; `build_stage2_cracked_mesh_for_theta.m:39`; `build_stage3c_polyline_exterior.m:324` |
| F23 | important | open | Qualified fingerprints are duplicated in qualifier, P1 solver, tip solver, study guards and tests. Current values agree. A new qualified family requires coordinated updates, not replacement of only one count. | `qualify_incremental_crack_candidate.m:219`; `main_stage2_embed_scaled_core_full_domain_theta0.m:160`; `solve_incremental_crack_tip.m:720`; `main_stage2_theta0_physical_solve.m:677` |
| F24 | confirmed safe | verified | CoreScale API limits are explicit: propagation supports 1/2; fixed-geometry qualifier/solver also support 0.5; display-only scale 4 is not physically admitted. No silent fallback to scale 1. Do not treat visualization H2 as a physically qualified propagation level. | `run_incremental_crack_path.m:43`; `qualify_incremental_crack_candidate.m:28`; `solve_incremental_crack_tip.m:752`; `verification/crack_path/main_tip_core_mesh_level_comparison.m:13` |
| F25 | cosmetic | open | H2 names mean scale 4 in mesh-level illustration, but h/2 refinement documentation uses another naming convention. Use numeric CoreScale and hTip, not level names, as identity. | `verification/crack_path/main_tip_core_mesh_level_comparison.m:6`; `verification/crack_path/main_fixed_geometry_tip_refinement_diagnostic.m:117`; `verification/crack_path/TIP_REFINEMENT_H2.md:1` |
| F26 | cosmetic | open | MeshFamilyLabel and output-directory tags are display/organization labels; the generic default can still say reference for nonreference numerical controls. Compatibility uses numerical controls, not labels. Do not infer provenance from filenames. | `run_incremental_crack_path.m:53`; `run_incremental_crack_path.m:397`; `verification/crack_path/main_increment_sensitivity_coarse.m:103` |
| F27 | confirmed safe | verified | P17/P21/P22/P23, 23 rows and 92 mm in canonical paper generators are deliberately tied to the audited 4-mm reference snapshot. Snapshot guards must not be weakened to accept 2/1-mm histories. These are not generic comparison tools. | `paper/export_pilot_evidence.m:12`; `paper/export_pilot_evidence.m:20`; `paper/sync_manuscript_data.py:16`; `paper/plot_existing_manuscript_figures.m:62`; `paper/verify_manuscript.py:16` |
| F28 | important | open | Fixed-geometry diagnostics pin reference indices 17/22 and 4*k lengths; custom source files are not a generic physical-length selector. For another increment derive the desired physical length/event from that history. Do not reuse these reference-index defaults. | `verification/crack_path/main_fixed_geometry_tip_refinement_diagnostic.m:24`; `verification/crack_path/main_fixed_geometry_tip_refinement_diagnostic.m:139`; `verification/crack_path/main_fixed_geometry_tip_coarsening_diagnostic.m:141`; `verification/crack_path/main_fixed_geometry_coarse_exterior_diagnostic.m:28` |
| F29 | confirmed safe | verified | The 4/2/1 workflow matches common physical lengths rather than equal segment indices and does not interpolate paths. Intersection uses absolute 1e-10 mm; nearest-row check retains 1e-8 mm. Ambiguous duplicate grids now reject. | `verification/crack_path/plot_increment_sensitivity_coarse_publication.m:75`; `verification/crack_path/plot_increment_sensitivity_coarse_publication.m:478` |
| F30 | confirmed safe | verified | The partial-level publication path uses 4/2 data until both 1-mm CSVs exist; a live atomic state is not treated as an exported complete study. Two- and three-level branches tested with explicit synthetic fixtures; no physical results fabricated. | `verification/crack_path/plot_increment_sensitivity_coarse_publication.m:45`; `verification/crack_path/main_increment_sensitivity_coarse.m:273` |
| F31 | important | open | Resume is continuation, not an idempotent completed-run load: default Resume=true rejects an already completed unchanged target. No silent extra solve occurs. Changing this continuation contract is outside the mechanical audit fixes. | `run_incremental_crack_path.m:608`; `verification/crack_path/main_increment_sensitivity_coarse.m:259` |
| F32 | important | open | Assembly depends on a BN_local helper resolved from MATLAB path outside this checkout. Record/pin MATLAB/PDE/runtime/path dependencies. No dependency or solver implementation changed here. | `stif_assem.m:44` |
| F33 | confirmed safe | verified | Canonical benchmark/staged reference fingerprints are explicitly reference-only, with an alternative-candidate branch for other qualified families. Do not use reference benchmark/regression entry points as generic increment tests. | `verification/crack_path/main_incremental_path_profile.m:20`; `verification/crack_path/main_incremental_path_profile.m:69`; `main_stage2_theta0_physical_solve.m:125`; `main_stage2_theta0_physical_solve.m:789` |
| F34 | cosmetic | open | Reference manuscript availability wording is not synchronized with the investigator-local exported 4/2-mm comparison. Author review must distinguish partial/coupled 4/2 data from a complete or isolated increment study. No manuscript interpretation was edited. | `paper/main.tex:818`; `verification/crack_path/plot_increment_sensitivity_coarse_publication.m:45` |
| F35 | confirmed safe | verified | Historical 8-mm audit constants and fixed Stage-I-source hashes are validation/snapshot identities, not production increment substitutions. No historical data were injected into another increment's seed or trajectory by this audit. | `build_stage2_scaled_audited_core.m:10`; `paper/audit_stage1_source_provenance.py:24`; `verification/sif_audit/AUDIT_CLOSURE.md:1` |

## Complete active parameter flow

The reference Stage-I source is loaded from
`paper/data/accepted_stage1_source.mat` by the named independent sensitivity
drivers. Stage I selects the mouth/normal/tangent and the hole-only elastic
record. Its reserved first-segment length is a Stage-II choice.

| Parameter | Source -> consumer | Result in the named coarse increment study |
|---|---|---|
| IncrementMM | study parser at `main_increment_sensitivity_coarse.m:51` -> `da=1e-3*IncrementMM` at line 68 -> in-memory `R0.summary.a0_reserved_m` at line 93 | 0.001 / 0.002 / 0.004 m |
| P1 carrier length | Stage-II embed reads summary at `main_stage2_embed_scaled_core_full_domain_theta0.m:87` -> explicit A0Override at line 118 -> adapter at `build_stage2_cracked_mesh_for_theta.m:34` -> geom's explicit a0 at `geom_hole_shortcrack.m:52` | Selected increment, not historical C.a0 |
| Subsequent carrier | `run_incremental_crack_path.m:94` -> coordinates incremented at line 342 -> qualifier's complete supplied polyline at `qualify_incremental_crack_candidate.m:133` | Full evolving polyline; last segment length validated |
| CoreScale | named study's fixed coarse scale 2 at `main_increment_sensitivity_coarse.m:96` -> P1 builder and runArgs -> qualifier at `qualify_incremental_crack_candidate.m:192` | Same scale 2 at all three increments |
| Base/tip size | `build_stage2_scaled_audited_core.m:65` sets hBase=HTipOverA0*a0; line 118 records hTip=Scale*hBase | hTip/Delta a = 2*0.00675308135 |
| Interaction/core radii | qualifier lines 210–215 and Stage-II embed lines 157–178 | ri=.10*Delta a, ro=.65*Delta a, rc=.75*Delta a; CoreScale does not rescale these radii |
| Transition / cap | named driver ratios 1 / 1.25 -> Stage-II embed lines 147–148 and qualifier lines 198–199 | Physical lengths scale with Delta a; custom absolute calibration can override effective transition |
| Exterior size | core design carries scale, hBase and near slope -> exterior constructors | Also depends on CoreScale; not determined by cap and transition alone |
| Count / target | study lines 68–76: target/da, exact divisibility guard, rounded integer, minimum two -> run MaxSegments | Target 92 mm gives 23 / 46 / 92 |
| P1 SIF seed | compatible P1 field/compact -> `main_increment_sensitivity_coarse.m:227` -> SeedKI/SeedKII at lines 252–253 | Own first-tip field per increment/family; no reference-SIF substitution |
| Output | default increment-specific da_1mm/da_2mm/da_4mm directories at study line 103; same directories for resumable files | Labels organize files; metadata/geometry establish identity |

### Lightweight measured propagation

These are geometric/parameter checks, not physical crack-path results.

| Delta a (mm) | CoreScale | hTip (mm) | ri (mm) | ro (mm) | rc (mm) | Transition (mm) | M1 cap (mm) |
|---:|---:|---:|---:|---:|---:|---:|---:|
| 1 | 2 | 0.0135061627 | 0.1 | 0.65 | 0.75 | 1 | 1.25 |
| 2 | 2 | 0.0270123254 | 0.2 | 1.3 | 1.5 | 2 | 2.5 |
| 4 | 2 | 0.0540246508 | 0.4 | 2.6 | 3 | 4 | 5 |

All three coarse core-only constructions had 3,318 T3 elements and native
sampling 19/28/23/18. The actual preserved accepted archives provide the
existing 2,976 EDI-support fingerprint. No new physical support/SIF result
was manufactured.

## Stage-I reuse census and isolation

The complete code occurrence list is `stage1_reuse_sites.csv`. The active
consumer groups are:

- Independent increment study: `main_increment_sensitivity_coarse.m:78`
  loads the accepted source and changes only its copied summary reservation.
  Resetting that one cell restores full R0 equality; C remains equal exactly.
- Independent M1 and 2h0 paths: `main_m1_independent_trajectory.m:91` and
  `main_tip2h0_independent_trajectory.m:92` pass the accepted state to their
  own P1 qualifier/solver and then pass their own SIF seed onward.
- Generic propagation: `run_incremental_crack_path.m:84`,
  `qualify_incremental_crack_candidate.m:50`,
  `solve_incremental_crack_tip.m:31`: the mouth/frame/material/loading are
  supplied from R0 while the reservation sets the Stage-II increment.
- P1 route: `main_stage2_embed_scaled_core_full_domain_theta0.m:69` and
  `main_stage2_theta0_physical_solve.m:105`: the explicit override prevents
  the retained C.a0=4 mm from determining a 1/2-mm carrier.
- Fixed-geometry diagnostics: refinement/coarsening/exterior entry points
  load the reference evidence and accepted frozen source, then qualify an
  existing path. Their fixed indices are reference-study selectors.
- Historical staged Stage-II/III and paper exports retain their explicit
  reference contracts. For example `main_stage2_embed_scaled_core_full_domain.m:83`,
  `main_stage3a_two_leg_tip_core_qualification.m:83`,
  `main_stage3c_kinked_two_leg_qualification.m:82`,
  `main_stage3d_kinked_two_leg_physical_solve.m:109`,
  `paper/export_pilot_evidence.m:20`. They are not called by the coarse
  increment orchestrator as generic increment-history loaders.

Unchanged Stage-I values include plate/hole geometry, material, loading,
anchors, selected mouth, normal/tangent, stress/fit record, and boundary
sampling. C.a0 retains its historical reserved value intentionally.
The Stage-II increment changes P1 geometry and every future segment, core
size, EDI radii, transition and cap in the dimensionless study. Existing
absolute source-carrier mesh parameters and hole-boundary discretization
remain fixed; consequently the entire plate mesh is not a homothetic
scaled copy at each increment.

## Critical open isolation issue: exterior dependence on CoreScale

At `build_stage3c_polyline_exterior.m:77` and
`build_stage2_scaled_audited_exterior.m:39`, both
`scale=design.scale` and `hRp=scale*(hBase+slope*rp)` are used. The exterior
function multiplies its radial increment by scale as well.

For the unchanged reference 4-mm settings, rc=3 mm, transition=4 mm,
far cap=2.5 mm, evaluating that same size law gives:

| CoreScale | h at core boundary (mm) | h at radius 12 mm (mm) |
|---:|---:|---:|
| 0.5 | 0.0555061627 | 0.4048754732 |
| 1 | 0.1110123254 | 0.7949096269 |
| 2 | 0.2220246508 | 1.4738318921 |

These are mesh-target evaluations, not stresses or new solutions.
Thus fixed cap/transition does not isolate the tip-core resolution from the
exterior target field. The diagnostic description/boolean
`exteriorMeshLawChanged=false` at
`main_fixed_geometry_tip_coarsening_diagnostic.m:272` and the tip-only
attribution in `paper/main.tex:730` / `paper/main.tex:832` require author
review against this coupling. Neither these claims nor the scientific law
was rewritten in the audit. Also, remeshing with the same numerical law
does not imply identical exterior connectivity.

## Checkpoint/resume conclusions

- Different increments are rejected by every segment length and seed geometry;
  full fields also compare the stored increment. The P1 shortcut now checks
  absolute dimensions and actual source geometry, rather than scale alone.
- Different cores use the existing qualified fingerprint table. Coarse
  accepted-result promotion is now possible with 2,976 support elements;
  using 11,316 indiscriminately was a dispatch bug, not a stricter valid gate.
- Different exterior families are compared by cap, transition, override
  structure and core scale. NaN identity values reject. New compact results
  carry explicit exterior controls and physical signatures; older compacts
  need their source candidate/field archives, loaded without U for promotion.
- Full checkpoint reuse compares actual T3/T6 geometry and current physical
  data/material. Existing coordinate tolerance remains 1e-12; physical
  parameters are identity values rather than fitted quantities.
- Declared historical schema/row order is checked. New atomic states record
  frozen physical identity. **Legacy fully embedded states remain trusted
  historical records if they lack this identity.** This audit cannot recover
  unrecorded provenance or certify arbitrary edited/copied old arrays.
- Old metadata-free reference-candidate inference remains a documented
  backward-compatibility assumption. Directory names alone are never proof.
- Identity rejection does not make output folders transactional: a fresh
  qualification may already have written files. Keep one directory per
  scientific configuration; archive transaction redesign was not performed.

## Publication and common-length matching

The 4/2/1 workflow is specialized to those increment values. It uses the
4-mm physical-length grid intersected with the finer grids, with absolute
`ismembertol(...,1e-10,'DataScale',1)` in mm and the unchanged 1e-8-mm
row-match check. It does not interpolate trajectories. An all-common
finer grid is intentionally not promised. Both two-level and three-level
availability paths were exercised with synthetic 12-mm straight-path
fixtures, obtaining common lengths [4,8,12] mm and the prescribed constant
turn density. These fixtures have no physical meaning and are not committed
as research results.

The loader now verifies the coarse-family study summary, increment, controls,
unique sequential finite rows, vertex lengths and geometry, and concordance
with the authoritative summary. Reordered segment-index comparisons in the
M1/2h0 reference helpers also require a common uniform increment.
Canonical committed-panel fallback no longer overrides explicit alternative
history requests.

A read-only snapshot at baseline checking found the 4-mm and 2-mm coarse
states completed through 23 and 46 segments, and the 1-mm state through
51 accepted segments with one appended geometry. This is a snapshot, not a
claim that a 1-mm 92-mm study is finished. The normal export protocol creates
the final CSV/summary only after the driver returns; this audit did not
complete or continue any run.

Canonical paper export/sync/verification scripts deliberately require the
4-mm P1–P23 reference, its 92-mm endpoint and known events. Their fixed
P17/P21/P22/P23 labels are not generic event selectors. Keep those guards for
the reference snapshot; derive physical-length/state-event selectors in any
new increment-study exporter instead of weakening reference checks.

## Safe fixes and validation

`test_crack_reproducibility` passed **15 portable checks and 22
investigator-archive checks (37 total)**. Tests cover:

- Stage-I copy equality and 1/2/4-mm core size/topology/native sampling;
- source-level A0Override forwarding, target divisibility and minimum length;
- rejection of different E, nu, ps, plate dimensions and load in cache identity;
- real compatible coarse P3 promotion at both 2 and 4 mm, retaining stored
  KI/KII exactly, without a new solve;
- rejected wrong increment, core, NaN family, unknown schema and reordered rows;
- incompatible P1 seed reuse, exterior provenance, public material/T6 checkpoint
  rejection before EDI/postprocessing;
- real CSV summary loading, malformed family metadata and duplicate grids;
- two-/three-level publication availability and exact common-length matching.

The test's archive-specific section requires investigator-local saved data
and reports a skip when unavailable. Portable checks and synthetic publication
fixtures remain runnable from the committed accepted Stage-I source.
Existing MATLAB Code Analyzer advisories remain; no syntax failure was found.

Numerical criteria remain unchanged: PCG 1e-10 / 5000, true residual 5e-10,
constraints 1e-14, 16-point EDI, geometric/mesh gates and the existing
family-specific sampling/support tables. Reuse/schema/input identity checks
were strengthened; no thresholds were relaxed. No physical interpretation
or manuscript/PDF was altered.

Reproduce the checks with an isolated writable directory:

```matlab
addpath(genpath(pwd));
Report = test_crack_reproducibility('WorkDir','<new audit scratch directory>');
```

Reproduce the lexical search:

```text
python verification/crack_path/reproducibility_audit/scan_repository.py --output-dir <audit-output>
```

## Required follow-up decisions, not performed

1. Decide whether the tip-resolution study is a core-plus-graded-exterior
   family or redesign a genuinely fixed-exterior control. This requires
   author/scientific approval and potentially new qualification/solves.
2. Retain/prove missing legacy origin metadata before trusting imported
   historical complete rows for another configuration.
3. Use isolated output roots and consider pre-write archive identity/atomic
   transaction design in a separate change.
4. Reconcile publication availability wording with actual completed levels
   and distinguish coupled increment/mesh-family convergence from a pure
   increment study.
