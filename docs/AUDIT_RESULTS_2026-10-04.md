# SIF audit results collected on 4 October 2026

Consolidation branch: `audit/2026-10-04-results`.

The branch starts from `sif-asymmetric-mesh-audit` at
`227004483f6f81e38081378037628beee6296686`. It therefore already contains all
audit source/documentation merges made on 4 October through Step60, including
the Step56 corrections, Step57 paired symmetric solve, Step58 matched EDI,
Step59 exact-nodal replay, and Step60 real-mesh visualization. This collection
adds the previously local results and Step61 implementation; it does not
repeat those experiments or change their historical provenance.

## Collected results

| Study | Included artifacts | Interpretation |
|---|---|---|
| Step56 | Exact reflected T3 candidate and small-data MAT | Symmetric paired mesh preflight |
| Step57 | Existing symmetric solved checkpoint and small-data MAT | Actual symmetric FEM/COD control, performed earlier today |
| Step58 | Matched EDI result, comparison and final progress MAT | Actual symmetric solved-field extraction |
| Step59 | Exact-nodal replay result, comparison and final progress MAT | Prescribed-field extraction control |
| Step60 | Both actual Step38 mesh PNGs and small-data MAT | Read-only visualization and primary support audit |
| Step61 | Driver, report, exact candidate T3 MAT, overview/tip PNGs, small-data MAT and seven CSV tables | Local pairing in unchanged asymmetric geometry; mesh and prescribed controls only |
| Scientific review | [Full review](SIF_asymmetry_scientific_review_2026-10-04.md) | Historical read-only assessment of Steps34--59 at its stated audited commit |

The review's evidence-scope statements describe what was available during
that review. Subsequent local checkpoint access and Step61 qualification are
documented separately in [STEP61.md](../verification/sif_audit/STEP61.md).

The exact collection inventory and SHA-256 hashes are in
[AUDIT_RESULTS_2026-10-04_manifest.csv](AUDIT_RESULTS_2026-10-04_manifest.csv).
The manifest covers 27 result/source files, this index and `.gitattributes`;
it excludes itself to avoid a self-referential checksum. Narrow path-specific
attributes preserve the original bytes of the collected text/CSV files
across checkout platforms, so their provenance hashes remain verifiable.

## Step61 outcome

- All mesh structural gates pass on the complete primary 0.8--5.2 mm
  FE-nodal-q support, including straddling elements and literal production
  participation.
- Paired radius 5.6 mm; collar approximately 0.650--0.842 mm; edge-aligned
  cavity maximum radius 6.441771 mm.
- 19,811 original triangles replaced and 19,630 exterior triangles retained
  exactly. Total T3 count is 38,987 versus 39,441; total T6 nodes 78,873
  versus 79,769.
- Tip target 0.05402465 mm matched; affected minimum angle improves from
  10.37 to 24.05 degrees. Shell-level resolution differences are tabulated.
- Affine, pure-I, pure-II and tiny mixed prescribed controls pass. Pure-I
  cross-mode leakage is `6.6035e-15`; unit-II recovery is `0.9999999802`;
  tiny mixed relative KII error is `1.9727e-8`.
- Ready to propose one asymmetric FEM solve. **No asymmetric solve was run.**

## Preservation and reproduction

The consolidation performed Git/file operations only. It did not assemble
stiffness, run FEM/EDI/COD, rerun MATLAB, or alter production code, Stage I,
the existing physical problem, or the Step61 candidate. The existing Step57
solved checkpoint is archived as an earlier symmetric control; Step61's MAT
artifacts contain no physical displacement field.

The original Step61 driver and recorded hashes remain unchanged. Its branch
guard intentionally requires `sif-asymmetric-mesh-audit`; review/merge these
sources into that study branch before reproducing the experiment there.
Step61 additionally requires the exact investigator-local
`verification/step38_tip_refined_solved.mat`, with checkpoint SHA-256 recorded
in STEP61. That older solved checkpoint is not duplicated in this collection.

Older local result files and document/compiler artifacts are left untouched.
`main` and the original `sif-asymmetric-mesh-audit` branch are not advanced by
this consolidation.
