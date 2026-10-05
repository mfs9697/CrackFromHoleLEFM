# SIF audit chronology and file map

The asymmetric SIF audit is scientifically closed. This file explains the
historical numbering and where each stage is documented so that apparent gaps
are not mistaken for lost work.

The immutable pre-cleanup snapshot is preserved on:

`archive/sif-audit-closed-2026-10-05`

at commit `becd2d549f7bcfd406849e38c98e01b82610b3f5`.

## Numbering map

| Audit stage | Canonical record in the cleaned branch | Note |
| --- | --- | --- |
| Steps 1--38 | `docs/SIF_ASYMMETRIC_MESH_AUDIT.md` plus the corresponding `main_step*.m` drivers | Original continuous forensic log |
| Steps 39--43 | Retrospective summary in `STEP44.md`; `main_step39_tip_cod_sensitivity.m` is retained | These were exploratory follow-ups, not five missing final result documents |
| Step 44 | `STEP44.md` | Exact Williams replay |
| Steps 45--47 | `STEP45.md`--`STEP47.md` | Symmetric-control chain |
| Step 48 | Result and interpretation are recorded at the start of `STEP49.md`; driver `main_step48_refined_matched_edi.m` is retained | No standalone `STEP48.md` was created |
| Step 49 | `STEP49.md` | Reflection-parity audit |
| Steps 50--51 | Measured results are recorded inside `STEP49.md`; the existing Step45/45c replay drivers were reused | No standalone STEP50/STEP51 drivers or documents were required |
| Steps 52--60 | `STEP52.md`--`STEP60.md` | Continuous documented chain |
| Step 61 | `STEP61.md` and `main_step61_local_paired_patch.m` | Restored from the historical Step62 branch during cleanup; mesh-only paired-patch experiment |
| Step 62 | `STEP62.md` | Structured graded reference family |
| Step 62B | `STEP62B.md` | Exterior calibration |
| Step 63 | `STEP63.md` | Level-0 physical solve |
| Step 63R | `STEP63R.md` | Recovery of the lost Level-0 checkpoint |
| Step 64 | `STEP64.md` | Matched physical EDI |
| Step 65 | `STEP65.md` | Level-1 mesh qualification |
| Step 66 | `STEP66.md` | Direct-solver memory preflight |
| Step 67 | `STEP67.md` | Fixed ICT attempt; failed before PCG |
| Step 67A | `STEP67A.md` | Parameter-free SGS-PCG qualification |
| Step 68 | `STEP68.md` | Level-1 physical convergence |
| Step 69 | `STEP69.md` | Additional (s=1/\sqrt2) physical point and three-scale analysis |
| Step 70 | `STEP70.md` and `AUDIT_CLOSURE.md` | Synthesis and scientific closure |

## Why Step61 matters but is not part of the final convergence proof

Step61 constructed a locally reflection-paired patch inside the real
asymmetric Step38 geometry. All structural and prescribed-field extraction
controls passed, including essentially exact recovery of a prescribed
(K_{II}/K_I=10^{-4}) signal.

The candidate was intentionally **not** used as the final physical convergence
family because pairing, tip connectivity, and shell-by-shell resolution all
changed together. Step62 therefore replaced it with a deterministic structured
family in which the global scale could subsequently be varied in a controlled
way.

Step61 remains valuable as a feasibility and design experiment and is retained
for provenance.

## Active versus historical material

Everything under `verification/sif_audit/` should now be regarded as the
**closed verification record** unless a file explicitly says otherwise.

The audit drivers are not the starting point for the next crack-path study.
The next research branch should begin from the hole-initiation problem and use
the qualified interaction-EDI methodology to determine the first-segment
direction from a new (K_{II}(\theta)) calculation.

No new physical FEM result should be added to this closed audit chronology as
"Step71".
