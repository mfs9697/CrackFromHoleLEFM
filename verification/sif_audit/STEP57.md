# Step 57: one explicitly authorized FEM solve on the exact Step56 reflected candidate

## Why this step is now scientifically justified

The mesh-only development phase has reached a clean stopping point. The investigator's Step56 run produced an exact, saved, **reflection-paired** T3 candidate that preserves the previously verified Step54/55 physical boundary and crack geometry while correcting the mesh-topology asymmetry observed in Steps49–53.

The candidate is not merely symmetric in element counts. Step56 explicitly verified:

- complete upper/lower T3 reflection;
- complete six-node T6 element reflection after the global T3→T6 upgrade;
- one shared node set on the intact horizontal ligament ahead of the crack tip;
- distinct upper/lower crack-face node IDs behind the tip, despite coincident post-collapse coordinates;
- shared T6 midside nodes on the intact ligament and distinct T6 midside nodes on the crack faces;
- symmetric 3/3 tip fans;
- zero native upper/lower crack-face radial mismatch;
- positive collapsed element areas and a minimum T3 angle of 20.095°;
- native crack-face sampling 72/72 with 19/16/13 points in the predeclared windows.

The candidate has **2578 T3 elements and 5427 T6 nodes**, with median tip edge 0.00014363 m. It therefore differs from the prior Step47 mesh (2391 T3, 5054 T6, median tip edge 0.00013535 m). A new solve would isolate the effect of using a deliberately reflection-paired topology much better than another uncontrolled remesh, but it would **not** be a pure h-refinement comparison.

## Guarded driver

`main_step57_reflection_paired_control.m` is deliberately **default-off**:

```matlab
[P57,O57] = main_step57_reflection_paired_control();
```

With no existing Step57 checkpoint this command validates the locally saved Step56 candidate and then raises `step57:ExplicitSolveApprovalRequired`. That is intentional.

Before any solve, the driver reloads the exact saved candidate and compact Step56 report and independently rechecks:

- Step56 solve-proposal gates;
- candidate/provenance identity;
- T3/T6 element/node counts;
- positive signed T3 areas;
- exact native radial grid and the 19/16/13 COD-window counts;
- the 3/3 tip fan and stored median tip edge;
- distinct coincident crack-face T3 nodes with only the tip shared;
- no duplicated coincident node IDs along the intact ligament;
- unchanged crack geometry, material, loading and half-domain symmetry conditions relative to Step47.

The FEM solver receives **exactly** `candidate.p` and `candidate.t`. The driver does not regenerate a PDE mesh. After the solve it verifies that the solver's T3 and T6 coordinate/connectivity arrays exactly match the accepted candidate before saving a unique checkpoint.

If, and only if, the investigator explicitly authorizes one new FEM solve, the command will be:

```matlab
[P57,O57] = main_step57_reflection_paired_control( ...
    'AllowSolve', true);
```

After that single solve, the driver automatically performs only the previously audited **native COD** postprocessing. It does **not** perform EDI. The intended stopping rule is to inspect the raw/fitted COD symmetry residuals first. Any matched EDI calculation remains a separate later decision.

The unique solved checkpoint, if authorized and successfully generated, will be:

`verification/step57_reflection_paired_symmetric_theta0_solved.mat`

and the COD-only compact output will use:

`verification/step57_reflection_paired_symmetric_field_leakage*`

An existing valid checkpoint is always reused; it is never overwritten or silently regenerated.

## Interpretation after a future authorized solve

Because the continuum control is exactly symmetric, physical `KII=0`. If the reflection-paired candidate sharply reduces the actual solved-field tangential COD/parity residuals, that would support mesh-topology asymmetry as one contributor to the prior numerical leakage. If it does not, the result would instead weaken that mechanism.

Neither outcome alone validates the tiny asymmetric physical Mode-II signal. The reflection-paired mesh has different connectivity from Step47 and a slightly different achieved tip scale. Therefore any later EDI comparison must use the **same absolute annulus** as Step45/47/48 and must remain a numerical-control comparison rather than an uncertainty correction.

## Current authorization state

**AUTHORIZED by the investigator on 2026-10-04:** exactly **one** new symmetric Step57 FEM control solve on the already saved Step56 reflection-paired candidate is approved.

The authorization is deliberately narrow:

- use the exact saved `verification/step56_reflected_mesh_only_candidate_T3.mat` candidate; **no remeshing or regeneration**;
- run the single guarded solver call through `main_step57_reflection_paired_control('AllowSolve',true)`;
- checkpoint the solved field only after the driver verifies the solver used the exact accepted T3/T6 candidate;
- run the driver's predeclared **native COD-only** postprocessing;
- **do not run EDI** in Step57;
- stop after printing the COD fit table and COD gates for interpretation.

Any later matched-domain EDI calculation requires a separate decision after reviewing this one solved control.



## First authorized attempt: pre-solve validator bug, no FEM solve occurred

On the first investigator-authorized Step57 invocation, the driver reproduced the saved Step56 native crack-face geometry (`72/72`, zero radial mismatch) and then stopped **before the solver call** with:

`Coincident upper/lower crack-face T3 nodes are not distinct as required.`

This was traced to a validator implementation error, not to the Step56 mesh. The validator applied `unique()` independently to upper and lower crack-face node IDs and then compared coordinates row-by-row. Because the reflected lower face intentionally uses **different node numbers**, sorting each ID set numerically destroys the physical upper/lower correspondence. Step56 had already established the intended topology.

The corrected validator now:

- requires no repeated IDs within either face list;
- requires equal upper/lower face-node counts;
- requires that the two ID sets share **exactly one node, the common crack tip**;
- projects each face's physical coordinates onto the saved crack tangent;
- sorts by that physical crack coordinate, not by node number;
- requires upper/lower collapsed coordinates and crack parameters to agree to tolerance.

No geometry, mesh recipe, candidate file, solver call, loading, boundary condition, or postprocessing rule changed.

**The authorized solve has therefore not yet occurred and the one-solve authorization remains available.** Rerun the same guarded command after pulling this fix. EDI remains prohibited in this step.


## Investigator's completed Step57 solve: reflection-paired FEM COD leakage collapses to roundoff

The investigator completed the **single authorized** Step57 symmetric FEM solve on the exact saved Step56 candidate. The driver verified and reused the saved reflected mesh **without remeshing**: `T3=2578`, `T6=5427`, native crack-face samples `72/72`, zero radial mismatch, COD windows `19/16/13`, symmetric `3/3` tip fan and median tip edge `0.000143632036 m`. The unique checkpoint was saved as `verification/step57_reflection_paired_symmetric_theta0_solved.mat`.

The genuine solved-field native crack-face ratios were all at floating-point-noise scale. Representative median raw `KII/KI`-like opening ratios were:

| r/a0 band | median raw tangential/normal opening ratio |
| --- | ---: |
| 0–0.04 | -9.3293e-14 |
| 0.04–0.08 | -1.1513e-13 |
| 0.08–0.12 | -1.0056e-13 |
| 0.12–0.20 | -8.5134e-14 |
| 0.20–0.30 | -7.9903e-14 |

All six predeclared COD fits were qualified. Their signed fitted ratios ranged from about `-9.38e-14` to `-1.84e-13`; the maximum finite fitted magnitude was **`1.83854063e-13`**, far below the existing `1e-6` numerical target. Relative to the Step47 maximum fitted COD leakage `6.27056963e-4`, this is approximately a **3.41e9-fold reduction**.

At the same time the Mode-I COD intercepts remained very close to the prior Step47 values. For the three linear fits, Step57 returned `KI=0.35212, 0.35499, 0.35618`, differing from the corresponding Step47 linear values by only about **0.068%, 0.020%, and 0.006%**, respectively. Thus the disappearance of the tangential opening is not accompanied by a comparable change in the dominant Mode-I response.

**Interpretation:** this is strong evidence that lack of reflection-paired finite-element topology was the dominant source of the previously observed **native-COD Mode-II leakage** in the symmetric benchmark. Because the Step57 mesh is not identical to Step47 (different connectivity and slightly different achieved tip scale), this is not a formal one-variable convergence proof, but the contrast is several orders of magnitude larger than those secondary mesh differences. The physical continuum expectation `KII=0` is recovered by the reflected solved field to numerical roundoff in the crack-face opening diagnostic.

This result does **not** yet determine whether the existing 16-point FE-nodal-q EDI extractor will also cancel to roundoff on the reflection-paired solved field. Earlier exact-field tests showed that EDI has its own interpolation sensitivity. A later EDI should therefore be treated as a **separate matched-annulus extractor control**, not as a correction factor or asymmetric-case uncertainty bound.

The final MATLAB line `disp(O57.rawTable)` produced an error **after all Step57 work and saves had completed** because the postprocessor's actual field name is `O57.rawBands`. No solve or data was lost. The correct display is:

```matlab
disp(O57.rawBands);
disp(O57.fitTable);
disp(O57.CODgates);
```

No EDI was performed in Step57.
