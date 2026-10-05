# Step 46: controlled same-geometry mesh-only preflight

## Why this step is needed

Our Step 45 control has a straight horizontal crack emerging from the exact rightmost point of a centered circular hole. Its **continuum** Mode-II SIF must be zero, but the **coarse actual-FEM** field gives domain-dependent nonzero EDI Mode II and too few native crack-face samples for COD extrapolation. Its crack-tip T3 fan is also not reflection-paired: two triangles above the crack and three below.

Step 45b/c gave a direct **same-mesh** mechanism test using prescribed analytical pure Mode-I nodal and Gauss-point fields:

| (r_o/a_0) | Actual FEM signed (K_{II}/K_I) | Exact nodal signed ratio | Exact at Gauss points signed ratio |
| ---: | ---: | ---: | ---: |
| 0.50 | `-1.563466231e-5` | approx. `-4.81123e-5` | approx. `-1.735158e-7` |
| 0.65 | `-8.766949318e-6` | approx. `-3.88687e-5` | approx. `-1.147788e-7` |

The exact Gauss-point tests had `KI=0.999999035` at 0.50 and `0.9999993666` at 0.65. Compared to **exact nodal** fields, exact Gauss-point evaluation reduces apparent Mode II leakage by roughly **277×** and **339×**, respectively. These results strongly implicate ordinary T6 interpolation of the prescribed singular displacement field on **this coarse mesh**. They do **not** quantify FEM-field error or independently prove small physical Mode II.

We should not spend more integration cycles on prescribed exact fields at further radii until we have a better-resolved actual FEM control.

## Single first objective: inspect a proposed refined mesh **without solving FEM**

The `main_step46_symmetric_mesh_preflight.m` driver constructs:
1. The **previously saved** Step 45 baseline control mesh from the local MAT checkpoint: no new computation of its displacement field.
2. A newly generated *baseline* control mesh, using the existing `build_stage2_centered_half_cracked_mesh_for_theta` procedure; this identifies the precise sharp-pencil vertex/edge IDs. **Geometry/mesh generation only; no solve.**
3. A proposed manually refined mesh from the **exact same PDE geometry-description struct `D`** and the same `Npoly=240`, `a0=0.004 m` horizontal crack. The old global `Hmax`, `Hgrad`, and `Hmin` options are kept unchanged. Only the requested local `Hedge` and `Hvertex` parameters are halved.

The old nominal scales are `Hface=hArc` and `Htip=0.5*hArc`, with `hArc=2*pi*(0.03)/240`. The proposed nominal values are `Hface=0.5*hArc`, `Htip=0.25*hArc`. **PDE Toolbox is not guaranteed to attain these local sizes** given unchanged global `Hmin` or other mesh constraints. The preflight measures actual generated edge medians, face-node counts, crack-face grid mismatch, and upper/lower crack-tip element counts before we authorize *any* new FEM solve.

Although the boundary polygon and prescribed refinement settings are held fixed aside from local mesh sizes, PDE `generateMesh` can retriangulate regions away from the crack tip. This is **not** a proven fixed-exterior-connectivity comparison or an isolated attribution of any future FEM error change to the crack tip alone. A genuinely reflection-paired FEM mesh is a separate possible later control; counting tip triangles is a diagnostic, not proof of reflection pairing.

**Safety checks:** the driver requires the original saved Step45 symmetric solved checkpoint, checks its documented `Npoly=240` and `a0=0.004` geometry, reidentifies the sharp-pencil edges/vertex after manual refinement, verifies exact `Pmid`, upgrades each mesh to T6 **solely for face classification**, rejects degenerate collapsed T3 triangles, and saves only a compact report without any mesh, FEM stiffness matrix, or computed displacements.

## One local run after GitHub Desktop pull

Select branch `sif-asymmetric-mesh-audit`, **Fetch origin → Pull origin**. From the repository root in MATLAB:

```matlab
addpath(genpath(pwd));
O46 = main_step46_symmetric_mesh_preflight();
disp(O46.meshTable);
disp(O46.gates);
```

The existing local checkpoint is expected at `verification/step45_symmetric_theta0_solved.mat`; if it lives elsewhere, pass `'ReferenceCheckpoint',absolutePath`. The driver does NOT need `P45`, `O45`, `R45b`, or `R45c` to remain in the MATLAB workspace. It writes the small data file `verification/step46_mesh_preflight_small_data.mat`.

Send the complete MATLAB log, particularly the **three mesh rows**, crack-face sampling, median `hTip/a0` and upper/lower tip-triangle counts. The printed `gates.readyForOneRefinedFEMProposal` is a *preflight criterion*, not authorization for the solve.

Do not generate another FEM solution yet. After interpreting the mesh-only results, decide whether the proposed local changes meaningfully improve both the tip region and native COD sampling. A subsequent explicitly authorized one-solve Stage-II control would need to use the verified proposed mesh and match an EDI integration annulus to the baseline (e.g., reusing the **exact saved** original `r_inner`, not inadvertently recomputing `2hTip`).


## Investigator's first mesh-only results: half face and half tip request

The investigator executed the original zero-solve Step 46 pilot on the saved Step 45 symmetry-control checkpoint. Both the original and regenerated baseline meshes were **exactly identical in sorted T3 coordinates (maximum discrepancy 0 m)**, including 1,721 T3 and 3,640 T6 nodes. Sharp-pencil identification required the existing temporary-background-mesh fallback; it successfully recovered tip vertex 65 and edges 64, 65, the same geometry used for the refined request.

With `FaceFactor=0.5`, `TipFactor=0.5` and fixed prescribed hole/crack geometry:

| Mesh | T3 | T6 | median tip edge [m] | hTip/a0 | tip T3 above/below | native upper/lower face points | max face abscissa mismatch [m] | COD points in 0.04–0.30 / 0.08–0.30 / 0.12–0.30 |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| Saved and regenerated baseline | 1721 | 3640 | 0.00043375 | 0.10844 | 2 / 3 | 12 / 12 | 0 | 4 / 3 / 2 |
| Proposed half/half mesh | 1837 | 3884 | 0.00021079 | 0.052697 | 3 / 3 | 22 / 22 | 0 | 7 / 6 / 4 |

The tip median decreased by ~51%, and the upper/lower element *counts* are equal, but count symmetry **does not demonstrate nodal reflection pairing**. The 0.04–0.30 COD window has seven native points, **one short of the driver's prespecified eight-node linear-fit minimum** and below the twelve-node quadratic minimum. The Step 46 gate `readyForOneRefinedFEMProposal=true` only recognizes mesh improvement; it does not mean the original COD fitting goals are met.

### Completed follow-up: quarter-size face request, NO FEM solve

For the next documented trial, we held the proposed tip setting at `TipFactor=0.5` and requested denser crack-face meshing with `FaceFactor=0.25`. This changes only the *requested mesh sizes*, not the physical geometry, baseline `Npoly=240`, crack length, or FEM displacement checkpoint. The mesher can retriangulate other regions, so inspect the **achieved** tip radius and face sample count rather than assuming they are fixed.

The existing preflight driver already supports these parameters; no MATLAB code change is required. From the local repository root:

```matlab
addpath(genpath(pwd));
O46face = main_step46_symmetric_mesh_preflight( ...
    'FaceFactor',0.25, ...
    'TipFactor',0.5, ...
    'SavePath',fullfile(pwd,'verification', ...
        'step46_mesh_preflight_face025_tip05.mat'));
disp(O46face.meshTable);
disp(O46face.gates);
```

Save to this distinct filename to retain the original half/half preflight result. The main comparison is whether the refined mesh reaches **at least eight native COD points** in 0.04–0.30 without sacrificing the measured tip improvement or the exact upper/lower radial match. Twelve points would also allow the original quadratic fitting threshold. If the two proposed meshes both lack adequate sampling, decide whether to undertake a more deliberate local crack-face mesh strategy; do not infer that an additional FEM solve would establish COD convergence.

A later one-solve refined control would require a separate explicit authorization and an implementation that preserves the chosen mesh and its provenance. The present pilot remains mesh-only.

## Investigator's second preflight: quarter face / half tip

The investigator executed the previous documented mesh-only follow-up with `FaceFactor=0.25` and `TipFactor=0.5`, still using the exact original centered-hole `Npoly=240` geometry, horizontal crack length `a0=0.004 m`, and the previously saved Step45 baseline checkpoint. Automatic pure-geometry pencil-ID identification again used the established temporary-mesh fallback, which recovered edge IDs [64 65] and tip vertex 65. The freshly generated baseline T3 node coordinates reproduced the saved baseline **exactly** (maximum sorted-coordinate discrepancy 0 m; both 1,721 T3 / 3,640 T6 nodes).

| Measured property | Baseline | Half-face/half-tip trial | Quarter-face/half-tip trial |
| --- | ---: | ---: | ---: |
| T3 triangles | 1721 | 1837 | 1961 |
| T6 nodes | 3640 | 3884 | 4152 |
| Tip median edge [m] | 0.00043375 | 0.00021079 | 0.00021708 |
| `hTip/a0` | 0.10844 | 0.052697 | 0.054271 |
| Tip-adjacent T3 upper/lower | 2 / 3 | 3 / 3 | 3 / 3 |
| Native crack-face nodes, upper/lower | 12 / 12 | 22 / 22 | 36 / 36 |
| Upper/lower face abscissa mismatch [m] | 0 | 0 | 0 |
| COD nodes, `0.04–0.30 a0` | 4 | 7 | 9 |
| COD nodes, `0.08–0.30 a0` | 3 | 6 | 8 |
| COD nodes, `0.12–0.30 a0` | 2 | 4 | 6 |

**Interpretation:** the newest mesh maintains approximately **half the baseline tip-edge median**, retains a 3/3 tip-adjacent *count* and matched opposite-face native abscissae, and now qualifies for the original **linear COD extrapolation** in both the widest `0.04–0.30` and middle `0.08–0.30` windows. It does **not** meet the prespecified twelve-point minimum for any **quadratic** COD fit, or the eight-point linear minimum for the narrowest `0.12–0.30` window. Identical upper/lower **counts** do not establish exact reflection symmetry, and no FEM stress/displacement results have been calculated for these candidate meshes.

The local MATLAB `SavePath` argument used `fullfile(pwd,'verification',...)`. Because the investigator's MATLAB current folder was already `verification`, the resulting compact file was harmlessly saved under a **nested** `verification/verification` folder. In subsequent runs use `fileparts(O46face.baselineCheckpointPath)` for the actual verified output folder, independent of MATLAB's current folder.

### Final proposed mesh-only check before considering a solve

Hold the **nominal tip refinement unchanged** at `TipFactor=0.5` and request `FaceFactor=0.125` to test whether the widest window gains at least **12 actual native points** (per the prespecified quadratic threshold), while preserving a small measured tip median, upper/lower abscissa match and at least 8 samples in the middle window. PDE meshing may change the tip fan and elements outside the refined region: inspect the achieved topology; do not assume that nominal `TipFactor=0.5` preserves it.

The existing driver already supports these settings, so **no new MATLAB code or GitHub pull** is needed:

```matlab
dataDir = fileparts(O46face.baselineCheckpointPath);
O46face125 = main_step46_symmetric_mesh_preflight( ...
    'FaceFactor',0.125, ...
    'TipFactor',0.5, ...
    'SavePath',fullfile(dataDir, ...
        'step46_mesh_preflight_face0125_tip05.mat'));
disp(O46face125.meshTable);
disp(O46face125.gates);
```

Run this **mesh-only** trial once and return its complete output. Its default reported gate `readyForOneRefinedFEMProposal` permits further planning but **does not** imply that quadratic COD fits qualify. If the achieved twelve-point target is met without degrading mesh quality, stop preflight experiments and plan one explicitly authorized refined FEM solve on this *exact verified geometry and mesh recipe*. Preserve the selected mesh provenance. Later EDI comparison must use the same **absolute** annulus as the old actual-FEM baseline; do not silently recompute `r_inner=2*hTip` on the new mesh.



## Final preflight achieved: one-eighth face / half tip (SELECTED)

The investigator ran the last **mesh-only** Step46 trial with \`FaceFactor=0.125\` and \`TipFactor=0.5\`, using the saved original Stage-II symmetric checkpoint. The geometric edge-ID recovery fallback returned tip vertex 65 and face edges 64, 65. The regenerated baseline once more reproduced the original 1,721-triangle / 3,640-T6-node mesh **exactly in sorted T3 coordinates (0 m difference)**. It did not run FEM or EDI.

| Measured property | Saved original | Half face / half tip | Quarter face / half tip | **Selected one-eighth face / half tip** |
| --- | ---: | ---: | ---: | ---: |
| T3 elements | 1721 | 1837 | 1961 | **2391** |
| T6 nodes | 3640 | 3884 | 4152 | **5054** |
| Median tip edge [m] | 0.00043375 | 0.00021079 | 0.00021708 | **0.00013535** |
| hTip/a0 | 0.10844 | 0.052697 | 0.054271 | **0.033837** |
| Tip-adjacent T3 above/below | 2 / 3 | 3 / 3 | 3 / 3 | **3 / 3** |
| Native upper/lower crack-face nodes | 12 / 12 | 22 / 22 | 36 / 36 | **72 / 72** |
| Upper/lower face abscissa mismatch [m] | 0 | 0 | 0 | **2.77556e-17** |
| Native COD samples in 0.04–0.30 a0 | 4 | 7 | 9 | **19** |
| Native COD samples in 0.08–0.30 a0 | 3 | 6 | 8 | **16** |
| Native COD samples in 0.12–0.30 a0 | 2 | 4 | 6 | **13** |

All preflight gates returned true, including the gate for twelve native COD samples in the widest window. More strongly, **each** of the three independently predeclared fitting windows now qualifies for **both** linear (>=8) and quadratic (>=12) COD extrapolation. The characteristic tip-edge median is approximately 31% of the original value; the upper/lower face radial mismatch is consistent with coordinate roundoff. Equal upper/lower tip-triangle *counts* do not imply a reflection-paired T3 mesh.

The selected compact preflight was saved successfully at \`verification/step46_mesh_preflight_face0125_tip05.mat\` (not the previously encountered nested \`verification/verification\` directory). No more mesh-size trial runs are necessary before examining a refined FEM response.

**Next planned phase:** [Step 47](STEP47.md) prepares a **single explicitly authorized** refined Stage-II symmetric FEM solution using this selected mesh recipe, verifies its achieved geometric/radial characteristics against the saved preflight and writes a separate refined checkpoint, then runs COD **only**. No new FEM solve has been executed by the assistant or claimed as completed.
