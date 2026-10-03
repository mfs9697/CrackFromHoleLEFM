# Step 56: reflect ONE upper triangulation, assemble the crack and ligament (mesh ONLY)

## Verified source geometry

The investigator's successful Step55 geometry-only MATLAB run established:

| Step55 measurement | Observed |
| --- | ---: |
| Original complete polygon | 127 vertices |
| Upper-half source polygon | 65 vertices |
| Source boundary edges after permitted subdivision | 128 |
| Reconstructed exterior segments | 128 |
| Upper material area | 0.014293 m² |
| Shortest upper prescribed polygon edge | 0.00078371 m |
| Artificial horizontal ligament cut | 0.116 m |
| Upper vertices on horizontal axis | Exactly 2: tip and right boundary |
| Maximum complete-edge reconstruction mismatch | 2.7972e-17 m |
| Unmatched original exterior edges | 0 |
| `upperSplitAndReconstructionPass` | true |

These values verify that the upper-half polygon retains the original external boundary, hole arc and original sharp-pencil upper face without a new physical crack. Step55 saved the exact upper polygon at `verification/step55_upper_half_boundary_small_data.mat`.

## One controlled mesh-only assembly experiment

`main_step56_reflected_mesh_preflight.m` loads the **already executed** Step55 report and the existing selected Step47 material/crack checkpoint. It reconstructs the original full 127-vertex polygon only as a consistency check. Then it:

1. Builds one PDE model for the verified **65-vertex upper polygon**. It generates a small temporary coarse upper-only mesh to identify its exact PDE edge labels for the upper pencil face and the **artificial tip-to-right intact-ligament seam**, and the vertex label of the shared sharp tip. It rejects missing/ambiguous labels. The temporary mesh provides geometric IDs only—no FEM solve.
2. Generates **one refined upper T3 mesh** with requested `FaceFactor=0.125` and `TipFactor=0.5` (same *nominal* factors as selected Step46/47). The achieved mesh is new and is **not** assumed to reproduce Step47 connectivity or achieved element sizes. Verifies that every original upper-polygon vertex remains present.
3. Creates the lower T3 mesh by reflecting the **same actual upper T3 coordinates and connectivity**, reversing each reflected triangle orientation. It keeps upper/lower nodes **identical by ID only on the artificial intact ligament ahead of the tip**, including the single shared tip. All other reflected nodes get distinct IDs.
4. Confirms that original upper and mirrored lower **sharp-pencil-face T3 nodes** are distinct behind the common tip. Only these pencil-face nodes are projected onto the existing original horizontal crack midline, using the original mouth and tip endpoints. On the final collapsed mesh the crack faces intentionally have coincident coordinates but different node IDs.
5. Validates that all reflected **T3 triangles**, and every vertex/midside of the corresponding **T6 triangles**, are mirror-paired. Every intact seam T3 edge must have exactly one upper and one lower incident triangle and its T6 midside node **must be shared**. Every consecutive upper/lower crack-face T6 edge must have **distinct midside node IDs**. Checks all collapsed T3 signed areas, the smallest internal angle and achieved tip-fan counts.
6. Runs the previously audited `native_COD_audit` with **zero displacements exclusively for geometry sampling**. It reports actual upper/lower T6 crack-face counts, radial mismatch, and population in the original three COD windows. It never interprets the zero field as a physical COD solution, and performs **no FEM solve or EDI**.

The driver distinguishes two gates: `meshReadyForReview` checks the structural reflection/topology and nondegenerate-area prerequisites; the stronger `readyForOneReflectedFEMProposal` additionally requires **at least twelve actual native COD samples in each of the three predeclared windows** and achieved tip refinement relative to original Step45. Neither gate is authorization to compute a new FEM solution. Inspect the actual mesh quality and all listed metrics before discussing any subsequent solve.

On a structurally passing result, it saves the **exact generated T3 coordinates/connectivity**, the explicitly split crack-node topology, and source-polygon provenance to `verification/step56_reflected_mesh_only_candidate_T3.mat`. A separate compact diagnostic file is `verification/step56_reflected_mesh_only_small_data.mat`. PDE models, stiffness matrices and field displacements are never saved or generated.

## Local MATLAB run

After pulling branch `sif-asymmetric-mesh-audit` through GitHub Desktop:

```matlab
addpath(genpath(pwd));

O56 = main_step56_reflected_mesh_preflight();

disp(O56.summary);
disp(O56.gates);
```

Return the complete output, particularly the identified upper-only PDE edge/vertex labels, the achieved min triangle area and angle, complete T3/T6 mirror errors, interface topology gates, actual tip-edge median and **all three native COD window counts**. MATLAB numerical execution is not yet verified; stop at the first error and report it rather than guessing new geometry IDs.

**Research interpretation:** even a fully mirrored mesh with a valid collapsed slit is a *new mesh* of the same physical problem, not evidence that existing FEM/EID Mode-II residuals are calibrated or that the separate much-finer asymmetric tiny-Mode-II signal is verified. If the mesh-only gates pass, the next decision is an explicit, separate authorization for **at most one** newly solved symmetric control on this exact candidate.
