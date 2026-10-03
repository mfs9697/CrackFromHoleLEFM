# Step 55: derive a mirror-assemblable upper-half boundary without changing the verified polygon

## Input established by the investigator

Step54 was executed on the previously solved centered right-half symmetry control: `Npoly=240`, `a0=0.004 m`, crack horizontal at `y=0`. It verified the **entire prescribed 127-vertex polygon** and every boundary edge under reflection about `y=0` at tolerance `1e-12 m`. The actual reported maximum mirrored vertex and whole-edge errors were both `2.7972e-17 m`. The upper/lower circular-arc vertex counts were **61/61**; the crack mouth and both sharp-pencil face descriptions reflected with **zero reported error**, and there were **zero unpaired vertices or boundary edges**. Its saved compact result is `verification/step54_geometry_symmetry_small_data.mat`.

This establishes that the **prescribed polygon** is symmetric, not that MATLAB PDE meshing or the computed FEM solution is symmetric. The two existing original/refined meshes lack interior reflection pairing (Step49); although common-grid displacement parity improved after refinement (Step53), actual EDI Mode II remained sign-changing and unresolved.

## Mesh-topology requirement before the next solve

A future mirror-assembled **T3 mesh** must satisfy two fundamentally different conditions on the horizontal crack line:

- **Crack faces behind the sharp tip:** upper and lower boundary nodes may occupy the same spatial coordinates **after pencil-face collapse**, but must retain **distinct node numbers and displacement degrees of freedom**. Coincident coordinates are intentional.
- **Intact ligament ahead of the sharp tip:** upper and reflected lower triangles must meet on a **single shared node set** along `[x_tip,A]`, including a **single shared tip node**. Sharing this artificial cut is essential: it is not a physical free surface.
- Every interior upper node/element must have a reflected lower counterpart, with consistent T3 element orientation and no degenerate elements. The complete original exterior boundary must be retained. **None of these mesh checks is implied merely by a symmetric boundary.**

Before PDE mesh generation, we perform one even narrower **polygon-only split**. This gives the later mesh assembler a checked upper-half source without inventing a different hole or crack.

## Precisely one Step55 geometry-only calculation

`main_step55_upper_half_boundary_preflight.m` loads the investigator's passed Step54 output and reconstructs the **same exact** original `build_domain_centered_half_pencil` polygon, without calling `generateMesh` or creating a PDE model.

It constructs the **upper material domain** by retaining the original uninterrupted exterior path from `(x_sym,+B)` through the **original upper circular arc**, upper pencil mouth and sharp upper pencil edge to the original tip. It adds exactly one **artificial straight horizontal segment** `[x_tip,A]` on `y=0` (the future shared intact ligament), then closes the upper domain using the original right/top rectangle exterior edges. The artificial cut is a bookkeeping boundary for upper-domain triangulation, **not** a new physical crack face.

As a stringent source-geometry check, the driver reflects *all non-cut upper polygon edges* and requires one-to-one complete **undirected segment** matching against the original full polygon's entire prescribed edge set. For this comparison only, it permits the original right-plate vertical edge to be subdivided at `(A,0)`; this changes an edge description, not the boundary curve. It also tests:

- The two and only two upper-domain polygon vertices on `y=0` are the original crack tip and the right-plate midpoint; the artificial cut lies exactly between them.
- The original sharp upper pencil flank from the original upper mouth to the original tip appears **unchanged**.
- The upper polygon has positive oriented area, no degenerate edges, all vertices on or above `y=0` and matching checkpoint geometry.
- No original exterior edge is lost, duplicated or replaced by a different segment.

The `upperSplitAndReconstructionPass` result only authorizes **planning a mesh-only upper triangulation and reflected assembly**, not an FEM solution. An expected upper-domain polygon and artificial-cut coordinates are retained in a small `O55` compact MAT report to avoid regenerating or hand-editing the boundary.

## Local run

After pulling `sif-asymmetric-mesh-audit` through GitHub Desktop:

```matlab
addpath(genpath(pwd));
O55 = main_step55_upper_half_boundary_preflight();
disp(O55.summary);
```

The compact data is saved to `verification/step55_upper_half_boundary_small_data.mat`. Please return the complete printed output, especially the reported `nOriginalSplitBoundaryEdges`, `nReconstructedExternalEdges`, `maxExteriorEdgeReconstructionError_m` and final gate.

**No PDE mesh, stiffness matrix, FEM solve, EDI or displacement-based COD analysis is performed.** We will not create a new FEM solution without your separate explicit authorization. A later mesh-only prototype can use the verified `O55.upperHalfPolygon` to create and mirror an upper T3 mesh, then independently test duplicate crack-face topology, shared ligament nodes, exact element reflection and fit-window sampling before any FEM solve.
