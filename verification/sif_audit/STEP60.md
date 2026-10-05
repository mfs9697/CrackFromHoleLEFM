# Step 60: visualize the real saved Step38 asymmetric mesh

## Purpose

Step60 is a **read-only topology/geometry visualization** of the actual investigator-local Step38 checkpoint:

`verification/step38_tip_refined_solved.mat`

It exists to replace schematic reasoning with a direct view of the real asymmetric mesh that produced the Step38 SIF result.

The driver does **not** reconstruct geometry. It loads the exact saved T3/T6 coordinates/connectivity and the saved crack topology from the checkpoint.

## What Step60 shows

`main_step60_visualize_step38_mesh.m` creates two PNG figures.

The four-panel overview contains:

1. **Complete saved Step38 T3 mesh** with all free T3 boundaries.
2. **Actual retained-hole/crack neighborhood** in global coordinates.
3. **Crack-tip topology in the local crack frame**, with the historical fixed EDI inner radius and all three outer radii:
   - `r_inner = 0.8 mm`;
   - `r_outer = 4.0, 5.2, 6.4 mm` (`r_outer/a0 = 0.50, 0.65, 0.80`).
4. **Exact FE-nodal-q element support** for the primary `r_outer/a0=0.65` domain.

A second, larger PNG isolates the primary `q`-support view for close inspection.

The crack-face markers are taken from the stored topological node sets. Upper and lower crack faces have coincident collapsed coordinates but distinct node IDs away from the shared tip.

## Exact q-support definition

The highlighted support is not based on centroid inclusion in a geometric annulus.

Step60 reproduces only the **element-participation logic** used by `SIF_LEFM_interaction_EDI.m` for:

- `WeightFunction='fe_nodal'`;
- `QuadratureRule=16`;
- the same radial nodal `q`;
- the same T6 interpolation of `grad q`;
- the same `norm(qgrad)>1e-14` Gauss-point participation condition;
- the same outer-radius quick element test.

It evaluates no displacement field, auxiliary fields, interaction density, or integral. Consequently the highlighted elements are the actual T6 elements in which the production EDI can have nonzero FE-nodal-`q` gradient support.

This distinction matters because T6 elements crossing the nominal inner/outer circles can contribute even when their centroids do not lie strictly between the circles.

## Strong provenance checks

The default file is accepted only if it matches the documented Step38 checkpoint:

- `a0 = 8 mm`;
- crack mouth approximately `[0.19998887220, -0.020817033905] m`;
- `20,164` T3 nodes;
- `39,441` T3 triangles;
- `79,769` T6 nodes;
- measured tip-edge median `5.4024650785e-5 m`;
- fixed EDI inner radius `0.8 mm`;
- stored outer ratios `[0.50,0.65,0.80]`.

A different checkpoint is rejected rather than silently plotted as “Step38.”

## Safety / computational scope

Step60 intentionally loads only:

- `mesh`;
- `crack`;
- `baseline`;
- `actualTip`;
- `a0`.

It does **not load the stored displacement vector `U`**.

It contains:

- no call to `solve_cracked_LEFM`;
- no call to `generateMesh`;
- no call to `SIF_LEFM_interaction_EDI`;
- no stiffness assembly;
- no COD extraction;
- no geometry modification;
- no remeshing.

The only T6 quadrature work is the inexpensive evaluation of `grad q` required to identify support elements.

## Outputs

Default output prefix:

`verification/step60_step38_real_mesh`

The driver saves:

- `step60_step38_real_mesh_overview.png`
- `step60_step38_real_mesh_primary_q_support.png`
- `step60_step38_real_mesh_small_data.mat`

The compact MAT file contains only summary/topology information and support-element IDs; it does not duplicate the mesh or displacement field.

The reported support table gives, for each historical outer radius:

- number of T6 elements with nonzero FE-nodal-`q` gradient support;
- minimum and maximum radial distance of support nodes;
- maximum support-element centroid radius;
- minimum support-node distance to the physical crack mouth.

These measurements are descriptive geometry data, not SIF uncertainty estimates.

## Local run

After pulling `sif-asymmetric-mesh-audit`:

```matlab
addpath(genpath(pwd));

O60 = main_step60_visualize_step38_mesh();

disp(O60.summary);
disp(O60.supportTable);
disp(O60.files);
```

Keep both PNG windows/files. The first scientific question after viewing them is whether a locally reflection-paired region can cover the complete primary `q)-support while preserving the actual retained-hole/mouth geometry and allowing a conforming transition to the unchanged exterior mesh.

Do **not** create a new mesh or run another FEM solve merely from Step60. First inspect the real topology and support extent.
