# Step 54: prescribed polygon reflection audit — geometry ONLY

## Why this small check comes next

The investigator's Step53 common-grid comparison removed both sampling confounders from Step52: all **42/42 COD point pairs** and **49/49 matched EDI-annulus point pairs** were successfully evaluated on both saved FEM meshes at identical physical locations, with one shared actual opening denominator. The affine T6 checks were within ~3e-17 m and direct pure Mode-I symmetry checks were zero. Measured actual FEM displacement parity residuals dropped with refinement by approximately **74%/58%** in the COD region (even ux / gauge-removed odd uy), and **43%/12%** in the matched EDI annulus, but the measured EDI Mode-II ratio had nevertheless changed sign and grown in magnitude. These results **do not establish convergence of a tiny physical Mode-II SIF**, and an additional FEM solve is not justified solely by improved displacement-parity RMS.

Step49 also established that neither original nor refined FEM mesh was exactly reflection-paired **off the crack faces**. Before constructing a fully reflection-paired new mesh, we must check the **input polygon itself**. A mirror-paired triangulation of an *asymmetric* polygon would change the mathematical problem and invalidate a comparison with the already solved checkpoints.

## Single zero-mesh calculation

The new `main_step54_geometry_symmetry_preflight.m` builds **only the geometry-description polygon** used by Steps45–47, with the existing repository functions `cfg_centered_half_domain` and `build_domain_centered_half_pencil`. It verifies that the original and refined saved checkpoint crack endpoints coincide and that the rebuilt `Npoly=240`, `a0=0.004 m` crack coordinates match the previously solved fields.

It then reflects the **entire outer boundary vertex list** and all its consecutive **boundary-edge segments** about `y=0`, testing whether each reflected vertex and complete **undirected edge** has an actual match to a stated tolerance of `1e-12 m`. Checking segments as well as vertices matters: different polygons could have the same vertex set but different boundary connectivity. The check also separately verifies matched quarter-hole upper/lower arc vertex counts, reflection pairing of both sharp-pencil faces and mouth endpoints, and a tip lying on the symmetry axis.

**Limitations:** passing this test would establish geometric reflection of the **prescribed polygon** reconstructed from the current repository's known baseline recipe. It is not a verification of unsaved historical PDE geometry preprocessing, of the mesh topology, of numerical traction parity, or of EDI/SIF accuracy. Failure would prohibit silently replacing the polygon with a mirrored variant; first investigate the prescribed discretized hole and appended-pencil construction.

The script generates **no PDE mesh, stiffness matrix, solution, EDI or COD output** and saves only compact geometry diagnostics. It leaves both previously saved full FEM checkpoints unchanged.

## One local MATLAB run

After pulling branch `sif-asymmetric-mesh-audit` through GitHub Desktop:

```matlab
addpath(genpath(pwd));

O54 = main_step54_geometry_symmetry_preflight();

disp(O54.summary);
```

The result is saved in `verification/step54_geometry_symmetry_small_data.mat`. Please return the full printed summary, especially **max reflected boundary vertex error**, **max full-edge error**, **sharp-pencil mouth/face pairing**, and the final Boolean gate.

**Decision gate:** if the complete prescribed polygon and appendix passes, proceed to *design a small mesh-only upper-half construction and mirror-assembly prototype* using that exact existing polygon. Do not authorize an additional FEM solve until any proposed new reflection-paired triangulation preserves the verified polygon, has correct duplicated crack faces/shared intact-ligament nodes, passes nondegenerate-element tests, and is evaluated separately for COD sampling and EDI-domain coverage.


## Investigator's completed Step 54 results

The investigator executed the geometry-only preflight successfully using the **existing Step45 and Step47 solved checkpoint crack endpoints** and the repository's original `Npoly=240`, `a0=0.004 m` polygonal geometry recipe. The requested reflection tolerance was **`1.000e-12 m`**.

| Independently measured geometry diagnostic | Actual MATLAB output |
| --- | ---: |
| Prescribed outer polygon vertices | 127 |
| Upper/lower circular-arc vertices | 61 / 61 |
| Maximum reflected polygon vertex error | 2.7972e-17 m |
| Maximum reflected whole-boundary-edge error | 2.7972e-17 m |
| Sharp-pencil mouth reflection error | 0 m |
| Sharp-pencil entire-face reflection error | 0 m |
| Sharp-tip distance from horizontal reflection axis | 0 m |
| Unpaired prescribed exterior vertices | 0 |
| Unpaired prescribed boundary edges | 0 |
| `prescribedPolygonReflectionPass` | **true** |

This **closes the prescribed-source-polygon reflection question** at the stated tolerance. It does **not** establish an exactly symmetric PDE-generated mesh or symmetric computed FEM solution. In particular, Step49 had observed **zero** exact reflected off-face node matches in both tested T6 meshes, and neither tip-adjacent triangle fan had complete mirror pairing. A later mirror-assembled mesh must retain the passed source polygon, intentionally duplicate crack-face nodes **behind** the tip and share intact-ligament nodes **ahead** of the tip.

The next Step55 geometry-only gate, described in [STEP55.md](STEP55.md), extracts an original-vertex upper material domain and creates **one artificial bookkeeping cut from the original sharp tip to the right plate boundary**. Reflecting all its non-cut exterior edges must reconstruct the **entire original verified source polygon** (with only a harmless subdivision of the plate's right vertical edge at `y=0`). This checks the upper-half domain description before any new mesh generation. No additional FEM solve has been authorized.
