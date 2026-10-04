# Step 58: one matched-annulus EDI on the solved reflection-paired Step57 field

## Purpose

Step57 established a striking result in the **actual solved symmetric FEM field**: once the mesh was constructed as an exactly reflection-paired T3/T6 discretization, the native crack-face COD Mode-II leakage collapsed to floating-point scale. The maximum fitted `|KII/KI|` over all six predeclared COD fits was `1.83854063e-13`, compared with `6.27056963e-4` on the refined but non-reflection-paired Step47 mesh.

That result strongly implicates mesh-reflection topology as the dominant source of the previous **COD** leakage. It does **not** answer whether the existing interaction-EDI extractor retains its own residual on the paired mesh. Earlier prescribed-field controls showed that FE-nodal interpolation can contaminate tiny EDI Mode-II signals even when the underlying analytical field is pure Mode I.

Step58 therefore asks one narrowly defined question:

> On the same solved Step57 reflection-paired FEM field, what does the unchanged 16-point FE-nodal-`q` EDI extractor return on the exact same physical annulus previously used for Steps45, 47, and 48?

## Fixed experimental domain

The driver reads the completed Step48 compact comparison and reuses its exact absolute annulus:

- `r_inner = 0.00086750243077 m`
- `r_outer = 0.0026 m = 0.65 a0`

For Step57, the measured median tip edge was `hTip = 0.000143632036 m`, hence `2 hTip ≈ 0.000287264072 m`, safely below the fixed Step48 inner radius. The same annulus is therefore admissible without changing any domain parameter.

## What the driver does

`main_step58_reflection_paired_matched_edi.m`:

1. requires the completed Step48 comparison file, the solved Step57 checkpoint, and the saved Step57 COD-only result;
2. verifies that the checkpoint is exactly the reflection-paired Step57 symmetric control (`T3=2578`, `T6=5427`, `Npoly=240`, zero crack angle, exact saved candidate reused);
3. requires the Step57 COD gate to have passed at near-roundoff level;
4. checks that Step48 defines one common absolute annulus and that it remains admissible for the Step57 tip scale;
5. calls the existing audited `main_step45_symmetric_field_leakage` postprocessor with:
   - `RunEDI=true`,
   - exactly one `ROuterOverA0=0.65`,
   - the exact Step48 `InnerRadius`,
   - the existing unchanged 16-point FE-nodal weight-function implementation;
6. verifies that exactly one EDI result was produced on the requested annulus;
7. prints a three-row comparison:
   - Original Step45 actual FEM,
   - Refined Step47 actual FEM,
   - Reflection-paired Step57 actual FEM.

The driver also reports `|q57|/|q47|` and `|q57|/|q45|`. These are descriptive symmetry-control ratios only. They are **not** corrections, uncertainty factors, or transferable bounds for the asymmetric crack calculation.

## Important scope

Step58 performs **no FEM solve, no mesh generation, no remeshing, no stiffness assembly, and no radius sweep**. It evaluates one EDI domain on an already solved field.

If Step57 EDI also collapses toward numerical roundoff, then both independent actual-field diagnostics—COD and EDI—will support the conclusion that reflection-paired topology removes the symmetric Mode-II residual. If EDI remains appreciably larger than COD, that would isolate a residual tied to the EDI extraction/interpolation rather than the solved displacement symmetry itself.

Either outcome remains a symmetric numerical-control result. It must not be subtracted from, or directly converted into an error bar for, the separate asymmetric physical `KII` signal.

## Local run

After pulling `sif-asymmetric-mesh-audit` through GitHub Desktop:

```matlab
addpath(genpath(pwd));

R58 = main_step58_reflection_paired_matched_edi();

disp(R58.comparison);
```

Return the full console output. Stop after this single-domain EDI result; do not start an EDI radius sweep.
