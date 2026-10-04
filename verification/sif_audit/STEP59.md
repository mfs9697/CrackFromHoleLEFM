# Step 59: prescribed exact nodal pure-I replay on the reflection-paired Step57 mesh

## Why this is the next forensic control

Step58 established that the **actual solved FEM field** on the exactly reflection-paired Step57 mesh gives roundoff-level Mode-II leakage in both independent diagnostics:

- native COD: maximum fitted `|KII/KI| = 1.83854063e-13`;
- matched 16-point FE-nodal-`q` EDI: `KII/KI = -4.686526172e-14`.

This strongly implicates asymmetric discrete representation as the origin of the earlier symmetric leakage. However, Steps45b/50 had shown a second numerical phenomenon on **unpaired** meshes: even a known exact leading Williams field with `KI=1, KII=0`, sampled exactly at T6 nodes and then interpolated inside the EDI calculation, produced artificial signed Mode-II ratios of order `1e-5`:

- original Step45 mesh: approximately `-3.88687e-5`;
- refined Step47 mesh: approximately `+3.37937e-5`.

Direct analytical evaluation at Gauss points reduced those values dramatically, showing that T6 interpolation of the singular prescribed field can itself create mesh-dependent leakage.

Step59 asks whether **exact reflection pairing of the T6 topology removes that synthetic nodal-interpolation leakage as well**.

## Fixed experiment

`main_step59_paired_exact_pureI_nodal.m` loads the completed Step58 `O45` result, which references the already solved Step57 checkpoint. It requires:

- Step57 stage metadata;
- exact Step56 candidate reuse;
- `T3=2578`, `T6=5427`;
- passing Step57 COD symmetry result;
- the completed Step58 actual-FEM EDI result at `r_outer/a0=0.65`;
- the same historical absolute `r_inner=0.00086750243077 m`.

It then calls the previously audited prescribed-field replay routine `main_step45_coarse_exact_pureI_replay` for **exactly one domain**. Despite its historical function name, that routine operates on the symmetric checkpoint supplied through the `O45` structure and uses the saved mesh exactly as stored.

The prescribed field is the leading Williams displacement with `KI=1`, `KII=0`, sampled at the actual T6 nodes. The routine first verifies exact native crack-face COD normalization and then evaluates one 16-point FE-nodal-`q` interaction EDI on the **same Step58 annulus**.

## Interpretation

If the exact-nodal synthetic ratio also collapses to floating-point scale on Step57, the evidence chain becomes unusually strong:

1. unpaired symmetric meshes produced nonzero actual-FEM COD and EDI Mode-II leakage;
2. those meshes also produced `~1e-5` exact-nodal synthetic pure-I EDI leakage, with sign changing under remeshing;
3. exact reflection pairing drove the **actual solved field** COD and EDI leakage to roundoff;
4. Step59 would show whether the **known pure-I nodal interpolation control** also becomes roundoff-small on that same paired topology.

Such a result would strongly support mesh-pairing symmetry as the mechanism that restores cancellation of interpolation-induced Mode-II contamination in the symmetric benchmark.

Even then, this remains a **known-field numerical control**. It does not prove the tiny `KII` in the physically asymmetric crack is exact, and its value must never be subtracted as a correction.

## Local run

After pulling `sif-asymmetric-mesh-audit`:

```matlab
addpath(genpath(pwd));

R59 = main_step59_paired_exact_pureI_nodal();

disp(R59.table);
```

This performs **no FEM solve, no mesh generation, no remeshing and no radius sweep**. It executes one prescribed-field EDI replay on the already saved paired mesh and exact Step58 annulus.

Return the complete console output and stop after this result.
