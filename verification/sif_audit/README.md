# SIF audit harness

This folder is the verification layer for comparing the historical
mirror-based circular J/mode-separation extractor with the interaction
equivalent-domain integral (EDI) extractor.

The central rule is simple: **solve the FEM problem once, then pass the
identical mesh and displacement vector to both SIF extractors.**

## First control case

`cfg_crack_path_two_leg_control.m` reproduces the current non-perforated
Crack-Path geometry/material/mesh scale as a two-leg *elastic* crack
control. The second leg is traction free here; this is intentionally not a
CZM solve.

Run the first gate directly:

```matlab
addpath(genpath(pwd));
R = main_step1_same_field_compare();
```

or configure the control explicitly:

```matlab
C = cfg_crack_path_two_leg_control(2.0);
R = run_crack_path_old_vs_edi(C);
```

The output contains:

- `R.solution`: the single FEM field used by both methods;
- `R.old`: `SIF_LEFM_circle2_debug` result and diagnostics;
- `R.edi`: `SIF_LEFM_interaction_EDI` result and diagnostics;
- `R.difference`: signed method-to-method differences.

## Scientific status

The interaction EDI implementation is still a prototype. Its normalization,
mode-II sign, and auxiliary-field derivatives must be independently verified
before EDI is used as the reference method.

The older published two-segment benchmark is a separate reproduction gate.
This control harness does not claim that reproduction yet.
