# Proposed compact figure set

Figures are redrawn as vector PGFPlots directly from the audited snapshot.
They do not alter the scientific records or run the numerical solver.
Existing plot identifiers below name their scientific source, not a claim
that every original plot image is currently present on disk.

| Manuscript figure | Source plot/routine and stored data | Scientific purpose | Placement |
|---|---|---|---|
| 1, panels a–b | `13_trajectory_late_detail`, `plot_final_clean_run_additional_results.m`; `State.vertices`, CP2 `C`, stage-I summary | Hole/path context plus a true-scale late detail. Solid P0–P23; orange dashed P23–P24, labeled unsolved. Mark P17, the interpolated ratio zero, and P23. | Main |
| 2, panels a–c | `09_KI_KII_panels` and original `02_theta`; additional/original plotting routines; accepted state rows | Align increasing KI, signed KII, and cumulative absolute direction. Show that changing turn sign does not mean the tip's y coordinate increases. | Main |
| 3, panels a–b | `10_late_mode_mixity_MTS`, supported by `05_mode_mixity`/`03_delta_theta`; accepted rows P15–P23 | Relate the mode-mixity reversal to the MTS-turn sign reversal. Mark the interpolated zero without inserting a new physical point. | Main |
| 4, four panels | `11_COD_turn_sensitivity`, `12_COD_mode_mixity_sensitivity`; all 176 stored `R.fitTable` rows | Compare all four windows and both degrees. Plot COD minus EDI in the late region so small but important differences are visible. EDI is the production reference, not assumed ground truth. | Main |
| S1, four panels | Original `07_numerical_quality`, `plot_final_clean_run_results.m`; accepted solver rows and qualification summaries | Minimum triangle angle, neighboring size ratio, iterations, and true residual normalized by its gate. Verify acceptance without interpreting diagnostics as path-error bounds. | Supplementary/appendix |

## Selection rationale

The 13 original views are not reproduced as 13 manuscript figures. The
trajectory views are combined; the absolute-angle plot is incorporated
with the aligned intensity histories; the late mode-mixity/turn panels
provide the mechanism; and the two COD sensitivity views become one
difference figure. Numerical quality stays outside the central physical
results. The original selected-window COD comparison is replaced with all
eight definitions to avoid selecting the apparently best fit.

The same shared axes and physical units are used within aligned panels.
The two spatial panels have equal x/y scales independently; their zoom
levels differ. The late crack shape is not vertically exaggerated.
All SIF/MTS/quality physical curves stop at P23. The turn predicted at P23
has index 24 but belongs to the P23 field; it does not supply P24 SIFs.

## Before final figure preparation

- Revisit the four-figure selection after increment/boundary checks; do not
  add every diagnostic merely because it exists.
- Retain vector export and review at the eventual journal column width.
- Preserve the qualified-unsolved line style and interpolation wording.
- A Stage-I boundary-stress figure may be added to verification material
  from the restored exact Stage-I source if initiation-site selection needs
  more visual support; it is not required to explain the present main result.
