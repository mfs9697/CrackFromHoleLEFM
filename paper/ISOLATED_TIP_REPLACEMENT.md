# Minimal isolated-tip replacement

The calculation archive `isolated_tip_resolution_study_20261009T200424169` is authoritative; the resumed `20261009T163513088` archive corroborates it exactly in checked mechanical records. Neither archive is modified. There are eight fixed-reference-geometry calculations (P17/P21/P22/P23 at core scales 2 and 0.5) and an independently seeded/propagated core=2, exterior=1 P1-P23 history.

The paper varies `CoreScale` while `ExteriorScale=1` holds the requested exterior sizing law fixed. Conforming exterior remeshing is permitted: identical exterior connectivity is not claimed. The accepted h0 reference is unchanged. Fixed-geometry and independent-path zero-crossing estimates are kept separate and remain interpolation diagnostics, not solved states or micrometre-accurate physical locations.

## Before/after numerical replacements

| Independent-path metric | Historical coupled | Corrected isolated |
| --- | ---: | ---: |
| Maximum absolute dy, micrometres | 0.12257 | 0.12760 |
| Maximum corresponding-tip distance, micrometres | 0.12295 | 0.12800 |
| Maximum absolute-direction difference, millidegrees | 0.25475 | 0.24539 |
| Maximum absolute mode-mixity difference | 6.72e-7 | 2.65e-7 |
| Interpolated crossing, mm | 84.145910853 | 84.146276379 |
| Crossing difference from unchanged reference, micrometres | -0.459 | -0.094 |

Fixed-geometry crossings replace 84.140247 / 84.146370 / 84.149515 mm with 84.140272 / 84.146370 / 84.149410 mm for core scales 2 / 1 / 0.5. Successive interpolation differences are 6.10 and 3.04 micrometres, giving the retained diagnostic order approximately 1.00 rather than 0.96. At fixed P17, the fine-core KI change is +0.00134%, q changes by +1.05e-6 and the MTS turn changes by -1.20e-4 degrees.

## Portable reproduction

From the repository root in MATLAB:

```matlab
addpath(genpath(pwd));
plot_tip2h0_vs_reference_publication();
```

This uses `paper/data/isolated_tip_resolution.json`, its byte-provenance manifest, and the unchanged committed exact reference snapshot. Investigator-local archives are not required. `isolated_tip_fixed.csv` and `isolated_tip_states.csv` accompany the canonical JSON. Historical `tip2h0_states.csv` is preserved as historical evidence and is never a fallback input to this figure.

The plot writes the four existing panel filenames, the plotted-data CSV and `isolated_tip_resolution_metrics.json`. Fonts, panel dimensions, palette, markers and axes are preserved. To re-export publication records from available investigator archives, call `export_isolated_tip_publication_data`; this reads/validates existing accepted results and writes portable data only.

Validation:

```matlab
test_isolated_tip_publication();
```

```text
python paper/verify_isolated_tip_correction.py
```

The numerical before/after values are machine-checkable from the metrics JSON. The manuscript's characteristic reference table, M1 row, increment tables and COD evidence are unchanged. The general synthetic-error statement is now explicitly scoped to reference meshes, matching its unchanged reference evidence macro.

See [FIGURE_CONSISTENCY_AUDIT.md](FIGURE_CONSISTENCY_AUDIT.md) for the source, units, verification status, action and limitation of every figure/subfigure. All figure/table labels, section order and figure inclusion order are preserved. No new experiment, model change, tolerance change or branch merge is part of this correction.
