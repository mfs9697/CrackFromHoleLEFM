# Google Drive snapshot preserved 2026-10-06

This directory preserves the two files from the historical Google Drive
`CrackFromHoleLEFM` snapshot (created/synchronized 2026-07-30) that did not
have the same repository path in the audited GitHub tree.

- `sample_hole_boundary_node_stress.m` — historical Stage-I helper for
  extrapolating integration-point stresses to actual T6 hole-boundary nodes
  and averaging adjacent-element contributions.
- `main_stage1_centered_hole_stress_check.m` — old root-level copy of the
  Stage-I diagnostic driver. The maintained GitHub workflow lives under
  `Stage1_Hole_Diagnostics/`; this copy is retained only for provenance.

These files are archival and are not part of the current production path.
The formerly external `BN_local.m` dependency is preserved separately at
the repository root as the exact historical helper.
