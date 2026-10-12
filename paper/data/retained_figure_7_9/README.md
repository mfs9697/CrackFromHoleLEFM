# Immutable vector sources for the Figures 7--9 layout pass

These are byte-identical copies of the accepted publication PDFs present at
repository commit `53b3d1d` before the layout-normalization pass. `profiles.json`
pins their complete SHA256 identities and the PDF command spans containing
the plotted graphics and legends. A changed source is rejected before output.

| Copy | Original repository asset |
|---|---|
| `m1_trajectory.pdf` | `paper/figures/m1_mesh_sensitivity/trajectory.pdf` |
| `m1_vertical_deviation.pdf` | `paper/figures/m1_mesh_sensitivity/vertical_deviation.pdf` |
| `m1_direction_deviation.pdf` | `paper/figures/m1_mesh_sensitivity/direction_deviation.pdf` |
| `m1_mode_mixity.pdf` | `paper/figures/m1_mesh_sensitivity/mode_mixity.pdf` |
| `increment_trajectory.pdf` | `paper/figures/increment_sensitivity/trajectory.pdf` |

The last asset is the already verified wide graphical export of the accepted
4/2/1-mm study, with its earlier provenance in
`paper/figures/increment_sensitivity/retained_trajectory_provenance.json`.
The original M1 full history and complete native increment vertex histories
remain unavailable. These copies are **graphical evidence, not replacements
for numerical histories**.

The wrapper embeds the original vector drawing commands. It creates shared
blank MATLAB axes decorations and relocates the original encoded legend
glyphs; it neither digitizes the paths nor synthesizes physical states.
Comparison graphics are resized with the plotting rectangle, as in an axes
resize. Spatial graphics are transformed uniformly using the retained circular
hole as a graphical calibration, preserving a round hole and the existing
physical coordinate scale. The common spatial viewport is 135--300 mm by
-55--15 mm; the increment panel gains blank viewport space rather than new
trajectory data.

Do not overwrite these source copies with normalized output: repeated redraws
must start from these original PDFs, never recursively wrap an already wrapped
publication panel.
