# Notation and scientific presentation review

The machine-readable register `notation_register.csv` covers 47 symbol
families/operators, with meanings, units, first source occurrence, subsequent
line references and distinctions. It excludes implementation identifiers from
the article and groups indexed variants without treating every state index as
a new physical variable.

## Implemented clarifications

- Plate notation is width `2A`, full height `B`. The software uses `C.A=2A`
  and `C.B=B/2`; no FEM geometry is changed. The schematic's lower-left
  center ordinate maps to the computational ordinate `y_c-B/2`.
- `theta_k` is absolute direction relative to the frozen initiation normal;
  `Theta_k=phi_*+theta_k` is the global direction. `Delta theta_(k+1)` is
  the local MTS update, whereas lowercase `delta theta` denotes a difference
  between trajectories at a matching physical crack length. This distinction
  also applies to the sensitivity tables and plot labels.
- Applied traction has one unsigned scalar name `sigma` in the sketch and
  equations. Bold/indexed stress denotes a tensor/component. The tensile
  initiation threshold `sigma_c` is not a fracture-growth criterion.
- The Lamé coefficient `lambda` has stress units; the positive initiation
  multiplier `lambda_ini` is dimensionless. Full notation and units make the
  distinction explicit without inventing another unnecessary coefficient.
- `q_K=K_II/K_I` is defined before use and differs from the radial weight `q`.
  `d Omega` denotes area measure, avoiding confusion with geometric half-width
  `A`. Apparent COD intensities and signed displacement jumps are identified.
- `eta_c=h_tip/h_0` and `eta_e` describe independent core/exterior sizing.
  `eta_e=1` fixes the requested exterior law, not connectivity. In the software,
  these remain `CoreScale` and `ExteriorScale`; their option names are unchanged.
- `a_k=k Delta a` is accumulated crack length, not horizontal coordinate.
  Lowercase `delta r` is the distance between matching tips. The order
  diagnostic uses a successive-difference measure `d`, not an asserted true
  discretization error. The COD constant `kappa` is not a curvature variable.

## Article editing and precision

The abstract is 190 words, self-contained, with no numbers, citations, figure
references or implementation names. The Introduction preserves existing
literature citations and states the contribution as a reproducible assessment
of recursive direction selection with established LEFM, interaction-integral
and MTS methods. Its final paragraph explains the scientific role of each
subsequent section.

Operational audit details are removed from the article while essential method
and acceptance information remains. The exact MTS root-selection safeguards
are unchanged in `kink_angle_LEFM_MTS.m`: the tensile stationary root, bounded
search, small-mode handling and compressive safeguard remain software details,
not a new rule. Model limits remain direction selection with prescribed
quasi-static increments, rather than growth kinetics, arrest or branching.

`update_presentation_numbers.py` changes displayed values only. Full-precision
JSON, MAT, CSV records and numerical verification are not overwritten. Geometry
and initiation displays use meaningful precision; state-table columns use
consistent rounding. Additional digits in mesh-comparison crossings resolve
numerical differences only and are explicitly not physical location accuracy.
The exact mesh calibration constant is retained as a reproducibility parameter.

## Scientific interpretation and proposals

No material change to the scientific interpretation is made. The supported
findings remain a mode-mixity maximum and sign reversal, opposite-signed MTS
turning, contraction of path differences under tested increment refinement,
and small mesh-family effects. Neither a rigorous continuum error bound,
micrometre-accurate local-symmetry location, nor an isolated free-boundary cause
is asserted.

For eventual journal preparation, the strongest contribution is the separation
of recursive trajectory sensitivity from extraction/mesh sensitivity, not a
new LEFM or MTS algorithm. A stronger causal boundary claim or formal convergence
claim would require new evidence and is **not implemented**. The prose follows
conventional research-article structure and US English; journal-specific class,
author information and submission requirements remain editorial decisions.

## Reproduction limitations

The raw M1 history and original hierarchy coordinates remain absent. Complete
native increment trajectories are also unavailable in this checkout, although
accepted plotting records survive for its three comparison panels. Verified
graphics are retained where regeneration would require reconstruction.
Their exact font/interpreter targets cannot be certified without raw sources.
This limitation is explicit in the graphical workflow and final audit.
