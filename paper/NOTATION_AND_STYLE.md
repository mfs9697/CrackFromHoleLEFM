# Scientific notation and editorial style — Crack Trajectory

Status: **Batch 1 working convention** (2026-10-11). This is an editorial guide, not a source of numerical evidence. The authoritative numerical results remain in the accepted saved records and the manuscript's cited verification reports. The approved Figures 1–9 and their graphical provenance are not changed by this guide.

## Principles

1. Use italic symbols for physical quantities and running mathematical indices. Use **upright (roman) descriptive subscripts**, e.g., \(r_{\mathrm{in}}\), \(r_{\mathrm{out}}\), \(h_{\mathrm{tip}}\), \(a_{\mathrm{LS}}\), and \(K_{\mathrm{II}}\).
2. Define **every symbol at its first substantive use**, in the introducing sentence or immediately following its defining equation. State the physical meaning, sign convention, and units where applicable.
3. Introduce an abbreviation or mathematical alias only when it is used repeatedly or substantially simplifies the equations. Avoid defining several symbols for the same quantity.
4. Keep manuscript notation distinct from MATLAB configuration-field names. Mathematical typography and descriptive subscripts do not require changes to code variables or numerical archives.
5. Use US English. Favor connected academic exposition over procedural audit language. Place numerical qualifications near the claims they qualify, rather than repeating generic disclaimers.
6. State whether a comparison holds **physical geometry fixed**, **independently propagates an alternative discretization**, or **changes the increment together with mesh and extraction scales**. These are different questions.
7. Reserve **local-symmetry point** for a solved state where \(K_{\mathrm{II}}=0\); in the current evidence, \(a_{\mathrm{LS}}\) is only a **linearly interpolated sign-change estimate**.
8. State clearly that the COD and EDI SIF estimates are different postprocessing methods applied to the same FEM displacement field; they are not independent physical solutions.

## Geometry and propagation

| Symbol | Meaning / convention |
| --- | --- |
| \(2A\), \(B\) | Full plate width and height in the manuscript; respectively 300 and 200 mm in the reference problem |
| \(R\) | Hole radius |
| \((x_c,y_c)\) | Hole-center coordinates in the midheight computational frame; the approved schematic measures the ordinate from the lower edge |
| \(\phi_*\) | Stress-selected initiation angle, counterclockwise from the positive global \(x\)-axis |
| \(P_0\), \(P_k\) | Hole-boundary initiation point and current tip after \(k\) segments |
| \(\Gamma_k\) | Complete existing crack polyline at step \(k\) |
| \(\Delta a\), \(a_k=k\Delta a\) | Prescribed straight-segment length and accumulated crack length |
| \(\bm n_{\mathrm{mat}},\bm t\) | Frozen outward material-side normal at \(P_0\) and its counterclockwise transverse vector |
| \(\theta_k\), \(\Delta\theta_{k+1}\) | Absolute segment angle in frozen frame and next incremental MTS turn; positive is counterclockwise |
| \(\bm e_1^{(k)},\bm e_2^{(k)}\) | Local last-segment forward direction and counterclockwise normal at the tip |
| \(\alpha_{\mathrm{MTS}}\) | MTS turning-angle function; **not** a second independently evolving crack angle |
| \(\kappa(a)\) | Signed curvature of a *hypothetical smooth limiting trajectory*, not the classical curvature of a polygon |
| \(a_{\mathrm{LS}}\) | Linearly interpolated crack length where \(K_{\mathrm{II}}/K_{\mathrm{I}}\) changes sign; not a solved state |

The current figure uses the accepted geometry labels \(x_c,y_c\). Preserve its artwork. The MATLAB configuration describes full width with \`C.A=0.300\` m and half-height with \`C.B=0.100\` m; the manuscript uses \(2A=0.300\) m and \(B=0.200\) m. These conventions must not be conflated.

## Elasticity and fracture mechanics

| Symbol | Meaning |
| --- | --- |
| \(E,\nu\) | Young's modulus and Poisson's ratio |
| \(\lambda,\mu\) | First and second Lamé constants |
| \(\sigma\), \(\sigma_c\) | Applied tensile traction magnitude and stress-based initiation threshold |
| \(\sigma_{tt}\), \(\phi_*\) | Tangential stress around the uncracked hole and its selected maximum-tension angle |
| \(\Lambda_{\mathrm{ini}}\) | Initiation load multiplier (distinct from \(\lambda\)) |
| \(K_{\mathrm{I}}^{(k)},K_{\mathrm{II}}^{(k)}\) | Signed opening- and sliding-mode stress-intensity factors at \(P_k\) |
| \(q_{\mathrm{K}}^{(k)}=K_{\mathrm{II}}^{(k)}/K_{\mathrm{I}}^{(k)}\) | Signed mode-mixity ratio, defined for \(K_{\mathrm{I}}^{(k)}>0\) |
| \(I^{(m)}\) | Mode-\(m\) interaction integral |
| \(w\) | Scalar radial weight used for the equivalent domain integral; not to be confused with \(q_{\mathrm K}\) |
| \(E'=E/(1-\nu^2)\) | Plane-strain effective modulus |
| EDI, COD | Equivalent domain integral and crack-opening displacement; first expand in prose before abbreviation |

It is acceptable to write \(K_{\mathrm{II}}/K_{\mathrm{I}}\) at the point of definition and when physically helpful. Otherwise use \(q_{\mathrm K}\) in repeated derivations rather than alternating between redundant notations.

## Discretization

| Symbol | Meaning |
| --- | --- |
| \(h_{\mathrm{tip}}\), \(h_0\) | Characteristic structured-tip spacing and production reference spacing |
| \(r_{\mathrm{in}},r_{\mathrm{out}}\) | Inner and outer radii of the interaction-integral annulus |
| \(r_{\mathrm{core}}\) | Radius of the structured near-tip mesh |
| \(h_{\mathrm{far}}\) | Requested far-field exterior mesh size cap |
| \(L_{\mathrm{tr}}\) | Grading transition length |
| \(2h_0,h_0,h_0/2\) | Physically evaluated structured-tip resolution family on prescribed crack geometries |
| M1 | Coarser exterior-sizing law with the production structured core retained |
| \(\delta y,\delta r,\delta\theta\) | Differences of matched-tip position and absolute segment direction between propagated trajectories; define their direction at first use |

Descriptive scale factors such as \(\eta_c\) and \(\eta_e\) should appear only if they meaningfully simplify the relevant numerical family definition; otherwise state the actual changed and fixed quantities.

## Required scientific distinctions

- **One-step centered-hole angular response:** the fitted slope \(F'(0)\approx-0.11415\) is a restricted angular response of a fixed family of two-segment geometries, not proof of attraction or stability of the full evolving crack.
- **Mode-mixity/curvature relation:** the MTS small-mixity expansion is analytic; the relation \(q_{\mathrm K}/\Delta a\to-\kappa/2\) is *conditional* on an appropriately smooth limiting path, convergent discrete turning densities, and extraction accuracy on the relevant scale.
- **Three-level increment family:** \(4/2/1\) mm calculations independently propagate their paths while absolute near-tip, EDI, transition, and exterior scales change together. They provide a coupled increment/discretization refinement study, not a pure increment experiment or formal convergence theorem.
- **Characteristic events:** the positive mode-mixity maximum, a bracketed sign change, and reversed MTS turning direction are observed in the accepted discrete solutions. The interpolated zero is not an independent physical solution.
- **Near free boundary:** shrinking ligament and rising \(K_{\mathrm I}\) are recorded; attributing redirection uniquely to boundary proximity requires a separate control.
- **Accepted states:** the reference physical history ends at \(P_{23}\). The qualified \(P_{24}\) geometry has no accepted physical solution.
- **First segment:** the initiation site is selected from the uncracked-hole stress; the initial outward-normal segment is prescribed, not selected by MTS. Only one active crack is considered.

## Prose conventions

Write each subsection as a scientific argument: the mechanical question, the mathematical/numerical method essential to answering it, the evidence, then the narrow conclusion. Reduce duplicated descriptions across Results, Discussion, and Conclusions. Prefer *indicates*, *is consistent with*, *is preserved under the tested change*, or the exact measured quantity over categorical claims of *proof*, *mesh independence*, *stability*, or *convergence*. Use meaningful numerical precision; historical full-precision evidence is retained separately. Keep approved figure/subfigure environments inline in \`paper/main.tex\`.
