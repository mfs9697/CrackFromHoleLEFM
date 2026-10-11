# Scientific notation and editorial style — Crack Trajectory

Status: **Editorial convention through Batch 6** (2026-10-11). This is an editorial guide, not a source of numerical evidence. The authoritative numerical results remain in the accepted saved records and the manuscript's cited verification reports. The approved Figures 1–9 and their graphical provenance are not changed by this guide.

## Principles

1. Use italic symbols for physical quantities and mathematical component or running indices, including \(\sigma_{tt}\) (the tangential--tangential stress component) and \(\sigma_{12}\). Use **upright (roman) descriptive subscripts**, e.g., \(r_{\mathrm{in}}\), \(r_{\mathrm{out}}\), \(h_{\mathrm{tip}}\), \(a_{\mathrm{LS}}\), \(x_{\mathrm c}\), and \(K_{\mathrm{II}}\).
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
| \((x_{\mathrm{c}},y_{\mathrm{c}})\) | Hole-center coordinates in the schematic, measured from the lower-left corner; the computational center has ordinate \(y_{\mathrm{c}}-B/2\) |
| \(\phi_*\) | Stress-selected initiation angle, counterclockwise from the positive global \(x\)-axis |
| \(P_0\), \(P_k\) | Hole-boundary initiation point and current tip after \(k\) segments |
| \(\Gamma_k\) | Complete existing crack polyline at step \(k\) |
| \(\Delta a\), \(a_k=k\Delta a\) | Prescribed straight-segment length and accumulated crack length |
| \(\bm n_{\mathrm{mat}},\bm t\) | Frozen outward material-side normal at \(P_0\) and its counterclockwise transverse vector |
| \(\theta_k\), \(\Delta\theta_{k+1}\) | Absolute segment angle in frozen frame and next incremental MTS turn; positive is counterclockwise |
| \(\bm e_1^{(k)},\bm e_2^{(k)}\) | Local last-segment forward direction and counterclockwise normal at the tip |
| \(\alpha_{\mathrm{MTS}}\) | MTS turning-angle function; **not** a second independently evolving crack angle |
| \(\mathrm{d}\theta/\mathrm{d}a\) | Signed curvature of a *hypothetical smooth limiting trajectory*, not pointwise curvature of a polygon |
| \(a_{\mathrm{LS}}\) | Linearly interpolated crack length where \(K_{\mathrm{II}}/K_{\mathrm{I}}\) changes sign; not a solved state |

The approved figure labels the center coordinates as \(x_{\mathrm{c}},y_{\mathrm{c}}\) in a lower-left frame; preserve its artwork. The computational origin is at the left-edge midpoint, and the hole-center ordinate there is \(y_{\mathrm{c}}-B/2\), not a separate variable. Thus the schematic center \((170,80)\) mm corresponds to \((170,-20)\) mm computationally. The MATLAB configuration describes full width with \`C.A=0.300\` m and half-height with \`C.B=0.100\` m; the manuscript uses \(2A=0.300\) m and \(B=0.200\) m. These conventions must not be conflated.

## Elasticity and fracture mechanics

| Symbol | Meaning |
| --- | --- |
| \(E,\nu\) | Young's modulus and Poisson's ratio |
| \(\lambda,\mu\) | First and second Lamé constants |
| \(\sigma\), \(\sigma_{\mathrm{ini}}\) | Applied tensile traction and prescribed tensile-stress threshold for initiation |
| \(\sigma_{tt}\), \(\phi_*\) | Tangential--tangential stress component around the uncracked hole (italic component indices) and its selected maximum-tension angle |
| \(\Lambda_{\mathrm{ini}}\) | Initiation load multiplier (distinct from \(\lambda\)) |
| \(K_{\mathrm{I}}^{(k)},K_{\mathrm{II}}^{(k)}\) | Signed opening- and sliding-mode stress-intensity factors at \(P_k\) |
| \(q_{\mathrm{K}}^{(k)}\) | Dimensionless signed mode-mixity ratio, defined by the numbered Eq. \`eq:mode-mixity\` for \(K_{\mathrm{I}}^{(k)}>0\) |
| \(I^{(m)}\) | Mode-\(m\) interaction integral |
| \(S(\alpha),T(\alpha)\) | Angular amplitudes of circumferential normal and radial--circumferential shear stresses, respectively |
| \(\kappa=3-4\nu\) | Dimensionless plane-strain Kolosov constant in COD extraction |
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
| Native crack-face node pair | Corresponding existing upper and lower T6 nodes at the same distance behind the current crack tip; COD uses their displacement jump without inserting additional sampling locations |
| \(\delta y,\delta r,\delta\theta\) | Differences of matched-tip position and absolute segment direction between propagated trajectories; define their direction at first use |

The manuscript uses actual tip spacings and the requested exterior sizing law rather than additional scale-factor aliases. This reduces notation without changing the implementation variables or saved numerical records.

## Numerical verification hierarchy

The manuscript distinguishes three experimental designs. **Fixed-geometry** comparisons evaluate SIF extraction and prescribed mesh changes without advancing the crack. **Independently propagated** exterior (M1) and structured-tip (\(2h_0\)) trajectories solve their own \(P_1\) fields and recursively recompute all subsequent directions; no accepted reference vertices after \(P_1\) are imposed. The **4/2/1-mm family** varies \(\Delta a\), the structured-core size, and the interaction-domain and exterior scales together. It is a coupled increment/discretization refinement study, not isolated crack-increment convergence.

The centered-hole control checks reflection symmetry and the restricted one-step map \(F(\theta_2)=\theta_3\); its fitted slope does not establish full evolving-trajectory stability. The table is intentionally placed after the defining equation in the manuscript source.

**COD sampling:** The four overlapping intervals in normalized rear-tip distance \(r/\Delta a\) are \([0.04,0.20]\), \([0.04,0.30]\), \([0.08,0.30]\), and \([0.12,0.30]\). Counts refer to native crack-face node pairs in that order; overlapping counts are not additive. The existing qualification requires at least 12 pairs in each interval. The illustrative \(4h_0\) core has 10/15/12/9 pairs and is excluded from the physical SIF/MTS comparisons. Keep the authoritative counts in \`paper/main.tex\` and the historical verification records; do not regenerate or interpolate missing samples.

## Required scientific distinctions

- **One-step centered-hole angular response:** the fitted slope \(F'(0)\approx-0.11415\) is a restricted angular response of a fixed family of two-segment geometries, not proof of attraction or stability of the full evolving crack.
- **Mode-mixity/curvature relation:** the MTS small-mixity expansion is analytic; the relation \(q_{\mathrm K}/\Delta a\to-\frac12\,\mathrm{d}\theta/\mathrm{d}a\) is *conditional* on an appropriately smooth limiting path, convergent discrete turning densities, and extraction accuracy on the relevant scale.
- **Three-level increment family:** \(4/2/1\) mm calculations independently propagate their paths while absolute near-tip, EDI, transition, and exterior scales change together. They provide a coupled increment/discretization refinement study, not a pure increment experiment or formal convergence theorem.
- **Characteristic events:** the positive mode-mixity maximum, a bracketed sign change, and reversed MTS turning direction are observed in the accepted discrete solutions. The interpolated zero is not an independent physical solution.
- **Near free boundary:** shrinking ligament and rising \(K_{\mathrm I}\) are recorded; attributing redirection uniquely to boundary proximity requires a separate control.
- **Accepted states:** the reference physical history ends at \(P_{23}\). The qualified \(P_{24}\) geometry has no accepted physical solution.
- **First segment and initiation stress:** the initiation site is selected from the uncracked-hole tangential stress. The adopted threshold \(\sigma_{\mathrm{ini}}=300\) MPa determines the initiation load multiplier \(\Lambda_{\mathrm{ini}}\); it is not a CZM cohesive strength. The initial outward-normal segment is prescribed, not selected by MTS. No cohesive traction--separation law or subsequent strength-based propagation criterion is introduced; only one active crack is considered.

## Introduction, Abstract, and bibliography conventions

The Introduction should move from the distinction between stress-selected initiation and finite-crack propagation, through established LEFM direction criteria and numerical path methods, to the need to assess recursively accumulated extraction/discretization differences. State the specific scientific question and contribution without claiming a new fracture criterion or general continuum convergence.

The Abstract should describe the method, characteristic mode-mixity sign reversal, and the tested numerical evidence. The 4/2/1-mm calculations vary the crack increment *and* the numerical length scales; do not call them an increment-only experiment. Do not describe COD and EDI as independent elastic solutions, a local-symmetry interpolation as a solved tip, or the nearing free boundary as an isolated causal mechanism.

Citations must match populated `paper/references.bib` entries and retain documented DOI metadata. Use the published **article number**, not an invented page range, where the journal uses article identifiers. For the 2026 Bogdanov--Selivanov IJF reference, Springer identifies **volume 250, article 35** (DOI 10.1007/s10704-026-00930-1). An article number may be recorded in the BibTeX note for compatibility with the current bibliography style. Check all first-use abbreviations and retain US English.

## Discussion and conclusion conventions

Separate the **direction-selection mechanism** (the signed MTS turn determined by current-tip mode mixity) from its **physical cause** (how hole geometry, total crack length, and the approaching free boundary contribute to the stress-intensity factors). The current calculations have not independently isolated free-boundary proximity.

Interpret fixed-geometry extraction, independently propagated mesh perturbations, and coupled 4/2/1-mm refinement only within their respective numerical designs. The latter changes several length scales together and cannot establish an isolated increment-error estimate, a formal convergence order, or a unique smooth limiting curve. The centered-hole angular slope refers only to the prescribed one-step map. The accepted path ends at \(P_{23}\), not at the qualified but unsolved \(P_{24}\) geometry.

Discussion paragraphs should synthesize explanations and limitations rather than reproduce tables of values. Conclusions should identify the method, the observed directional phenomenon, the scope of numerical support, and the essential outstanding questions without repeating the Results section.

## Main-text and caption distinction

Keep a concise methodological statement establishing that alternative mesh families are propagated independently from the common first segment. Do not rehearse each solver update in the Results or in the captions; describe the observed geometrical, angular, and mode-mixity responses instead. Numerical-acceptance and provenance details belong in the numerical methodology or archived verification evidence when they do not affect interpretation.

Figure captions state the plotted quantities, the compared configurations, and qualifications necessary to avoid a false physical inference. Avoid repeating the paragraph introducing the figure or narrating how the code generated it. Table captions identify the comparison conditions and the status of interpolated quantities; the table body retains the authoritative numerical values.

As an editorial starting point, graphical groups in Figures 2--9 and supplementary Figure S1 are set within `0.8\\textwidth` minipages; Figure 1 already has an included width of `0.78\\textwidth`. Figure and table captions use `\\small` and subfigure captions `\\footnotesize`. The three four-panel figure layouts retain their original panel ordering, with reduced fixed vertical boxes. These are **provisional display dimensions** that require checking a freshly compiled PDF for legibility, float placement and whitespace. Adjust locally if dense axes become too small; never reconstruct or alter the approved vector graphics.

## Prose conventions

Write each subsection as a scientific argument: the mechanical question, the mathematical/numerical method essential to answering it, the evidence, then the narrow conclusion. Reduce duplicated descriptions across Results, Discussion, and Conclusions. Prefer *indicates*, *is consistent with*, *is preserved under the tested change*, or the exact measured quantity over categorical claims of *proof*, *mesh independence*, *stability*, or *convergence*. Use meaningful numerical precision; historical full-precision evidence is retained separately. Keep approved figure/subfigure environments inline in \`paper/main.tex\`.
