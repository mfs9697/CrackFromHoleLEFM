# Scientific forensic sweep

This file records scientific checks of the pilot manuscript against the
implemented mechanics and the preserved numerical evidence.  It is an audit
record, not manuscript prose.  No numerical solver is rerun and no scientific
source is modified.

## Batch 1 — governing problem, initiation, MTS, and path recurrence

Status: **PASS WITH ITEMS FOR AUTHOR REVIEW**

### 1. Geometry, material, and loading — PASS

The manuscript problem statement agrees with
`cfg_first_segment_asymmetric.m` and the physical solver:

- plate: `A=0.30 m`, `B=0.10 m`, so
  `[0,A] x [-B,B]`;
- circular hole: center `(0.17,-0.02) m`, radius `0.03 m`,
  480 polygon segments;
- isotropic material: `E=210e3 MPa`, `nu=0.30`;
- plane strain: `ps=1`;
- normalized remote-y traction: `sig0=1`.

For the production crack solves, the top and bottom signed resultants are
explicitly gated to `+A` and `-A` at unit traction.  No load is applied to
the vertical sides, hole boundary, or crack faces.  The manuscript statement
that those boundaries are traction free is consistent with the implementation.

The minimal anchoring is also reproduced correctly.  The code fixes

- `ux=0`, `uy=0` at `(0,-B)`;
- `uy=0` at `(A,-B)`.

This removes rigid translation and rotation without imposing a symmetry
condition on the asymmetric full-domain problem.

**Notation item, not a mechanics error:** the manuscript uses `lambda` both
for the Lamé constant in the constitutive law and `lambda_ini` for the
initiation load factor.  These are distinguishable but unnecessarily easy to
confuse.  Rename the Lamé constant during prose/equation editing.

### 2. Stage-I initiation rule and frozen frame — PASS, with one physical-scope issue

The manuscript accurately describes the production Stage-I procedure:

1. solve the hole-only problem at unit remote-y loading;
2. recover nodal T6 stresses;
3. interpolate them through the actual containing T6 element at three
   material-side offsets;
4. linearly extrapolate the stress components to the hole boundary;
5. find the maximum positive extrapolated tangential stress;
6. refine the angular location by the accepted local quadratic fit.

The implemented selection is
`sigPos=max(sig_tt_eff,0)`, followed by the discrete maximum and the local
quadratic refinement.  The reported point and frame are therefore consistent
with the manuscript equations.

The frame is also correct:

`n_mat=[cos(phi_star),sin(phi_star)]`

is radial from the hole center into the material, and

`t_hat=[-sin(phi_star),cos(phi_star)]`

is the counterclockwise tangent.  The first segment is then prescribed as
`P1=P0+Delta a*n_mat` with `Delta a=4 mm` and `theta_1=0`.

The manuscript correctly says that the earlier first-angle sweep did not
select the production `theta_1`.

#### Author-review item A1 — initiation strength is not needed for the path-direction result

The configuration contains `sig_c=300 MPa`, and the code computes

`lambda_ini=sig_c/sigma_tt,max(unit)`.

This is implemented correctly, but the subsequent trajectory uses unit
traction and depends on SIF ratios/directions, not on `lambda_ini`.
Therefore the 300-MPa threshold is not part of the central directional
prediction.  Before final drafting, either:

- justify this strength value physically and retain the initiation load, or
- present it as a normalization/example and avoid making it appear to be a
  calibrated fracture-initiation parameter essential to the trajectory.

#### Author-review item A2 — near-competing Stage-I tensile maximum

The accepted Stage-I record has a secondary tensile maximum only about
1.336% below the selected primary maximum and separated by about 176 degrees.
The single-crack model is mathematically well defined because the primary
maximum is unique under the accepted gate, but the small gap is physically
important.  At the primary threshold, the competing site is close to the
same tensile level.

The manuscript already states that only one active crack is retained.  For a
journal paper this needs to remain an explicit scope limitation: the model
does not determine whether a competing crack would initiate after the first
crack changes the stress field.  This is not a numerical error in the present
single-crack trajectory.

### 3. Current-tip SIF frame — PASS

The manuscript defines the current frame from the last segment,
not from the mouth-to-tip chord.  This is exactly how
`SIF_LEFM_interaction_EDI.m` constructs

`e1=(P_k-P_{k-1})/|P_k-P_{k-1}|`,
`e2=(-e1_y,e1_x)`.

This convention is also used by the polyline-safe COD wrapper.  The
manuscript's signed definitions

`KI = lim sqrt(2*pi*r) sigma_22(r,0)`

and

`KII = lim sqrt(2*pi*r) sigma_12(r,0)`

are consistent with the implemented local Cartesian auxiliary fields and the
stored sign history.

No contradictory mouth-to-tip SIF frame was found in the production path
algorithm.

### 4. Implemented MTS rule — EQUATIONS PASS; criterion wording needs precision

The manuscript reproduces the implemented functions exactly:

`S(alpha)=KI*cos(alpha/2)^3 - 3*KII*sin(alpha/2)*cos(alpha/2)^2`

and

`T(alpha)=KI*sin(alpha/2)*cos(alpha/2)^2
          +KII*cos(alpha/2)*(1-3*sin(alpha/2)^2)`.

Direct differentiation gives

`dS/dalpha = -(3/2) T`,

so solving `T(alpha)=0` is exactly the stationary circumferential-stress
condition used by the code.

The manuscript also correctly records the actual numerical root procedure:

- zero turn when `|KII|<1e-12`;
- small-mixity initial estimate `-2*KII/KI`;
- clipped search interval `[-80 deg,80 deg]`;
- initial `+/-10 deg` bracket and expansion to the full interval;
- the nearly pure-II initial guess;
- the opposite-angle safeguard when the computed circumferential amplitude
  is compressive.

No production angular sweep is performed.

#### Author-review item A3 — distinguish "stationary tensile branch" from a general maximization algorithm

The implementation solves the shear-zero/stationarity equation near the
expected tensile branch and applies a compressive-root safeguard.  It does
not enumerate every stationary root and explicitly compare all local maxima
of `S(alpha)`.

For the accepted history, `KI>0` and the turns stay on the small tensile
branch, so this does not invalidate any stored result.  But the methods
section should not imply that the code performs a global maximization of
`sigma_theta theta` over all angles.  The present wording
"tensile stationary direction by solving T=0" is appropriately faithful and
should be preserved or sharpened rather than replaced by a stronger claim.

### 5. Small-mode-mixity asymptotics — PASS

Expanding the implemented `T(alpha)=0` relation for small `alpha` gives

`alpha ~= -2 KII/KI`

in radians.  The manuscript uses this only to explain the sign of the turn
and explicitly says it does not replace the numerical root calculation.
That is correct.

### 6. Incremental recurrence — PASS

The manuscript recurrence agrees with `run_incremental_crack_path.m`:

`Delta theta_(k+1)=alpha_MTS(KI_k,KII_k)`;

`theta_(k+1)=theta_k+Delta theta_(k+1)`;

`P_(k+1)=P_k+Delta a [cos(theta_(k+1)) n_mat
                       +sin(theta_(k+1)) t_hat]`.

The first segment is prescribed.  The accepted P1 SIF pair seeds the MTS
prediction of `theta_2`.  Thereafter each accepted fixed geometry is
qualified, solved once, postprocessed, and at most one new 4-mm segment is
appended.

The production driver does not optimize `theta_k` by scanning candidate
directions.  The manuscript's "one physical solve per fixed current
geometry" statement is correct.

### 7. Direction prediction versus propagation law — PASS and essential scope statement

The current computation prescribes the finite crack advance
`Delta a=4 mm`.  MTS determines its direction but there is no toughness,
energy-balance, load-increment, or crack-growth threshold deciding whether
that advance should occur.

The manuscript correctly calls the study a direction-selection calculation
under normalized loading rather than a load-versus-time growth history.
This distinction is scientifically essential and should remain prominent.

### 8. First-batch verdict

No contradiction was found between the manuscript's governing problem,
Stage-I selection, frozen-frame convention, signed MTS equations, or
incremental recurrence and the production implementation.

The first batch therefore has **no blocking mechanics error**.

Items to carry into later manuscript editing:

1. avoid the dual use of `lambda`;
2. decide how much physical significance to assign to `sig_c=300 MPa`;
3. retain the competing-initiation-site limitation caused by the 1.336% peak
   gap;
4. describe the implemented MTS selector as a tensile stationary-root
   procedure on the accepted branch, not as an exhaustive global angular
   maximizer.

Next forensic batch: interaction/domain-integral formulation, auxiliary
fields, COD extraction, mesh/extraction verification, and what the existing
controls do and do not establish.


## Batch 2 — interaction EDI, auxiliary fields, COD, mesh controls, and solver evidence

Status: **PASS WITH IMPORTANT INTERPRETATION LIMITS**

This batch checks the manuscript equations and claims against
`SIF_LEFM_interaction_EDI.m`,
`native_COD_polyline_audit.m`,
`verification/sif_audit/native_COD_audit.m`,
`qualify_incremental_crack_candidate.m`, and
`solve_incremental_crack_tip.m`.

### 9. Interaction-integral density — PASS

The manuscript reproduces the implemented interaction density correctly.
With the physical field labeled 1 and the auxiliary field labeled 2, the
code evaluates

`A_j = -W^(1,2) delta_(1j)
       + sigma^(1)_(ij) u^(2)_(i,1)
       + sigma^(2)_(ij) u^(1)_(i,1)`

and integrates `A_j q_,j` over the finite-element domain.

The implemented mutual-energy term is

`W^(1,2)=sigma^(2):epsilon^(1)`.

For the same isotropic linear-elastic constitutive tensor this equals
`sigma^(1):epsilon^(2)` by reciprocity, so the apparently asymmetric
coding of the scalar term is not a mechanics inconsistency.

The engineering-shear convention is also handled correctly in the code:

`sigma:epsilon = sigma11*epsilon11 + sigma22*epsilon22
                  + sigma12*gamma12`.

No missing factor of two was found in that contraction.

### 10. EDI normalization — PASS

For isotropic LEFM the implementation uses

`I^(1,2) = (2/E') [KI^(1) KI^(2) + KII^(1) KII^(2)]`

and therefore, for a pure auxiliary field of amplitude
`K_aux`,

`K_m = E' I^(m)/(2 K_aux)`.

For plane strain the code and manuscript both use

`E'=E/(1-nu^2)`.

This matches the actual implementation.  The historical factor-of-two
normalization error documented in the SIF audit is not present in the
current production extractor.

**Notation item B1:** `K_aux=1` is numerically one but physically has SIF
units.  The final manuscript should avoid language that might make it look
dimensionless.

### 11. Weight function and extraction support — PASS, and stronger than the prose currently emphasizes

For the production option `WeightFunction='fe_nodal'`, nodal values are

- `q=1` for `r<=r_i`;
- linearly decreasing nodally for `r_i<r<r_o`;
- `q=0` for `r>=r_o`.

The gradient used in the EDI is obtained from the same T6 interpolation as
the displacement field.  The production values are

`r_i=0.4 mm = 0.10 Delta a`,
`r_o=2.6 mm = 0.65 Delta a`.

A useful geometric fact is hard-gated but is only implicit in the pilot:
the previous crack vertex is one full increment, 4 mm, behind the current
tip.  Since `r_o=2.6 mm < 4 mm`, the complete primary EDI support remains
inside the straight current segment and cannot reach the previous kink.

Likewise the paired core radius is 3 mm, still smaller than the 4-mm segment.
The qualifier requires the EDI support to lie wholly inside that untouched
paired core.  The support fingerprint is 11316 T3 elements for every
qualified P2--P24 geometry.

This is an important positive verification point: the current-tip
interaction domain never straddles a historical kink.

### 12. Free-boundary clearance of the local extraction region — PASS

The late path does not fail because the EDI annulus itself touches the right
boundary.

Recorded minimum physical-boundary clearances are approximately:

- P17: 32.071 mm;
- P21: 16.128 mm;
- P22: 12.145 mm;
- P23: 8.160 mm;
- qualified P24: 4.171 mm.

All exceed the 3-mm paired-core radius and the 2.6-mm EDI outer radius.
Therefore the right free boundary influences the physical solution
globally, but it does not geometrically cut through the local extraction
annulus in any accepted state or in the qualified P24 geometry.

This distinction should be retained when discussing "boundary influence":
the late change is not a trivial contour/domain-intersection artifact.

### 13. Auxiliary Williams fields — PASS, with a dependency caveat

The auxiliary stresses are analytical leading-order isotropic crack-tip
fields in the current-tip frame.  Auxiliary displacement derivatives are
computed by centered finite differences with

`h=max(1e-7 m,1e-5 r)`.

For the production annulus, the first term controls the step throughout
because `r<=2.6e-3 m`.

The physical interaction density uses:

- analytical auxiliary stress;
- centered-difference auxiliary `u_,1`;
- FEM T6 physical displacement gradients and stresses.

The auxiliary strain calculated from finite-difference displacements is not
needed in the production interaction density itself; it exists mainly for
diagnostic/analytic-replay work.  The manuscript wording is compatible with
this, but a final derivation should not imply that a separately computed
auxiliary strain enters the production mutual-energy term.

#### Author-review item B2 — synthetic Williams tests are not an independent theory implementation

The prescribed Williams replay is a very strong topology/interpolation/
extractor test, but its imposed displacement field uses the same standard
Williams sign convention as the EDI auxiliary fields.  It should therefore
not be described as an independent theoretical derivation of the crack-tip
fields.  Its proper evidential role is exactly what the pilot mostly says:
modal recovery, leakage, superposition, and mesh/extractor qualification.

### 14. Synthetic qualification — PASS and correctly scoped

For every qualified path state the code prescribes three nodal fields:

- pure I: `(KI,KII)=(1,0)`;
- pure II: `(0,1)`;
- tiny mixed: `(1,1e-4)`.

No physical boundary-value problem is solved for these fields.  They are
replayed through the same T6 interpolation, current-tip topology, FE-nodal
q, 16-point quadrature, and EDI extractor.

The gates verify pure-mode recovery, cross leakage, the 2-by-2 modal
recovery matrix, tiny-mixed recovery, superposition, and the same 11316
element support.

The manuscript correctly states that these controls qualify extraction
capability and topology and **do not establish the error of the physical FEM
solution**.  This limitation is scientifically essential.

### 15. 16-point quadrature — PASS, citation still required

The production EDI explicitly requests the 16-point triangular rule.  The
implementation identifies it as the Dunavant degree-8 rule.  The manuscript's
"16-point triangular rule" is correct.

A verified quadrature citation is still required before submission.  This is
a bibliography gap, not a numerical inconsistency.

### 16. COD formula — PASS

For plane strain, the code uses

`kappa=3-4 nu`

and converts the upper-minus-lower crack-face displacement jump through

`KI_app  = mu/(kappa+1) sqrt(2 pi/r) Delta u_2`,
`KII_app = mu/(kappa+1) sqrt(2 pi/r) Delta u_1`.

This follows directly from the same leading Williams displacement
convention at the upper and lower crack faces.  The manuscript equation is
therefore consistent with the implementation and with the sign convention
used by the EDI.

### 17. COD sampling geometry — PASS, wording can be more exact

The polyline wrapper intentionally replaces only the extraction descriptor
by the last two path vertices.  Thus COD distances and components are
measured in the last-segment frame, never along the mouth-to-tip chord.

The four windows end at at most `0.30 Delta a = 1.2 mm`, well inside the
current 4-mm segment and also inside the 3-mm paired core.

The historical COD routine classifies native T6 upper/lower face nodes,
sorts them by tip distance, and evaluates the lower-face displacement on
the upper-face abscissae with shape-preserving cubic interpolation.  In the
qualified production meshes the upper/lower abscissae are additionally
gated to match to numerical tolerance, so this interpolation does not hide
a mismatched face grid.

#### Author-review item B3 — "independent COD" should mean independent postprocessor, not independent solution

COD and EDI use the same physical displacement vector.  They differ in
postprocessing principle and finite-distance behavior, but they share all
physical FEM discretization error.

The manuscript already says this explicitly.  A final subsection title such
as "Alternative COD-based SIF extraction" or "Independent postprocessing
route" would be harder to misread than simply "Independent COD
postprocessing."

### 18. What the COD comparison actually establishes — PASS

Across all 176 stored COD fits, every definition gives positive mode mixity
at P21 and negative mode mixity at P22.  Thus the **sign-change bracket**
between those solved states is robust to the recorded COD window/order
choices.

The comparison does not establish the precise zero-crossing length.

At P23 the largest absolute COD--EDI turn difference is about
0.010236 degrees, roughly 1.64% of the EDI turn of 0.624347 degrees.
The four quadratic COD turns lie about 0.30--0.63% below the EDI value.
Near P21 one linear rear-window ratio differs from EDI by 26.3% relatively,
even though the absolute turning difference is only 0.002453 degrees.

The pilot interprets these numbers correctly: relative errors become
ill-conditioned near local symmetry, and the 84.15-mm zero is a linear
visualization estimate between 4-mm-spaced solved states, not a precision
measurement.

### 19. Reflection-paired core — PASS, but do not confuse mesh symmetry with physical symmetry

The 12678-element structured tip core is reflection paired about the current
crack line.  Opposite crack faces remain distinct; the intact ligament
remains shared; T3/T6 coordinate pairing and topology are hard-gated.

No displacement symmetry constraint is imposed on the physical asymmetric
problem.  The paired core is a **numerical anti-leakage mesh design**, not a
symmetry reduction of the physical solution.  Consequently a nonzero
physical KII remains fully admissible.

This distinction is worth making explicit if reviewers could otherwise read
"reflection-paired mesh" as a symmetry boundary condition.

### 20. Mesh-quality gates — PASS as integrity checks, not convergence proof

The qualifier checks positive T3 areas, positive T6 Jacobians, preserved
physical boundaries, seam conformity, no duplicate triangles, maximum edge
incidence, minimum angle, neighboring-size ratio, exact core coordinates,
and EDI containment.

The current acceptance limits include minimum angle 20 degrees and maximum
neighboring size ratio 1.8.  Accepted late states remain inside those gates,
although the size-ratio values approach 1.8 closely.

These checks demonstrate that the intended mesh construction remains valid.
They do **not** demonstrate path convergence.  The pilot makes this
distinction correctly.

### 21. Linear solver evidence — PASS as algebraic accuracy only

For each accepted physical state the current solver:

1. assembles the unclamped symmetric stiffness;
2. checks its symmetry;
3. removes the three minimal anchor DOFs;
4. symmetrizes the free matrix by `(Kff+Kff')/2`;
5. applies `symamd`;
6. uses parameter-free SGS preconditioning;
7. solves with PCG at reported tolerance `1e-10`, maximum 5000 iterations;
8. independently recomputes the free-DOF true relative residual;
9. requires that residual to be at most `5e-10`.

The accepted P2--P23 runs require 2442--3011 iterations and have true
residuals below the stated gate.  The physical field is checkpointed after
the algebraic solve passes and before COD/EDI/MTS postprocessing.

This establishes tight linear-system accuracy.  It does not bound FE
discretization or fracture-parameter error.  Again, the pilot says this
correctly.

### 22. Historical 8-mm convergence audit — useful but nontransferable

The historical straight 8-mm problem has a scientifically closed
three-level paired-mesh audit and provides strong evidence for the EDI/COD
architecture and the removal of mesh-induced mode-II leakage.

It does **not** provide increment-size convergence or mesh convergence of the
present evolving 4-mm polyline path.  The pilot correctly refuses to promote
that historical convergence result into a present-path error bound.

This is conservative and should not be weakened during editing.

### 23. Second-batch verdict

No blocking error was found in the manuscript interaction-integral
normalization, auxiliary-field convention, COD conversion, last-segment
frames, mesh-support containment, synthetic controls, or solver claims.

The most important positive finding is that the late free boundary never
intersects the paired core or EDI annulus.  Hence the observed mode-mixity
reversal cannot be dismissed as an extraction domain physically colliding
with the boundary.

The most important limitations are also clear:

1. EDI-domain sensitivity for the evolving path has not yet been performed;
2. physical mesh refinement at the late tips has not yet been performed;
3. COD is an alternative extraction route from the same field, not an
   independent physical solution;
4. the synthetic Williams replay validates the extractor/topology, not the
   global FEM solution;
5. the precise 84.15-mm local-symmetry location is not converged with respect
   to the 4-mm crack increment.

Next forensic batch: inspect the **physical interpretation of the trajectory
itself** — monotonic KI, the P17 maximum, the P21--P22 sign change, the
turning reversal, the right-boundary attribution, P24 termination, and
whether any causal or predictive claim in the abstract/discussion goes
beyond what the computed sequence actually proves.


## Batch 3 — physical interpretation of the accepted trajectory and terminal state

Status: **PASS, WITH ONE CENTRAL CAUSAL CLAIM DELIBERATELY UNPROVEN**

This batch checks the abstract, Results, Discussion, Conclusions, figure plan,
and preserved state tables against the accepted P1--P23 sequence and the
qualified-unsolved P24 geometry.

### 24. Monotonic opening-mode intensity — PASS

The accepted P1--P23 record has strictly increasing `KI` at unit reference
traction:

`0.36648 -> ... -> 0.78679 (P17) -> 0.96631 (P21)
           -> 1.06672 (P22) -> 1.24436 (P23)`.

The manuscript statement that opening-mode intensity increases throughout
the accepted history is therefore exactly supported by the preserved data.

This is a **geometry-dependent SIF history under the same normalized
traction**, not an experimentally prescribed load history.

#### Author-review item C1 — replace "unloading" by a mechanics-exact statement

Because every geometry is solved at the same unit reference traction and the
MTS direction is invariant under positive load scaling, external unloading
is not a competing process represented in this calculation.

The scientifically precise observation is:

> the mode-mixity reversal occurs while the normalized opening-mode SIF
> continues to increase.

The manuscript's phrase "rather than unloading" is understandable, but it
can be read as if an external load-history mechanism had been compared.
That comparison was not performed.

### 25. P17 maximum — PASS

Both the positive `KII` and the positive ratio `KII/KI` reach their
discrete accepted maxima at P17 (68 mm):

`KII(P17)=2.54836873974682e-3 MPa sqrt(m)`;

`(KII/KI)(P17)=3.23893080298227e-3`.

The MTS turn also reaches its largest negative magnitude there:

`Delta theta_18 = -0.371145045509677 deg`.

After P17, `KI` continues to rise while `KII`, `KII/KI`, and the
negative-turn magnitude all decrease through P21.  The manuscript's stated
sequence is therefore correct.

A useful interpretive point is that the coincidence of the P17 maxima in
`KII` and `KII/KI` is an observed property of this discrete path, not an
identity that must hold generally.

### 26. The MTS turn is exactly a function of mode mixity — IMPORTANT POSITIVE RESULT

For `KI != 0`, dividing the implemented stationarity equation by `KI`
shows that the MTS angle depends on the SIF pair only through

`q=KII/KI`.

Thus, under the present LEFM/MTS model, changes in the magnitude of `KI`
alone do not determine the turn.  The sign and magnitude evolution of the
incremental turn are an algebraic consequence of the changing mode mixity.

This makes the manuscript's mode-mixity interpretation stronger than the
phrase "primarily driven by mode II" might suggest: **within the chosen MTS
law, the local direction selector is controlled by the SIF ratio, not by
their common positive scale**.

At the same time, this does not identify which global geometric feature
causes the ratio itself to evolve.

#### Figure-design consequence C2

The plotted MTS-turn history is not an independent numerical observable
relative to the `KII/KI` history; it is computed deterministically from it.
Showing the two aligned is useful pedagogically for sign conventions, but it
is not two independent pieces of physical evidence.  A final figure caption
should make that relationship clear.

### 27. Local-symmetry bracket — PASS

The accepted solved states give

`q(P21)=+8.13426737346537e-5` at 84 mm,

`q(P22)=-2.14159160711531e-3` at 88 mm.

Since `KI>0`, this brackets `KII=0` and hence the MTS zero-turn state.

The manuscript's linear interpolation of `q` gives

`a_LS = 84.146369911941 mm`.

This is correctly labeled as an interpolation rather than a solved crack
state.

A forensic sensitivity check illustrates why the extra digits are not
physical precision.  If one linearly interpolates `KII` itself instead of
the ratio, the estimate is approximately

`84.133049974822 mm`,

only 0.01332 mm lower but already different in the third decimal place.
Linear interpolation of the MTS turn gives approximately

`84.146371418828 mm`,

essentially the same as ratio interpolation near zero.  All of these
differences are tiny compared with the 4-mm state spacing; they reinforce
the manuscript's decision to report the crossing only as a visualization
estimate.

### 28. Turning reversal and accepted geometry — PASS, but the amount of post-reversal evidence is limited

At P21 the turn is still slightly negative:

`Delta theta_22=-0.00932118365470 deg`.

At P22, after the sign change in mode II,

`Delta theta_23=+0.245405694839291 deg`.

That positive turn generates the segment P22--P23, and P23 is an accepted
physical state.  Therefore the **first geometrical manifestation of the
turning reversal is physically accepted**, not merely extrapolated.

At P23 the physical solve gives a still more negative mode mixity and
predicts

`Delta theta_24=+0.624347038002271 deg`.

That value belongs to the accepted P23 field, but the resulting P23--P24
geometry has no accepted P24 physical SIF solution.

#### Author-review item C3 — avoid implying a long established post-reversal branch

There is one accepted new segment after the turn changes sign
(P22--P23).  The next positive turn is a valid prediction from P23, but its
endpoint is only qualified geometrically.

Thus the paper can robustly claim **onset of turning reversal** and a
continued positive next-turn prediction.  It should be more cautious with
language suggesting that a long return branch or asymptotic reorientation
has already been demonstrated.

### 29. Absolute direction and vertical motion — PASS

The most negative accepted local absolute direction occurs at P22:

`theta_22=-3.64596382938318 deg`.

The accepted P23 segment has the less negative direction

`theta_23=-3.40055813454389 deg`.

Because the frozen material normal itself is inclined downward in global
coordinates by about 1.5606 degrees, both segments remain downward-sloping
globally.  The tip y-coordinate decreases monotonically through P23.

The manuscript correctly warns that a positive incremental turn does not
mean that the crack tip has started moving upward.

### 30. Right-boundary interpretation — PLAUSIBLE AND CAREFULLY WORDED, BUT NOT CAUSALLY ESTABLISHED

The geometrical evidence for increasing proximity to the right boundary is
unambiguous.  The horizontal/right-boundary clearance decreases from about

- 32.07 mm at P17,
- 16.13 mm at P21,
- 12.14 mm at P22,
- 8.16 mm at P23,

while the crack advances toward `x=A`.

The local EDI/core region remains well clear of that boundary, so the
observed reversal is not a direct extraction-domain collision artifact.

However, along this one trajectory several quantities change together:

- distance to the right free boundary decreases;
- distance from the original hole increases;
- total crack length increases;
- the complete crack shape/curvature changes;
- the current-tip orientation changes.

Therefore the present sequence cannot uniquely attribute the evolution of
`KII/KI` to the right free boundary.

The manuscript's current formulation -- "consistent with increasing
free-boundary influence" and "does not prove that the boundary is its sole
cause" -- is scientifically defensible.

#### Author-review item C4 — do not promote "boundary-induced redirection" to the main claim before a control

A stronger phrase such as "boundary-induced redirection" would require the
control already identified in `REMAINING_WORK.md`: move the right boundary
or compare matched fixed crack geometries while controlling what is held
fixed.

A particularly clean experiment would evaluate the **same late accepted
polyline geometries** in plates with progressively more distant right
boundaries.  That would separate the boundary-distance effect from changes
in the crack path itself better than simply recomputing a wholly different
trajectory in a wider plate.

### 31. Hole influence is also a covarying physical mechanism

As the crack approaches the right boundary it is simultaneously moving away
from the circular hole from which it initiated.  The present calculation
does not separately quantify the decay of hole-induced asymmetry versus the
growth of right-boundary influence.

Therefore a future causal discussion should not reduce the alternatives only
to "boundary distance versus accumulated curvature"; changing hole influence
is another geometric contribution that can be separated only by controlled
comparisons.

This does not invalidate the current interpretation, but it broadens the
mechanism that a reviewer may ask about.

### 32. P24 classification — PASS

The preserved scientific archive establishes:

- P24 geometry exists at total path length 96 mm;
- its qualification gates pass;
- its physical-boundary clearance is about 4.171 mm, still larger than the
  3-mm paired core;
- no accepted P24 physical result exists.

The paper correctly stops every physical SIF/mode-mixity curve at P23 and
draws P23--P24 only as a qualified-unsolved extension.

The repository evidence available to the pilot does not preserve an accepted
P24 displacement/SIF result and does not independently preserve the complete
failed residual trace.  Therefore the manuscript is right not to quote a
P24 SIF or reinterpret the strict solver gate after the fact.

#### Author-review item C5 — "linear-solver limit" is stronger than the preserved evidence

The abstract currently refers to the "eventual linear-solver limit."  What
is strictly established by the manuscript evidence layer is that the P24
attempt did not produce an accepted physical state under the prescribed
solver acceptance rule.

Unless the original failure trace is archived with the paper evidence,
prefer "P24 did not satisfy the strict linear-solver acceptance criterion"
to language implying a fundamental or asymptotic solver limit.

### 33. P24 failure is not physical arrest — PASS

No toughness condition, energy-release criterion, propagation threshold,
dynamic instability calculation, or arrest criterion is evaluated at P24.
The absence of an accepted numerical solve therefore contains no evidence
that the physical crack arrests there.

The Discussion and Conclusions explicitly preserve this distinction.

### 34. Predictive scope of the trajectory — PASS, with title-level caution

The method predicts the **direction** of a sequence of imposed 4-mm advances.
It does not predict when an advance occurs or whether a particular external
load is sufficient for the next advance.

Calling the result a crack trajectory is reasonable provided that assumption
is explicit, as it is in the present manuscript.  Claims of a complete
fracture-growth law would not be supported.

This point will matter when choosing the final title and abstract language:
"incremental LEFM crack-trajectory prediction" is acceptable; "crack-growth
prediction" without qualification would be too broad.

### 35. Third-batch verdict

The central physical sequence is real and internally coherent:

`positive q -> P17 maximum -> decay -> P21/P22 zero bracket
 -> negative q -> positive MTS turn -> accepted P22--P23 redirection`.

The strongest defensible scientific statement is:

> In the accepted fixed-increment LEFM sequence, the current-tip mode mixity
> evolves nonmonotonically and changes sign while the normalized mode-I SIF
> continues to increase; because the MTS direction depends on `KII/KI`,
> this sign change reverses the incremental turning sense.

What is **not yet established** is why the mode mixity changes sign.  The
approaching right free boundary is a physically plausible contributor and
the local extraction domain remains cleanly separated from it, but a causal
boundary-effect claim requires a controlled boundary comparison.

The present manuscript is notably disciplined on this distinction.  The
scientific priority before submission is therefore not to rewrite the
interpretation more aggressively, but to perform the control that would
allow the interpretation to become stronger.

Next forensic batch: assess the paper as a potential journal contribution
rather than only as a correct computation -- novelty versus classical MTS
literature, which verification material is essential in the main paper,
what additional calculations most efficiently raise the evidence level, and
what a skeptical fracture-mechanics reviewer is most likely to challenge.
