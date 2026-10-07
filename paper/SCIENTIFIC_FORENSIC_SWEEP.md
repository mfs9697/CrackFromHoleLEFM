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
