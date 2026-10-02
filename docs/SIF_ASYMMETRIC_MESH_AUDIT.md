# SIF extraction on asymmetric FEM meshes: audit log

## Purpose

The audit asks a narrow question: how much does the historical
mirror-based mode separation distort (K_I) and (K_{II}) when the
discrete FEM field is not mirror symmetric, and does an interaction
equivalent-domain integral provide a stable alternative on arbitrary
asymmetric meshes?

No attempt will be made to preserve or justify old numerical values.  Old
results are retained as historical/control values until independently
verified.

## Source-code audit

Two repositories are involved:

- **mfs9697/Crack-Path** supplies the established Crack-Path geometry,
  T6 FEM, loading, and the historical SIF workflow.
- **mfs9697/CrackFromHoleLEFM** contains both the historical
  `SIF_LEFM_circle2` family and the prototype
  `SIF_LEFM_interaction_EDI`.

The old `SIF_LEFM_circle2.m` present in both repositories implements the
same conceptual procedure: mirrored contour points (P,Q) are evaluated
independently in the discrete T6 field, symmetric/antisymmetric fields are
formed, separate (J_I,J_{II}) values are integrated, and those are
converted to SIFs.

The key numerical vulnerability is therefore not continuum mode separation
itself. It is the assumption that the two independently interpolated FEM
fields at mirrored locations have the parity properties of the continuum
solution. An asymmetric triangulation need not satisfy that assumption.

## Step 1 implemented on branch `sif-asymmetric-mesh-audit`

The first code-extraction step deliberately leaves both SIF extractors
unchanged.

A verification-only pipeline has been added:

1. `cfg_crack_path_two_leg_control.m` defines a two-leg elastic control
   case using the current non-perforated Crack-Path scales.
2. `build_crack_path_polyline_LEFM_mesh.m` adapts the pencil-channel
   geometry from Crack-Path and collapses *all* crack segments to a
   traction-free polyline while retaining distinct crack-face topology.
3. `solve_crack_path_polyline_field.m` solves one T6 displacement field
   with the Crack-Path loading and rigid-body constraints.
4. `run_crack_path_old_vs_edi.m` sends that exact same FEM field to
   `SIF_LEFM_circle2_debug` and `SIF_LEFM_interaction_EDI`.

The helper `subdivide_last_leg.m`, which was used by existing geometry
code but absent from CrackFromHoleLEFM, has also been restored from
Crack-Path.

Two repository-level issues found during extraction are intentionally kept
outside the production path for now:

- `kinking_LEFM_1leg.m` in CrackFromHoleLEFM calls
  `geom_pencil_1leg.m`, but that geometry file is absent from the target
  repository;
- the existing `mesh_pencil_domain.m` computes its `Hmin` from an
  already subdivided last segment and then divides by `ncoh` again,
  effectively introducing an extra factor of `ncoh` in the local target
  size. The verification builder therefore implements the Crack-Path
  `h_last = L_last/ncoh` rule directly instead of using that routine.

Neither production issue is changed in this PR because the first goal is to
isolate the SIF comparison from unrelated refactoring.

## What this step does *not* establish

This first control case is **not yet the historical published two-segment
benchmark**. The exact mapping of the older publication geometry and its
plane-state convention is still to be closed before a published numerical
baseline is declared reproduced.

Likewise, `SIF_LEFM_interaction_EDI.m` remains a prototype. At this stage
old-vs-EDI differences are method-to-method differences, not EDI-based
errors.

## Verification gates

The next gates are:

- independently verify the EDI normalization and mode-II sign using exact
  analytical crack-tip fields;
- establish a deliberately mirror-symmetric FEM mesh and require old/EDI
  agreement there;
- construct controlled asymmetric mesh families while keeping geometry,
  material, loading, and crack path fixed;
- sweep old contour radius and EDI annulus radii;
- correlate old-method changes with explicit mesh-asymmetry measures;
- only after those gates, decide which extractor is suitable for production
  crack-path calculations.

Incremental crack growth is outside the scope of this audit until these
gates are closed.


## First same-field numerical result

The first local MATLAB run of `main_step1_same_field_compare` produced:

| Method | KI | KII |
|---|---:|---:|
| old mirror/J | 0.33973 | 0.0033961 |
| interaction EDI | 0.68007 | 0.0073326 |

with the same FEM field, old contour radius (8.0\times10^{-4}), and EDI
annulus ([8.0\times10^{-5},8.0\times10^{-4}]).

The mode-I ratio is approximately (0.68007/0.33973 \approx 2.0018). This
is a strong diagnostic signature of a possible factor-of-two normalization
issue in the EDI conversion, but it is not by itself sufficient to alter the
production implementation. The mode-II value is small and correspondingly
more sensitive; its ratio is not used as the normalization diagnostic.

## Step 2: synthetic Williams-field normalization test

The branch now includes:

- `validate_EDI_Williams_fields.m`;
- `main_step2_validate_edi_normalization.m`.

This test constructs an independent polar annulus with duplicated crack-face
nodes, prescribes exact leading-order Williams mode-I, mode-II, and mixed-mode
displacements, and asks the unchanged EDI implementation to recover the
imposed SIFs. It reports both raw/input and (0.5\,\mathrm{raw}/\mathrm{input})
ratios across mesh refinements.

Decision rule:

- if raw/input tends to 1, the current EDI normalization is consistent;
- if raw/input tends to 2 while (0.5\,\mathrm{raw}/\mathrm{input}) tends
  to 1, the conversion requires an explicit factor (1/2);
- pure-mode cross leakage and the sign of recovered mode II are checked at
  the same time.

No EDI source-code correction has been made before this gate is run.


## Step 2 result: normalization error identified and corrected

The independent exact-field test was run locally on 2026-09-27.

For the coarse polar annulus mesh (Nr=8, Nth=64), the pre-correction
implementation returned approximately:

- pure mode I: KI_recovered/KI_input = 1.9722;
- pure mode II: KII_recovered/KII_input = 1.9722;
- mixed mode: both modal recovery ratios = 1.9722.

For the refined mesh (Nr=16, Nth=128), the ratios were:

- pure mode I: 2.0088;
- pure mode II: 2.0088;
- mixed mode: 2.0088 for both modes.

The non-imposed modal components converged toward zero, and pure mode II
was recovered with the same sign as the imposed field.

This establishes that the interaction integral itself was producing the
standard interaction quantity

\[
I^{(1,2)} = \frac{2}{E'}\left(
K_I^{(1)}K_I^{(2)} + K_{II}^{(1)}K_{II}^{(2)}
\right),
\]

whereas the conversion in `SIF_LEFM_interaction_EDI.m` had omitted the
factor (1/2).

The conversion has therefore been corrected to

\[
K = \frac{E'}{2K^{aux}} I.
\]

This correction is based on the independent analytical Williams-field
experiment, not on fitting the EDI result to the historical mirror/J value.

The Williams-field driver has now been converted into a regression test:
after the correction, recovered/input must approach one, pure-mode cross
leakage must approach zero, and mode-II sign must remain positive.


## Post-correction Williams-field regression

The corrected EDI implementation was rerun locally on 2026-09-27.

Fine mesh (Nr=16, Nth=128):

- pure mode I: KI_recovered/KI_input = 1.0044;
- pure mode II: KII_recovered/KII_input = 1.0044;
- mixed mode: both modal recovery ratios = 1.0044;
- pure-I cross leakage: KII = -2.2813e-7;
- pure-II cross leakage: KI = 6.8440e-7.

The coarse mesh (Nr=8, Nth=64) recovered approximately 0.9861 for both
modes, so the refinement trend brackets unity and the fine-mesh error is
about 0.44%.

The post-fix regression therefore passes the configured 2% tolerance.
Mode-II sign is confirmed independently by the pure-II exact field.

Status of the EDI normalization/sign gate: **PASSED**.


## Corrected same-field comparison

After the independently verified factor-1/2 correction, the same FEM field
was reevaluated locally:

| Method | KI | KII |
|---|---:|---:|
| old mirror/J | 0.33973 | 0.0033961 |
| corrected interaction EDI | 0.34004 | 0.0036663 |

At the matched outer radius 0.0008:

- mode-I difference relative to EDI is about 0.091%;
- mode-II difference relative to EDI is about 7.37%;
- the vector SIF difference relative to the EDI SIF norm is about 0.120%.

Because KII is only about 1% of KI in this control case, its componentwise
relative difference is much larger than the vector discrepancy. This result
does **not** yet identify mesh asymmetry as the cause. Radius/domain
sensitivity must be separated first.

## Step 3: same-field radius/domain sweep

The next gate solves the FEM field once and then performs:

1. a matched outer-radius sweep with r/lastLeg = 0.2, 0.3, 0.4, 0.5, 0.6
   for both the old contour and EDI outer domain;
2. an EDI inner-radius sweep at fixed r_outer/lastLeg = 0.5.

The old-method table also records JII/abs-integral cancellation, integrand
sign changes, a P/Q element-size proxy from detJ, and barycentric-margin
mismatch. This is intended to determine whether the residual small-mode
difference is primarily contour/domain sensitivity or is already correlated
with discrete P/Q asymmetry before deliberately asymmetric meshes are built.


## Step 3 result

The local same-field radius/domain sweep was completed on 2026-09-27.

Across matched outer radii r/lastLeg = 0.2 to 0.6, the old KI stayed near
0.3397 while the corrected EDI KI stayed within roughly 0.3384 to 0.3406.
The full-vector old-versus-EDI difference remained below 0.4%.

The small KII component was much more sensitive. The old-versus-EDI
componentwise KII difference varied non-monotonically from about 0.8% to
8.9%. The available P/Q detJ and barycentric mismatch indicators did not
track that difference monotonically.

At fixed EDI outer radius r_outer/lastLeg = 0.5, increasing
r_inner/r_outer from 0.05 to 0.30 changed KII from about 0.00377 to 0.00452,
while KI changed only mildly.

Therefore the residual KII discrepancy is not yet evidence of a mirror-mesh
error. Extraction-domain sensitivity of the small mode must first be
separated from the numerical FEM-field approximation.

## Step 3C: exact-field EDI domain sweep

The branch now contains
`verification/sif_audit/main_step3c_exact_field_edi_domain_sweep.m`.

It reuses the exact Williams-field machinery from Step 2 on a fixed fine
polar T6 mesh. Two EDI annulus sweeps are performed without solving a
physical boundary-value problem:

- outer-radius sweep: r_outer = 0.06, 0.08, 0.10, 0.12, 0.16 with
  r_inner/r_outer = 0.20;
- inner-radius sweep: r_inner/r_outer = 0.10, 0.20, 0.30, 0.40, 0.50 at
  fixed r_outer = 0.12.

Pure mode I, pure mode II, and mixed mode KI=1, KII=0.35 are evaluated.
If recovered/input ratios remain close to one with weak annulus dependence,
the much stronger KII sensitivity seen in Step 3 can be attributed mainly
to the finite-element approximation of the physical crack-tip field rather
than to the EDI formulation itself.


## Step 3C result

The exact-field EDI annulus sweep was completed locally on 2026-09-27.

For the outer-radius sweep at fixed r_inner/r_outer = 0.20, the recovered
pure-mode ratios ranged approximately from 0.978 to 1.001. The total range
was about 2.18% for mode I and 2.19% for mode II. The mixed-mode KI and KII
ratios followed the same scalar trend.

For the inner-radius sweep at fixed r_outer = 0.12, the recovered/input
ratios ranged approximately from 0.977 to 0.996. The total range was about
1.87% for both modes.

Cross-mode leakage remained negligible throughout, at roughly 1e-7 to 1e-5.
Pure mode I, pure mode II, and the mixed field all showed essentially the
same annulus dependence.

Interpretation:

- the corrected EDI formulation does not show a special mode-II instability
  on the exact Williams field;
- the remaining annulus dependence in this exact-field test is almost
  mode-independent and is consistent with interpolation/quadrature effects
  on the fixed T6 mesh;
- this approximately 2% exact-field annulus sensitivity is far smaller than
  the roughly 20%+ spread of the small KII component seen in the numerical
  FEM crack-tip field in Step 3B;
- therefore the strong KII-only sensitivity in Step 3B is primarily a
  property of the numerical FEM field and its interaction with the
  extraction domain, not an intrinsic mode-II defect of the EDI
  formulation.

One caveat remains: the exact Williams displacement field is sampled at T6
nodes and differentiated through the finite-element interpolation. A short
mesh-refinement check should confirm that the approximately 2% annulus
oscillation decreases with refinement before the deliberately asymmetric
mesh experiment.


## Adoption of canonical S0 mirror-reflected mesh

After reviewing the Step 3C topology, the synthetic Williams-field mesh was
changed to a strict mirror-reflected construction.

The canonical S0 mesh is now generated by:

1. constructing only the upper half-annulus, 0 <= theta <= pi;
2. triangulating its polar cells;
3. reflecting every upper T3 node through x2=0;
4. reflecting the complete upper T3 connectivity;
5. reversing reflected triangle orientation to keep all elements CCW;
6. sharing only the theta=0 radial line;
7. keeping theta=+pi and theta=-pi crack-face node IDs distinct;
8. converting the completed mirror-symmetric T3 mesh to T6.

Therefore S0 is mirror symmetric in node coordinates, T3 connectivity, and
the resulting T6 midside-node geometry. The previous Step 3C result obtained
with the same-diagonal synthetic topology is retained as historical audit
evidence but is superseded for the canonical reference calculation.

`validate_EDI_Williams_fields` now uses `MeshTopology='mirror_reflected'`
by default, while `legacy_same_diagonal` remains available only for
controlled comparison/reproduction.

Step 3C has been updated to use the canonical S0 mesh explicitly and prints
the maximum mirror-coordinate error plus confirmation that the two negative-x
crack faces have distinct node IDs.


## Canonical S0 rerun: Step 2 and Step 3C

The mirror-reflected S0 mesh was run locally on 2026-09-27.

Step 2 regression:

- maximum mirror-coordinate error: 0.000e+00;
- upper/lower crack faces remain distinct;
- fine mesh (Nr=16, Nth=128): KI recovery ratio = 1.0044 for pure mode I;
- fine mesh: KII recovery ratio = 1.0044 for pure mode II;
- pure-mode cross leakage is approximately 1e-14;
- regression status: PASS at 2% tolerance.

Compared with the previous same-diagonal synthetic topology, exact mirror
reflection reduces pure-mode cross leakage from roughly 1e-6--1e-5 to
machine precision, while leaving the imposed-mode recovery essentially
unchanged.

Step 3C on canonical S0:

- outer-radius sweep recovered/input range: about 2.18% for mode I and
  2.19% for mode II;
- inner-radius sweep recovered/input range: about 1.87% for both modes;
- pure mode I, pure mode II, and mixed mode share essentially the same
  annulus dependence;
- cross-mode leakage is at machine precision, roughly 1e-14.

Therefore the approximately 2% annulus oscillation observed in the exact
Williams-field test is not caused by the former non-mirrored connectivity.
It persists on a strictly mirror-reflected mesh and is therefore attributed
to finite T6 interpolation/quadrature of the nodally sampled singular field
on this fixed radial/angular discretization. The topology asymmetry mainly
manifested as tiny cross-mode leakage, which disappears on S0.

This strengthens the Step 3 conclusion: the much larger KII-only annulus
sensitivity of the numerical FEM crack-tip solution is a property of that
discrete physical field and its extraction, not an intrinsic mode-II defect
of the corrected EDI formulation.


## Step 3D: canonical S0 mesh-refinement study

A refinement driver has been added for the strictly mirror-reflected S0 mesh.
It uses three levels: (Nr,Ntheta) = (8,64), (16,128), and (32,256).

For every level the driver prints the complete mesh size: T3 vertices,
T3 elements, T6 nodes, T6 elements, displacement DOFs, mirror-coordinate
error, and whether the two crack faces remain distinct. It also opens a
full-mesh figure for each level by default.

The exact mixed Williams field KI=1, KII=0.35 is then reused for the same
outer- and inner-annulus sweeps as Step 3C. Pure-mode recovery and cross-mode
leakage are still checked at the baseline annulus for each level. The main
convergence quantity is the range of recovered/input ratios over each
annulus sweep. If the roughly 2% oscillation is a discretization effect,
these ranges should decrease with refinement while S0 cross-mode leakage
remains at machine precision.


## Step 3D result: S0 refinement

The canonical S0 refinement study was run locally on 2026-09-27.

Total meshes:

- L1 (Nr=8, Ntheta=64): 585 T3 vertices, 1024 T3/T6 elements,
  2193 T6 nodes, 4386 displacement DOFs;
- L2 (Nr=16, Ntheta=128): 2193 T3 vertices, 4096 elements,
  8481 T6 nodes, 16962 DOFs;
- L3 (Nr=32, Ntheta=256): 8481 T3 vertices, 16384 elements,
  33345 T6 nodes, 66690 DOFs.

All three meshes have zero reported mirror-coordinate error and distinct
upper/lower crack-face node IDs.

Baseline exact-field recovery at annulus [0.024,0.12] is not monotone with
refinement: mixed-mode vector error is about 0.42%, 1.80%, and 0.86% for
L1, L2, and L3 respectively. Cross-mode leakage remains at machine
precision.

Outer-annulus recovery ranges decrease overall from about 4.25% at L1 to
2.19% at L2 and 1.90% at L3. The inner-annulus range is non-monotone:
about 0.94%, 1.87%, and 0.78%.

Therefore the expected clean monotone convergence of annulus sensitivity is
not yet demonstrated. The oscillatory behavior is consistent with the
analytic radial weight having a sharp piecewise gradient whose circular
boundaries cut elements and are integrated by Gauss-point inclusion/
exclusion. Refinement changes the relative alignment of those boundaries
with the T6/Gauss layout, producing aliasing-like oscillations.

## Step 3E: FE-consistent weight-function check

Before any deliberate mesh-asymmetry experiment, the EDI implementation is
now given an optional `WeightFunction='fe_nodal'` mode. The default remains
`analytic_radial` so prior audit results are reproducible.

For `fe_nodal`, radial q values are assigned at all T6 nodes and q-gradient
is computed using the same T6 shape-function derivatives as the displacement
field. Gauss points are no longer accepted/rejected solely because their
radius falls on one side of an analytic annulus boundary; elements that
straddle a q transition contribute through the FE-interpolated gradient.

`main_step3e_compare_edi_weight_functions` compares both weight choices on
all three canonical S0 refinement levels using the same exact mixed Williams
field and the same outer/inner annulus sweeps. A substantial reduction and
smoother refinement trend with `fe_nodal` would identify annulus-boundary
quadrature aliasing as the main source of the residual exact-field
oscillation.


## Step 3E result: FE-consistent q removes annulus aliasing

Step 3E was run locally on 2026-09-27 using the canonical S0 meshes.

For the legacy clipped analytic radial weight, the maximum annulus-sweep
range remains at percent level: about 4.25%, 2.19%, and 1.90% for the
three refinement levels.

For the FE-consistent nodal weight, the maximum ranges collapse to roughly:

- L1 (8,64): 4.58e-4;
- L2 (16,128): 6.25e-5;
- L3 (32,256): 7.97e-7.

On the finest mesh the recovered mixed field is effectively KI=1 and
KII=0.35 across every tested inner and outer annulus, to the displayed
precision. This is a reduction of the annulus sensitivity by more than four
orders of magnitude relative to the clipped analytic-q implementation at
the same refinement level.

Interpretation: the percent-level oscillations seen in Steps 3C--3D were
caused predominantly by cut-element/Gauss-point aliasing at the analytic
circular q-transition boundaries. They are not an intrinsic limitation of
the interaction-integral formulation and are not a mode-II instability.
The FE-nodal q behaves as the numerically consistent EDI weight on the T6
mesh.

## Step 3F: physical-field check of the EDI weight

Before changing the production/default EDI weight, a same-field physical
FEM comparison has been added. The two weight functions are applied to the
same solved two-leg crack field, so no mesh or displacement field changes
between them. The driver reports outer- and inner-annulus ranges and the
old-mirror/J versus EDI differences at r/lastLeg=0.5. If FE-nodal q also
substantially reduces domain sensitivity there, it can be promoted to the
canonical EDI weight for the later mesh-asymmetry audit.


## Step 3F result: FE-nodal q on the physical FEM field

Step 3F was run locally on 2026-09-27. Both EDI weight functions were
applied to the identical solved two-leg crack field.

Compared with the clipped analytic radial weight, the FE-nodal weight
substantially reduces domain sensitivity:

- outer KI range: 2.1684e-3 -> 3.8928e-4 (factor about 5.6);
- outer KII range: 4.2970e-4 -> 4.5627e-5 (factor about 9.4);
- inner KI range: 8.9880e-4 -> 6.2595e-6 (factor about 144);
- inner KII range: 8.5631e-4 -> 7.3387e-6 (factor about 117).

At r/lastLeg=0.5, the old mirror/J versus EDI comparison changes from
(analytic q) KI=0.33973 vs 0.34004 and KII=0.0033961 vs 0.0036663 to
(FE-nodal q) KI=0.33973 vs 0.33977 and KII=0.0033961 vs 0.0034604.
The relative vector difference decreases from 1.2019e-3 to 2.2939e-4,
about a factor of 5.2. The componentwise KII discrepancy decreases from
7.37% to 1.86%.

On the FE-nodal EDI inner-radius sweep, KI stays near 0.33978 and KII near
0.003461 over inner/outer ratios 0.05--0.30. Thus the strong KII variation
previously seen in Step 3B was predominantly an artifact of the clipped
analytic q integration, not of the physical crack-tip field itself.

Conclusion: FE-nodal q is adopted as the canonical EDI weight for all
subsequent asymmetry-audit experiments. `analytic_radial` is retained as a
legacy/reproduction option until the audit is complete; the production
default is not changed yet.


## Step 4: graded concentric-ring mesh experiment

A new synthetic crack-tip mesh family has been added to match the preferred
concentric-ring construction used for illustration. The inner circle is
uniformly subdivided, successive radial rings are geometrically spaced, and
the angular nodes are staggered by half a sector on alternating rings. This
produces a graded near-equilateral triangular mesh whose element size grows
smoothly away from the crack tip.

For nominal angular increment dtheta, the equilateral outward-triangle
radial ratio is

q_eq = cos(dtheta/2) + sqrt(3) sin(dtheta/2).

The number of radial intervals is selected from q_eq, then the actual
constant ratio q is adjusted slightly so the last ring lands exactly on the
prescribed outer radius. With the default Ntheta=64, r0=0.005 and r1=0.20,
the construction uses 46 radial intervals; q is approximately 1.0835, very
close to the equilateral target. Because an exact tiling of a circular
annulus by equilateral triangles is geometrically impossible while all ring
nodes remain on concentric circles, triangle quality is reported explicitly
instead of claiming every element is exactly equilateral.

The baseline S0 mesh is created by building the upper half and reflecting
its coordinates and connectivity exactly. Three controlled lower-half
asymmetry cases are generated from the same ring radii: a smooth angular
shift at unchanged nominal density, a roughly 1.5-times coarser lower
angular spacing, and a 2-times coarser lower angular spacing. The positive
x radial line is shared, while the negative-x upper/lower crack faces keep
distinct node IDs.

`main_step4_graded_ring_asymmetry_experiment` prescribes exact pure-I,
pure-II, and small-mixed Williams fields on each mesh. The historical
mirror/J extractor is swept over five contour radii, while the canonical
interaction EDI uses FE-nodal q on a fixed annulus. Since the exact SIFs are
known, both extractors are assessed against truth; neither is used as a
reference for the other. Full meshes are plotted by default and mesh size,
quality, grading ratio, extraction errors, and old-contour sensitivity are
printed.


## Step 4 mesh correction after visual inspection

The first graded-ring implementation was rejected after inspection of the
exported mesh. The staggered rings had been formed by inserting theta=0 and
theta=pi endpoints into a half-sector-shifted angular sequence. Consequently
the first and last angular segments on every staggered ring were only half
the nominal size. This violated the intended uniform ring subdivision and
created localized poor elements along both x-axes, including the negative-x
crack seam. The earlier Step 4 numerical extraction results are therefore
superseded and must not be used for conclusions about asymmetry.

The corrected construction now divides every half-ring uniformly. Adjacent
rings use M and M+1 equal angular intervals alternately, which interlaces the
nodes while keeping theta=0 and theta=pi as exact nodes of every ring. Thus
there are no forced half-size seam segments and the crack faces remain a
straight, clean negative-x boundary.

The radial growth target is now explicitly described as a near-equilateral
shape target, not as an exact equilateral construction. For Ntheta=64 the
default corrected mesh uses Nr=45 and q approximately 1.08543. Independent
geometry reproduction predicts S0 triangle quality min about 0.823, fifth
percentile about 0.865, median about 0.957, with no triangles below quality
0.8 and crack-seam minimum quality about 0.865.

Hard S0 validation checks were added: within-ring segment-length spread must
be at roundoff level, crack-face nodes must lie on x2=0, and crack-seam
triangle quality must be at least 0.8. Step 4 now prints these quantities
before reporting any SIF results.

## Step 4C: local crack cut in the approved literal parent lattice

The preceding M/M+1 half-ring construction is retained as an experimental
audit record. The approved parent for this step is the literal closed
annular lattice in `build_literal_ring_lattice.m`, with the same N=64 equal
angular segments on every ring and alternating angular phase
0, Delta-theta/2, 0, Delta-theta/2, and so on. Its radial spacing remains
r0=0.005, r1=0.20, Nr=46, q approximately 1.08349620. The parent builder and
its triangulation are unchanged.

The separate `build_literal_ring_crack_cut_mesh.m` inserts the ray
x2=0, x1<0 after this closed parent T3 mesh has been constructed. Triangles
whose interiors intersect the ray are split locally. New intersection
vertices lie on their original edges, with one shared intersection per
edge before the faces are separated. A phase-shifted ring intersects the
ray on its polygon chord, at radius r*cos(pi/64); moving that intersection
to the analytic circle would alter the approved parent geometry and is
therefore avoided.

Original seam vertices are retained, with their roundoff-sized `sin(pi)`
y coordinates set to exact zero in the cut mesh. The returned parent mesh
retains the unsnapped coordinates. Duplicate crack-face IDs have exactly
identical coordinates and belong exclusively to their respective sides.
Unsplit triangles that touch the lower seam only substitute lower-face
IDs. Away from the ray, original vertex IDs, coordinates and triangle
connectivity must remain unchanged. Parent-element provenance is recorded
for every final T3 element, distinguishing actual subdivision from seam-ID
substitution.

The automatic T3 geometry audit must pass before T6 midside nodes are
created. The resulting T6 mesh is then audited separately, including the
topological separation of face midside nodes. The gate checks:

- unchanged nodes and connectivity outside the ray-intersection neighborhood;
- parent/child coverage and positive signed T3 areas;
- face nodes on exact x2=0 with distinct, coincident upper/lower IDs;
- no elements crossing the ray and no edges connecting the two crack faces;
- near-cut quality Q >= 0.70 and minimum angle >= 25 degrees, using
  Q = 4 sqrt(3) area / (a^2+b^2+c^2).

`main_step4c_preview_literal_crack_cut` is a geometry-only driver. It shows
the full mesh, a negative-x seam zoom, and a zoom where the crack meets the
inner boundary. Split children are colored by side, and coincident face
nodes are shown using different markers at the same coordinates. Optional
`Visible='off'` and `OutputDir` arguments support headless PNG, MAT and JSON
exports. The MAT output preserves both the original parent and final cut
meshes for direct comparison.

This step makes no changes to existing experimental mesh builders or SIF
drivers, and performs no SIF calculation. Passing the geometry gate is a
prerequisite for any later SIF experiment; it does not by itself establish
an SIF result or approve a change to the SIF experiment.

### Step 4C geometry result (MATLAB R2023a, 2026-09-27)

Both T3 and T6 geometry gates passed. Of the 5,888 original triangles,
46 intersect the ray through their interiors and are bisected. Another
92 touch the ray only at an existing vertex; the 5,750 outside this
neighborhood retain their exact connectivity. All 2,984 original nodes
off the ray retain their exact coordinates. The largest roundoff-only
normalization of an original on-ray y coordinate is 2.45e-17.

The cut introduces 23 shared edge intersections and 47 lower-face copies,
giving 3,078 T3 nodes and 5,934 triangles. The 9,011 topological edges then
produce 12,089 T6 nodes. Each crack face has 47 T3 vertices and 93 T6
nodes, with distinct IDs and exactly equal coordinates across the seam.
The mesh has 220 boundary edges and Euler characteristic one.

Minimum near-cut quality is 0.7512661928 and minimum angle is
30.08398398 degrees. The minimum signed area is 5.46207368e-8; the maximum
relative parent/child area discrepancy is 2.18e-16. The local bisections
explain the quality reduction from the parent minimum of approximately
0.993576. The gate floors 0.70 and 25 degrees allow a modest margin below
the expected cut geometry; the historic 0.8 seam gate would reject these
necessary local bisections and is not reused.

All three exported views were inspected: the full mesh, negative-x seam,
and inner boundary entry. Existing mesh builders, the parent builder,
`T3toT6_fast`, and the SIF experiment were not modified. The geometry
driver calls the approved parent builder and `T3toT6_fast`; it does not
run the experimental builders or any SIF calculation.

`test_literal_crack_cut` also passed: it locks the approved ring geometry
and expected counts, then verifies rejection of 15 corruptions, including
off-cut edits, nonpositive areas, overlapping children, crack bridges,
welded T3/T6 face nodes, sliver elements, invalid midpoints and incomplete
face-node lists.


## Step 4D: exact SIF validation on the approved literal crack-cut mesh

The approved parent geometry is now frozen: 64 equal chord segments on
every circular ring, alternating 0/half-sector phase, 46 radial intervals,
and the validated local negative-x crack cut from Step 4C.

Before prescribing the analytical field, the Williams displacement helper
was made crack-face aware. This is essential because the upper and lower
crack-face nodes have identical coordinates after the cut. Coordinates
alone would map both to theta=+pi in atan2 on many MATLAB builds. The helper
therefore accepts explicit upper/lower T6 face-node IDs and enforces
theta=+pi on the upper face and theta=-pi on the lower face.

`main_step4d_literal_mesh_sif_validation` performs the first accepted SIF
calculation on this mesh. It prescribes four exact Williams fields:
pure I (1,0), pure II (0,1), mixed 5% (1,0.05), and mixed 1% (1,0.01).
The last case reflects the small-mode-II regime of the physical two-leg
control problem.

The historical mirror/J extractor is swept over circular radii 0.02--0.14.
The canonical interaction EDI uses FE-nodal q and is checked on four
annuli, with [0.02,0.12] as the reference domain. Both methods are compared
directly with the imposed exact KI/KII; neither extractor is treated as
truth. The driver reports pure-mode cross leakage, vector error, contour/
domain ranges, and the old P/Q interpolation mismatch diagnostics.


## Step 4D result: exact SIF recovery on the approved literal crack-cut mesh

Step 4D was run locally on 2026-09-27. The geometry gate passed before any
SIF extraction. The cut mesh has 3078 T3 nodes, 5934 T3/T6 elements, 12089
T6 nodes, 47 T3 crack-face pairs and 93 T6 crack-face pairs. The near-cut
minimum quality is 0.751266 and the minimum angle is 30.084 degrees.

On the exactly symmetric S0 mesh, the historical mirror/J method shows
machine-level cross-mode leakage in the pure-mode tests. At the reference
contour rI=0.08 it recovers KI=1.0002135 for pure mode I and KII=0.9995729
for pure mode II. The old-method contour ranges over rI=0.02--0.14 are
6.865e-4 for KI and 6.641e-4 for KII, i.e. modest sub-0.1% variation.

The P/Q diagnostic mismatch at rI=0.08 is at roundoff level: median detJ
mismatch 4.47e-15 and median barycentric-minimum mismatch 4.05e-15. This
confirms that the approved S0 lattice is not only geometrically symmetric
but also supplies effectively identical mirrored interpolation stencils to
the old Ishikawa--Kitagawa--Okamura decomposition.

The canonical FE-nodal EDI is substantially more accurate on the same exact
fields. At the reference annulus [0.02,0.12], pure-I and pure-II errors are
about 8e-7, with cross leakage about 2e-13--4e-13. Its domain ranges are
about 2.42e-7 for KI and 7.48e-7 for KII.

For the 5% and 1% mixed-mode cases, the old method preserves the imposed
small KII very well: at rI=0.08 the KII errors are -2.14e-5 and -4.27e-6,
respectively. FE-nodal EDI errors are roughly 4.14e-8 and 8.29e-9.

Scientific conclusion: on a truly mirror-symmetric discrete mesh, the old
mirror/J mode-separation method is internally sound and highly accurate.
Therefore any degradation observed after controlled lower-half asymmetry can
be attributed to loss of discrete mirror correspondence rather than to an
intrinsic flaw of the continuum parity decomposition. The FE-nodal EDI now
serves as the mesh-general comparison method, while exact KI/KII remain the
ground truth in synthetic tests.


## Step 4E: controlled asymmetry-to-modal-contamination curve

Step 4D established the exact symmetric baseline. Step 4E now perturbs only
the lower-half geometry of that same approved crack-cut mesh while freezing
the upper half, both crack faces, all node IDs, and all T3/T6 connectivity.
For each lower-half T3 corner with polar angle -pi<theta<0,

theta_new = theta + alpha*dtheta*sin(theta).

The radius is preserved exactly. The perturbation vanishes on both x-axes
and reaches a maximum magnitude alpha*dtheta near the lower vertical axis.
Default alpha values are 0, 0.05, 0.10, 0.20, 0.40, and 0.80 sector widths.
This creates a one-parameter family in which topology and resolution count
are unchanged; only discrete mirror correspondence is progressively lost.

The driver `main_step4e_literal_mesh_asymmetry_curve` repeats the exact
pure-I, pure-II, 5%-mixed, and 1%-mixed Williams tests. The old mirror/J
method is swept over rI=0.02--0.14 and the FE-nodal EDI is checked on three
annuli. A compact detective table reports: normalized geometric mirror
mismatch; old P/Q detJ and barycentric mismatch; false KII generated from
pure mode I; false KI generated from pure mode II; relative KII error for
KII/KI=0.01; the same quantities from EDI; and minimum mesh quality.

The intended causal test is therefore direct:

discrete mirror mismatch -> P/Q interpolation mismatch -> modal contamination.

Because the exact SIFs remain prescribed on every mesh, neither extraction
method is treated as truth. The expected discriminator is that the old
mirror/J false cross-mode terms grow with asymmetry while FE-nodal EDI
remains nearly invariant.


## Step 4E result: smooth geometric asymmetry

Step 4E was run locally on 2026-09-28. The lower-half angular perturbation
increases the normalized median mirror mismatch from roundoff to 0.5656
while preserving radii and connectivity. Mesh quality remains acceptable:
the minimum quality decreases only from 0.7513 to 0.7208 and the minimum
angle from 30.08 to 28.15 degrees.

The old mirror/J interpolation diagnostics respond strongly to the imposed
asymmetry. At rI=0.08, median detJ P/Q mismatch rises from roundoff to about
0.0498, while the barycentric mismatch reaches about 0.0877. Pure-mode
cross leakage also appears: false KII generated from pure mode I grows from
roundoff to roughly 4.1e-4 by alpha=0.4, and false KI generated from pure
mode II reaches roughly 2.0e-4. The increase is clear over small-to-moderate
asymmetry, although it saturates and becomes mildly non-monotone at the
largest perturbation.

The FE-nodal EDI remains effectively invariant: false cross-mode terms stay
at roughly 1e-9--4e-8 and the 1%-mixed KII relative error remains of order
1e-6 or below. This confirms that the contamination is specific to the
mirror-dependent extraction rather than a generic consequence of degrading
the mesh geometry.

An important nuance is that the 1%-mixed KII error of the old method does
not grow monotonically with alpha; it remains below about 6.3e-4 in relative
magnitude. Thus the pure-mode leakage metric is a more sensitive detector of
broken parity than the scalar KII error of this particular mixed field. The
result suggests that the old method is more robust to smooth coordinate
perturbations than the initial hypothesis implied, despite measurable modal
contamination.

## Step 4F: connectivity-only asymmetry (A_conn)

To isolate the effect most directly relevant to arbitrary unstructured FEM
meshes, Step 4F keeps every T3 vertex coordinate exactly identical to S0 and
changes only lower-half triangulation. Natural annular-cell diagonals are
flipped in deterministic nested fractions 0, 0.10, 0.25, 0.50, and 1.00.
The upper half and crack neighborhood are frozen. Therefore geometric mirror
symmetry is exact, but mirrored interpolation stencils are progressively
different.

`main_step4f_literal_mesh_connectivity_asymmetry` repeats the exact pure-I,
pure-II, 5%-mixed and 1%-mixed tests and reports old P/Q mismatch, false
cross-mode SIFs, mixed-mode KII error, EDI response, contour/domain ranges,
and mesh quality. This is the cleanest test of whether connectivity alone
can break the historical mirror-mode separation.


## Step 4F implementation correction after first run

The first Step 4F run aborted at the 25% flip level because three old/J
contour points were reported outside the mesh. The console diagnostics
revealed the cause before any scientific interpretation was attempted. At
10% there were 276 requested diagonal flips but only 543 distinct changed
triangles, whereas 552=2*276 are required if every flipped diagonal owns a
disjoint pair of triangles. At 25%, 690 flips changed only 1313 distinct
triangles instead of 1380. Thus some selected diagonals shared an owner
triangle. Sequentially applying such overlapping flips overwrote an earlier
child connectivity and created local holes/nonconforming topology.

The A_conn generator has been corrected by first constructing a deterministic
maximal matching of candidate lower-half annular diagonals: each T3 triangle
may participate in at most one flip. The requested fractions are now taken
from this disjoint pool, so the family remains nested while every selected
flip changes exactly two unique triangles.

Additional hard gates now require positive element areas, edge incidence at
most two, unchanged boundary-edge set, conserved total area, unchanged Euler
characteristic, frozen T3 coordinates, and exactly two changed triangles per
flipped diagonal. The aborted first Step 4F run is superseded and provides
no accepted SIF result.


## Step 4G result: direct mirrored-stencil diagnostic

Step 4G was run locally on 2026-09-28 using the same connectivity-only A_conn family. The new reflected-parent-triangle diagnostic identifies the topology loss much more directly than the earlier scalar detJ/barycentric medians.

For flip fractions 0, 0.10, 0.25, 0.50, and 1.00, the fraction of contour P/Q samples having exactly mirrored T3 parent triangles is respectively 1.0000, 0.8750, 0.78333, 0.5625, and 0.21667. The corresponding absolute false KII produced from an exact pure-mode-I field is approximately 0, 1.2874e-4, 1.4871e-4, 4.0903e-4, and 4.6878e-4.

Thus the useful topology metric is not the median detJ or median barycentric mismatch. Those medians remain near roundoff through the 50% case because more than half of the sampled P/Q pairs still occupy mirrored-equivalent local patches. The direct mismatch distribution reveals the altered pairs immediately: its 95th percentile jumps to about 0.50418 already at 10% flips, while the fraction of exactly mirrored T3 pairs decreases monotonically with connectivity asymmetry.

Across these five deterministic cases, the absolute false KII and the broken-pair fraction (1 - fraction_exact_mirror_T3) have a descriptive Pearson correlation of about 0.95. This is not treated as a universal law, but it strongly supports the mechanism: loss of mirrored parent-element correspondence is associated with increasing parity contamination in the historical mirror/J separation.

Step 4G therefore closes the synthetic topology mechanism test: exact node-coordinate mirror symmetry alone is insufficient; the old decomposition also relies on sufficiently mirrored discrete interpolation stencils. FE-nodal interaction EDI does not require this correspondence.


## Step 5: bridge to the physical two-leg FEM mesh

After Step 4G established the synthetic topology mechanism, a physical-mesh bridge was added. `main_step5_physical_mesh_stencil_audit` solves the existing two-leg Crack-Path-style FEM control once and evaluates the historical mirror/J and canonical FE-nodal EDI extractions over several contour radii.

For each radius the driver now reports the direct reflected-T3 parent-stencil diagnostics introduced in Step 4G, including median/p95/max mismatch and the fraction of exactly mirrored P/Q parent triangles. It also reports old-vs-EDI KI/KII differences and the combined vector difference on the same displacement field.

This step is intentionally descriptive rather than causal because exact physical-field SIFs are not known. The causal evidence remains the exact synthetic Steps 4D--4G; Step 5 asks whether the same broken-stencil signature is present in the actual Crack-Path-style mesh where the earlier small-KII discrepancy was observed.


## Step 5 result: physical two-leg mesh has no exact mirrored contour stencils

Step 5 was run locally on 2026-09-28 on the existing two-leg physical FEM control. Across r/lastLeg = 0.2--0.6, the fraction of exactly mirrored P/Q parent T3 pairs is zero at every tested radius. The direct mirrored-T3 mismatch median is about 0.31--0.52, the 95th percentile about 0.75--0.79, and the maximum about 0.89--1.11. Thus the actual Crack-Path-style mesh strongly violates the discrete mirrored-stencil condition identified in Steps 4F--4G.

On the same physical displacement field, old/J and FE-nodal EDI remain very close in KI, but KII differences are visibly contour-dependent. Relative to EDI KII, the old-minus-EDI KII discrepancy is approximately -3.41%, +1.93%, -0.92%, -1.40%, and +3.02% for r/lastLeg = 0.2, 0.3, 0.4, 0.5, and 0.6. For r/lastLeg >= 0.3, EDI KII is nearly constant (about 0.0034603--0.0034628; range about 0.072%), whereas old/J KII varies by about 4.4% over the same radii.

The reported correlation between broken-pair fraction and KII discrepancy is NaN for a simple reason: the broken-pair fraction equals one at every tested radius, so it has zero variance. This does not weaken the bridge; it means the physical mesh is already fully outside the exact mirrored-stencil regime throughout the tested contour family.

Because exact physical-field SIFs are unknown, Step 5 alone cannot assign the remaining KII difference to one extractor. It only shows that the physical mesh exhibits exactly the discrete condition under which the synthetic tests demonstrated old-method sensitivity.

## Step 6: exact Williams field sampled on the physical two-leg mesh

To remove the remaining ambiguity, Step 6 samples exact local Williams fields directly on the actual two-leg Crack-Path mesh. The extraction radii remain at or below 0.6 of the final crack-leg length, so the audit domain lies wholly inside the straight final leg and does not reach the earlier kink. Exact pure-I, pure-II, and 1%-mixed fields are prescribed in the local final-leg frame, with explicit upper/lower crack-face branch assignment.

This makes KI/KII known exactly while retaining the real physical-mesh node layout and connectivity. The historical mirror/J and FE-nodal EDI methods can therefore be compared against truth on the very mesh that produced the Step-5 discrepancy. This is the final bridge needed before deciding how the production crack-growth workflow should handle SIF extraction on asymmetric meshes.


## Step 6 result: exact Williams field on the physical two-leg mesh

Step 6 completed locally on 2026-09-30. The final physical crack leg has length 0.0016 and carries 41 T3 face corners and 81 T6 face nodes per face, with exactly one shared upper/lower T6 node at the mathematical tip. The stored face orientation agrees with the +e2/-e2 local-frame convention.

The exact-field result closes the remaining attribution gap. On the actual Crack-Path mesh, the historical mirror/J method produces appreciable false cross-mode content even though the prescribed field is exactly pure. For pure mode I, false KII ranges from about 5.0e-4 to 1.38e-3 over r/Llast=0.2--0.6. For pure mode II, false KI ranges from about 5.2e-5 to 1.93e-3. The direct reflected-T3 exact-pair fraction is zero at every tested radius.

For the exact 1%-mixed field, the old relative KII error is approximately -5.44%, +1.01%, -1.27%, -0.072%, and -1.60% for r/Llast=0.2, 0.3, 0.4, 0.5, and 0.6. FE-nodal EDI is more accurate at every radius: approximately +1.22%, -0.463%, -0.0342%, +0.00165%, and -0.00457%, respectively. At r/Llast=0.5, the old KII error is about 0.072% while EDI is about 0.00165%; at r/Llast=0.6 the old error is about 1.60% while EDI is about 0.0046%.

The EDI degradation at the smallest radius is consistent with under-resolution of the inner q transition: with ncoh=40, htip=Llast/ncoh=4e-5, while the Step-6 inner radius is 0.1*r. Thus r/Llast=0.2 gives r_inner=3.2e-5 (< htip), whereas r/Llast=0.5 gives r_inner=8e-5 (=2 htip). The best-resolved domains are therefore the mid/outer cases rather than the smallest contour.

Scientific conclusion: the full chain is now demonstrated. Exact discrete mirror correspondence gives excellent old/J recovery (Step 4D); controlled loss of geometric or connectivity symmetry creates modal contamination (Steps 4E--4G); the actual physical Crack-Path mesh has no exact mirrored parent-element pairs (Step 5); and exact Williams fields sampled on that physical mesh reproduce the old-method sensitivity directly (Step 6). FE-nodal EDI does not require mirrored stencils and is substantially more accurate on the same asymmetric physical mesh.

Before promoting EDI into the production crack-growth path, the remaining implementation-specific gate is auxiliary-field derivative sensitivity, because the present prototype still obtains auxiliary displacement gradients by centered finite differences.


## Step 7: auxiliary-derivative sensitivity gate

Step 6 establishes the physical-mesh extraction difference using exact truth. The remaining implementation-specific question for EDI is whether its accuracy depends materially on the centered finite-difference step used for auxiliary displacement gradients.

`SIF_LEFM_interaction_EDI` now accepts `AuxDerivativeScale`, a positive multiplier on the existing step h=max(1e-7,1e-5*r). The default remains 1.0 for backward compatibility. The extractor also reports median and maximum mismatch between auxiliary strain computed from finite-difference displacement gradients and strain recovered from the analytical auxiliary stress through the constitutive matrix.

`main_step7_edi_aux_derivative_sensitivity` uses the actual two-leg physical mesh and exact pure-I, pure-II, and 1%-mixed Williams fields. It fixes the EDI annulus at the well-resolved Step-6 choice [0.05,0.5] Llast, for which r_inner=2 htip at ncoh=40, and sweeps derivative-scale multipliers 0.01--30. A broad plateau in recovered KI/KII together with small strain-stress mismatch would show that the remaining EDI error is governed by FE interpolation/domain discretization rather than the numerical auxiliary derivative step.


## Step 7 result: EDI auxiliary derivative step is not controlling the SIFs

Step 7 completed locally on 2026-09-30. On the actual two-leg physical mesh, with exact Williams fields and the well-resolved FE-nodal annulus [0.05,0.5] Llast, the recovered SIFs form an extremely broad plateau versus auxiliary finite-difference scale.

For the exact 1%-mixed field, KI_EDI remains about 0.99996 and KII_EDI about 0.01000016 for derivative-scale multipliers from 0.01 through 10. The relative KII error is approximately 1.645e-5 (0.00165%) over this entire range. Even at scale 30 it changes only to about 1.714e-5. Pure-I and pure-II recovery are similarly invariant.

At the default derivative scale 1, the auxiliary strain-vs-analytical-stress mismatch is very small: mode-I median about 6.3e-8 and maximum about 1.64e-4; mode-II median about 1.4e-8 and maximum about 9.63e-5. The pointwise maximum grows with very large derivative scales and becomes pathological at scale 30, while the medians remain small and the integrated SIFs remain nearly unchanged. This is consistent with a small number of centered finite-difference probes crossing the crack displacement branch cut near a face; scale 30 is therefore a stress test, not a recommended operating point.

In the Step-7 domain, the floor in h=max(1e-7,1e-5*r) dominates, so the scale sweep effectively varies h from about 1e-9 to 3e-6, a factor of 3000. The resulting SIF plateau demonstrates that the remaining EDI error is governed by FE interpolation/domain discretization rather than the auxiliary finite-difference step.

Operational conclusion: `AuxDerivativeScale=1` is validated for the present workflow. Closed-form auxiliary derivatives remain a possible cleanup/refinement, but they are not required before using FE-nodal EDI as the canonical SIF extractor on asymmetric meshes.


## Step 8: production hole-crack preflight

The synthetic/physical-control audit is now sufficiently mature to enter the actual hole-crack production workflow, but the first production step is deliberately an exact-field preflight rather than an immediate crack-angle sweep.

`compute_SIF_for_stage2_compare.m` processes one existing Stage-II displacement field with both extractors without repeating the FEM solve. The historical branch uses the circular J / Ishikawa--Kitagawa--Okamura separation; the comparison branch uses FE-nodal interaction EDI with `AuxDerivativeScale=1`. The EDI outer radius defaults to the same `G2.tip.radiusJ` used by the historical contour. The inner radius is chosen as `max(0.1*r_outer, 2*h_tip)`, where `h_tip` is estimated from T3 edge lengths of elements incident on the crack-tip node. The adapter also returns the direct mirrored-parent-element diagnostics.

`main_step8_production_crack_path_preflight.m` builds the actual Stage-II hole+short-crack mesh at theta=0 by default, recovers the production crack-face node sets, assigns exact local Williams pure-I, pure-II and 1%-mixed displacement fields, and compares both SIF extractors directly with known truth. It also plots the full production mesh and a crack-tip zoom suitable for later thesis export.

This preflight is the gate before the full production theta sweep. It checks that the FE-nodal EDI domain is adequately resolved on the actual hole-crack mesh and that the conclusions from Steps 4--7 transfer to this specific production discretization. Only after this gate is accepted should KII(theta) and the zero-KII direction be compared between the two extractors.


## Step 8 result: exact-field preflight on the actual production hole-crack mesh

Step 8 completed locally on 2026-09-30 for the production Stage-II mesh at theta=0 deg. Stage I located the initiation point at phi=359 deg on the discretized hole boundary. The appended crack length is a0=Llast=0.004. The collapsed crack has 7 T3 face nodes and 13 T6 face nodes per face, with exactly one shared node at the mathematical tip.

The production mesh again has no exactly mirrored parent-element pairs on the historical contour: fraction_exact_mirror_T3=0, with reflected-T3 mismatch median 0.358, p95 0.726, and maximum 0.903.

The local T3 tip-edge scale is coarse relative to the earlier audit mesh: h_tip,min=4.09e-4, h_tip,median=4.70e-4, h_tip,max=5.39e-4. With r_outer=0.002=0.5 Llast, the automatic EDI rule therefore selected r_inner=9.40e-4=2 h_tip,median, leaving a transition width about 1.06e-3.

On exact pure Mode I, old/J returned KI=1.00273 and false KII=-3.353e-3, whereas FE-nodal EDI returned KI=1.00000096 and false KII=1.628e-4. On exact pure Mode II, old/J returned false KI=2.477e-3 and KII=0.998962, whereas EDI returned false KI=-1.899e-4 and KII=1.000026.

For the exact 1%-mixed field (KI=1,KII=0.01), old/J returned KII=0.010874, a relative KII error of +8.741%, while FE-nodal EDI returned KII=0.010163, a relative error of +1.631%. The full-vector error is about 2.87e-3 for old/J and 1.63e-4 for EDI. Thus the production mesh reproduces the same qualitative conclusion as the preceding audit: the mirror-based extractor is much more strongly contaminated by the asymmetric interpolation layout, while EDI is substantially more accurate.

The EDI result is not yet as close to exact as on the refined two-leg audit mesh. This is consistent with the much coarser production tip discretization and the relatively narrow q-transition annulus forced by the 0.5 Llast outer radius. Therefore the production theta sweep should not rely on a single EDI domain. Because post-processing is cheap compared with remeshing and solving, the next driver should evaluate several safe EDI outer radii on every identical FEM solution and compare the resulting zero-KII direction. This makes root stability, rather than one local KII percentage, the production acceptance criterion.


## Step 9 result: production theta sweep and a newly exposed KII-sign issue

Step 9 completed locally on 2026-09-30 for theta=-12:1:8 deg and matched extraction radii r/a0=0.50,0.65,0.80. The FE-nodal EDI KII(theta) curves are extremely stable across all three domains. Around the zero crossing, EDI gives KII(theta=0) about -1.287e-3 and KII(theta=1 deg) about +2.10e-3 for all three domains, leading to interpolated signed roots 0.37990, 0.37876, and 0.37871 deg. The total root spread across the three EDI domains is only 0.00119 deg.

The historical decomposed-J output behaves fundamentally differently. For theta<0 it remains predominantly positive even where EDI gives a stable negative KII, and near theta=0 its reported sign changes with contour radius. At r/a0=0.50 it has two sign changes; at 0.65 and 0.80 it has none. The previously printed old/J root range of 0 deg was therefore misleading because only one radius produced any finite sign-crossing root.

Inspection of the legacy conversion exposed a conceptual reason: after modal separation, JII is a Mode-II energy-release contribution and is quadratic in the physical KII amplitude. The code converts it using sign(JII)*sqrt(abs(JII)*Eeff). But in exact LEFM JII=KII^2/Eeff, so sign(JII) cannot recover the physical sign of KII. A negative JII can only arise from numerical integration error, not from a negative physical KII. The earlier exact-field audit used positive Mode-II amplitudes and therefore did not reveal this limitation.

This does not invalidate the historical method as a magnitude extractor: on symmetric meshes it accurately recovers |KII| and its modal energy. It does mean that its returned sign must not be used as a signed local-symmetry root indicator. For the historical decomposed-J branch, the meaningful production diagnostic is the minimum of |KII| (or JII), whereas the interaction EDI provides a genuinely signed KII because the interaction term is linear in the target mode amplitude.

Step 9 code has therefore been hardened without changing legacy returned values. The old/J sign crossings are now labeled legacy diagnostics only; the driver reports the theta minimizing |KII| for each historical contour. The EDI root stability remains the signed local-symmetry result of interest. `SIF_LEFM_circle2_debug` now also records explicitly that physical KII sign is not recoverable from JII alone.


## Step 10 result: explicit positive/negative Mode-II sign audit

Step 10 completed locally on 2026-09-30 and decisively separates sign loss from mesh-asymmetry contamination.

On the exactly symmetric S0 mesh, exact pure Mode-II fields with KII=+1 and KII=-1 produce exactly the same historical modal result: KII_legacy=+0.99957 and JII=2.2731e-4 in both cases. Likewise, exact mixed fields (KI,KII)=(1,+0.01) and (1,-0.01) both produce KII_legacy=+0.0099957 and the same JII=2.2731e-8. Thus the historical decomposed-J procedure recovers the Mode-II magnitude very accurately on S0 but contains no information about the physical sign of KII. FE-nodal interaction EDI recovers +1/-1 and +0.01/-0.01 with errors below about 1e-6 and 1e-8, respectively.

On the actual production mesh, exact pure KII=+1 and -1 again give identical old/J magnitude 0.99896 and identical JII=4.3243e-6, while EDI recovers +1 and -1 with signed error about 2.61e-5. For the mixed +/-1% cases, the asymmetric mesh additionally contaminates the historical magnitude: old/J gives +0.010874 for true +0.01 but +0.0076715 for true -0.01. EDI gives +0.010163 and -0.0098374. The EDI pair is nearly perfectly linear: their mean is about +1.63e-4, equal to the pure-Mode-I false-KII leakage on this mesh, while their half-difference is about 0.0100002, essentially the prescribed Mode-II amplitude.

Step 10 therefore establishes two independent effects. First, modal JII is quadratic and cannot encode sign even on a perfect symmetric mesh. Second, asymmetric interpolation stencils contaminate the modal magnitude. Interaction EDI is signed because the interaction term is linear in the target modal amplitude and it also remains substantially less sensitive to the asymmetric interpolation layout.


## Step 11: symmetry-enforced production parity benchmark

The centered-hole geometry and remote vertical tension are analytically symmetric, so the rightmost initiation point should be phi=0 deg and a normal crack at theta=0 should be a pure Mode-I benchmark. The current Stage-I numerical detector instead selects phi=359 deg, which mixes the SIF-extractor question with a separate initiation-point discretization bias.

`main_step11_production_symmetry_parity.m` therefore preserves the production Stage-II mesh/solver but overrides only the initiation geometry to the exact rightmost symmetry point: x_star=center+[R,0], n_mat=[1,0], t_hat=[0,1]. It runs paired +/-theta production meshes and evaluates both methods at r/a0=0.50,0.65,0.80. The EDI acceptance checks are KII(0) approximately zero, KII(-theta) approximately -KII(+theta), KI(-theta) approximately KI(+theta), one signed zero crossing near theta=0, and root stability across domains. The historical branch is evaluated only through |KII| evenness because Step 10 proved that decomposed JII cannot supply the physical KII sign.

This benchmark should be completed before refining the 0.379-deg EDI root from Step 9. If the symmetry-enforced benchmark returns a stable zero root, the Step-9 offset can be attributed primarily to the Stage-I 359-deg initiation bias plus ordinary mesh discretization rather than to the EDI extractor.


## Step 11 result: symmetry-enforced production benchmark passes

Step 11 completed locally on 2026-09-30. The numerically detected Stage-I maximum occurs at phi=359 deg, whereas the exact symmetry point is phi=0. The sampled effective tangential stress at phi=0 is 3.5588003469, while the numerical maximum is 3.5986892523, a relative bias of 1.1084%. This confirms that the one-degree initiation offset is not merely a tie-breaking artifact; the current Stage-I recovered-stress/interpolation procedure itself is slightly asymmetric.

After imposing the exact rightmost initiation geometry (phi=0, n_mat=[1,0], t_hat=[0,1]) while keeping the production Stage-II meshing and solver unchanged, FE-nodal EDI satisfies the expected symmetry very well. At theta=0, KII is only about -3.10e-5, -2.80e-5, and -2.92e-5 for r/a0=0.50,0.65,0.80. The corresponding interpolated zero-KII roots are 0.009026, 0.008180, and 0.008510 deg, giving a domain spread of only about 8.46e-4 deg.

The EDI odd-parity defect in KII is at most about 0.43--0.57% of the KII scale over theta=+/-3 deg, while the KI even-parity defect is only about 4.7e-5--1.0e-4 relative. Thus the production EDI extractor and Stage-II workflow preserve the expected centered-hole symmetry to high accuracy once the initiation geometry is fixed.

The historical magnitude-only branch is much less symmetric: the maximum relative mismatch in |KII(-theta)| versus |KII(+theta)| is about 4.0%, 12.4%, and 5.0% for r/a0=0.80,0.65,0.50, respectively. This is consistent with the earlier conclusion that asymmetric mirrored interpolation contaminates the decomposed-J modal magnitude even when sign is no longer considered.

Interpretation: the Step-9 EDI root near +0.379 deg was not a physical symmetry-breaking crack direction. It arose mainly because Stage I supplied a crack-start normal at phi=359 deg rather than the exact symmetric phi=0 point. With the initiation geometry corrected, the production Stage-II EDI root collapses to about +0.0085 deg, effectively zero at the present mesh resolution. The next remaining production issue is therefore Stage-I initiation-point accuracy/symmetry, not SIF extraction.


## Step 12: Stage-I recovered-stress / scattered-interpolation refinement audit

Source inspection exposed a more basic issue in the current Stage-I boundary sampler. For a circular hole, `n_out=[cos(phi),sin(phi)]` is explicitly used elsewhere as the material-side normal, pointing from the hole center into the solid. However, `sample_hole_boundary_stress.m` currently defines `xq = xb - eps_shift*n_out` while its comment says that the query is shifted into the body. The minus sign actually moves the query radially into the hole cavity. At the current Npoly=240 settings, hhole=2*pi*R/240≈7.854e-4 and eps_shift=0.25*hhole≈1.963e-4, so the query circle has radius R-eps≈0.029804 for R=0.03: it is well inside the polygonal cavity rather than in the FEM material domain.

The current use of `scatteredInterpolant` masks this domain error because it constructs an interpolation over the cloud of recovered nodal stresses rather than respecting the nonconvex FEM topology; in particular, it can return finite values at points lying inside the hole. This provides a concrete mechanism for the Stage-I 359-deg bias and means that refinement of the legacy procedure alone may converge to an interpolation artifact rather than to the material-side boundary stress.

Step 12 therefore keeps `StressExt` unchanged and separates three postprocessing variants on the same Stage-I solve: (A) the exact current behavior, cavity-side query plus scattered interpolation; (B) corrected material-side query plus the same scattered interpolation; and (C) corrected material-side query plus topology-respecting T6 interpolation of the same recovered nodal stresses. The default coupled refinement family uses Npoly=[120,180,240,360,480], keeps the boundary sampling fixed at 1440 angles, and scales h_arc=2*pi*R/Npoly, Hmin=Hhole=h_arc, Hmax=20*h_arc. For every case the driver reports whether the query points actually lie in the FEM domain, the right-hand peak angle, the stress at phi=0, peak bias, local +/-phi symmetry defect, left/right peak mismatch, and mesh counts.

The key acceptance question is no longer simply whether the selected peak approaches phi=0 with refinement. First, the material-side query must lie in the actual FEM domain. Then material-side scattered and topology-respecting T6 interpolation should converge toward one another and toward the symmetric peak. Only after that evidence should the production Stage-I sampler be changed.


## Step 12 result: Stage-I refinement isolates the dominant postprocessing defect

Step 12 completed locally on 2026-09-30 for Npoly=[120,180,240,360,480] with fixed angular sampling Nphi=1440. The legacy cavity-side query has fraction_query_in_actual_FEM_domain=0.000 at every refinement level. Therefore every legacy stress query point lies outside the actual material mesh, inside the hole cavity. `scatteredInterpolant` nevertheless returns finite values because it interpolates over the recovered-stress point cloud without respecting the hole topology.

The legacy sequence is correspondingly erratic and non-convergent as an initiation locator: the right-hand peak angle moves -8.5, -0.25, -1.0, +1.25, +0.5 deg as Npoly increases, while the local symmetry defect varies from 3.58% to 0.15% to 1.30% to 0.052% to 0.993%. At the current Npoly=240 level it reproduces the earlier phi=-1 deg bias, with a 1.108% peak-over-phi0 stress excess and 1.304% local symmetry defect.

Both corrected material-side variants place 100% of query points in the actual FEM domain at all refinement levels. They also agree closely with one another. At Npoly=240, material-side scattered interpolation gives a right-window peak at -0.5 deg with only 8.29e-5 relative peak-over-phi0 excess and 2.41e-4 local symmetry defect; topology-respecting T6 interpolation gives +0.5 deg with 7.40e-5 peak excess and 2.68e-4 symmetry defect. At Npoly=480, both place the discrete right-window maximum at +0.25 deg; the peak excess is only 3.30e-5 for scattered and 1.55e-5 for T6, while the local symmetry defects are 5.94e-5 and 7.15e-5.

The residual +/-0.25--0.5 deg movement of the discrete maximum should not be interpreted as a physical angular bias. Nphi=1440 gives a 0.25-deg sampling increment, and the stress maximum becomes extremely flat: by Npoly=480 the difference between the sampled peak and phi=0 is only O(1e-5) relative. Thus the discrete argmax becomes dominated by tiny remeshing/postprocessing noise even though the stress field itself is nearly symmetric.

Scientific conclusion: the dominant Stage-I defect is the sign of the radial query shift, not insufficient mesh refinement. The current legacy sampler evaluates stresses in the cavity and should not be used for production initiation. `scatteredInterpolant` is a secondary concern: once the query is moved to the material side, scattered and topology-respecting T6 interpolation are already close and converge toward the same nearly symmetric stress profile. For production, topology-respecting T6 interpolation is preferable because it cannot silently bridge across the hole topology. A final audit should compare this recovered-nodal-stress T6 interpolation with direct stress evaluation from the FEM displacement gradient at the same material-side query points before changing the production Stage-I implementation.


## Step 13 result: recovered-T6 versus direct-from-U Stage-I stresses

Step 13 completed locally on 2026-09-30 for the same Npoly=[120,180,240,360,480] coupled-refinement family and material-side offset eps=0.25 h_hole. The two estimators converge strongly toward the same stress field. The relative difference in sigma_tt(phi=0) decreases from 2.92e-3 at Npoly=120 to 2.18e-4 at Npoly=480. Over the +/-12-deg right-hole window, the relative Linf difference decreases from 5.05e-3 to 2.46e-4 and the relative L2 difference from 2.49e-3 to 1.26e-4. Thus there is no evidence of a distinct limiting stress introduced by StressExt recovery.

The recovered-T6 field is substantially smoother. The direct-U angular roughness indicator is about 6.9, 4.4, 6.0, 6.2, and 4.1 times larger than the recovered-T6 value over the five refinement levels. Nevertheless, direct-U roughness also decreases rapidly with refinement, and both estimators identify the same +0.25-deg discrete peak on the finest Npoly=480 mesh. The remaining +/-0.25--0.5-deg peak hopping is therefore best interpreted as grid/remeshing noise on an extremely flat maximum, not as a systematic estimator bias.

The traction-free diagnostics provide an independent consistency check. Recovered-T6 and direct-U give almost identical material-side residuals at every mesh: max|sigma_nn|/sigma_tt-scale falls from about 1.38e-2 to 3.28e-3, and max|sigma_nt|/sigma_tt-scale from about 7.0e-3 to 1.76e-3. Their near-linear reduction with h is consistent with evaluating at an offset eps proportional to h rather than exactly on the boundary.

At Npoly=480, recovered-T6 gives sigma_tt(0)=3.589569 and direct-U gives 3.588785, a relative difference of only 2.18e-4; the full-window Linf difference is 2.46e-4. Recovered-T6 has much smaller peak-over-phi0 bias (1.55e-5 versus 2.69e-4) and lower angular roughness (1.54e-4 versus 6.24e-4). This supports retaining recovered nodal stresses for the production initiation criterion, but only with material-side topology-respecting T6 interpolation. Direct-U should be retained as an audit/reference estimator rather than the default production field.

Before changing production, one remaining sensitivity check is useful: vary the material-side offset fraction eps/h_hole (for example 0.05,0.10,0.25,0.50) on one or two refined meshes. This will verify that the chosen offset does not materially bias sigma_tt(0) or the inferred initiation position and will quantify the approach of sigma_nn and sigma_nt to zero as the query approaches the traction-free boundary.


## Step 14: material-side query-offset sensitivity

After Step 13 showed convergence of recovered-T6 and direct-from-U stresses, Step 14 varies the material-side sampling distance eps/h_hole over [0.05,0.10,0.25,0.50] on two representative coupled-refinement meshes, Npoly=240 and 480. Each mesh is solved only once; both estimators are then evaluated at identical material-side query rings.

The experiment reports sigma_tt(0), the discrete right-window peak angle/value and peak-over-phi0 bias, +/-phi symmetry error, relative sigma_nn and sigma_nt residuals, recovered/direct differences, and within-mesh ranges over all offsets. The primary questions are whether sigma_tt(0) and the inferred initiation point are insensitive to the chosen offset, and whether the traction residuals decrease as eps approaches the exact-circle boundary. Because the FE hole boundary is polygonal, a nonzero residual floor at very small eps is possible until the geometry is further refined; this should be interpreted together with the Npoly dependence rather than as a failure of the stress estimator.


## Step 14 result: offset sensitivity and boundary-limit interpretation

Step 14 completed locally on 2026-09-30 for Npoly=240 and 480 with eps/h_hole=[0.05,0.10,0.25,0.50]. All material-side query points remain inside the FEM domain. The raw sigma_tt(0) value is not offset-insensitive: for recovered-T6 it ranges by about 2.69% over the four offsets at Npoly=240 and by 1.37% at Npoly=480. This is not a defect; it reflects the physical radial stress gradient because the query point moves farther into the body as eps increases.

The traction-free residuals behave exactly as desired. On Npoly=480, recovered-T6 max|sigma_nn|/sigma_tt-scale decreases from 6.48e-3 at eps/h=0.50 to 8.29e-4 at 0.05, while max|sigma_nt| decreases from 3.50e-3 to 3.90e-4. Npoly=240 shows the same trend at approximately twice the residual level. Linear extrapolation of the residuals versus eps/h gives small nonzero intercepts that drop by roughly a factor of four to five when Npoly doubles (for recovered-T6, sigma_nn intercept about 7.9e-4 -> 1.8e-4 and sigma_nt about 1.7e-4 -> 3.4e-5), consistent with a geometry/discretization floor rather than an estimator failure.

Recovered-T6 and direct-U continue to agree. At the finest Npoly=480 mesh, their relative sigma_tt(0) difference is 4.19e-4 at eps/h=0.05 and only 5.38e-5 at 0.50; the full-window Linf difference remains below 4.41e-4 for all tested offsets. Direct-U remains rougher near the boundary, so recovered-T6 is still the preferred production field.

Most importantly, sigma_tt(0) is nearly linear in eps/h over the tested range. A linear fit gives an estimated zero-offset boundary value of approximately 3.6154--3.6159 for recovered-T6 at Npoly=240 and 3.6172--3.6174 at Npoly=480 (depending on whether all four or only the smaller offsets are used). The corresponding direct-U zero-offset estimates are about 3.60877 and 3.61563. Thus boundary extrapolation is more principled than selecting an arbitrary fixed nonzero offset. The recovered-T6 boundary-limit estimate changes by only about 0.04% between Npoly=240 and 480 and agrees with the Npoly=480 direct-U boundary estimate within about 0.05%.

Production implication: do not freeze eps=0.25 h as the physical boundary stress. Instead, evaluate recovered-T6 stresses at several small material-side offsets and extrapolate sigma_tt(phi,eps) to eps=0. This preserves topology, reduces dependence on an arbitrary sampling distance, and provides a natural traction-free consistency check. A local angular fit can then be applied to the extrapolated boundary stress field to estimate the initiation angle without discrete-sample hopping.


## Stage-I production redesign started

Following Steps 12--14, the audited findings have now been translated into a guarded production redesign while retaining the historical path for reproducibility. The new root-level sampler `sample_hole_boundary_stress_v2.m` uses recovered nodal stresses from `StressExt`, locates several material-side query rings in the actual FEM topology, interpolates through the containing T6 element, and linearly extrapolates each stress component to eps->0. The default production offsets are eps/h_hole=[0.05,0.10,0.25].

`find_hole_initiation_point_v2.m` then finds the discrete tensile maximum of the extrapolated boundary field and refines its angular position with a periodic local quadratic fit (default five angular samples). The exact circular point and local frame are reconstructed from the fitted angle rather than copied from the nearest sampling point. `run_stage1_hole_initiation.m` centralizes the solve, postprocessing, and initiation-point selection. The historical cavity-side scattered-interpolation path remains available through `C.stage1.method='legacy_scattered'`.

`cfg_hole_initiation.m` now selects `boundary_extrapolated_t6` by default, with Nphi=1440, shift fractions [0.05,0.10,0.25], linear radial extrapolation, and five-point angular fitting. The principal production drivers have been wired through the centralized Stage-I workflow. No Stage-II SIF-selection logic was changed in this redesign step.

Before relying on the redesigned Stage-I output downstream, `verification/sif_audit/main_step15_stage1_redesign_preflight.m` should be run. It exercises the actual production wrapper, evaluates the legacy sampler on the same FEM solution for reference, reports discrete and fitted initiation angles, boundary-limit stress/load, radial-fit diagnostics, and centered-hole symmetry relative to the exact candidate peaks at 0 and 180 deg.


## Step 15 result: redesigned radial Stage-I chain passes; five-point angular locator is still too local

Step 15 completed locally on 2026-09-30 on the current Npoly=240 production mesh. The redesigned production chain gives a boundary-limit stress maximum of 3.6166977, consistent with the Step-14 zero-offset estimate, and lambda_ini=82.9486. All three material-side query rings have in-domain fraction 1.0, and the maximum relative radial-fit RMSE in sigma_tt is only 1.84e-5. Thus the material-side topology-respecting T6 sampling plus eps->0 extrapolation passes its production preflight.

The extrapolated boundary field is also highly symmetric: the local +/-phi stress defect near the right peak is 1.62e-4 and the left/right peak mismatch is only 3.13e-6. sigma_tt(0)=3.615716 and sigma_tt(pi)=3.615215. These diagnostics are far better than the historical cavity/scattered result, which still selects phi=-1 deg and gives the lower sampled maximum 3.598689 on the same FEM solution.

However, the current five-point angular quadratic refinement does NOT pass the strict centered-hole direction gate. The discrete maximum is at -0.5 deg and the fitted vertex is -0.4880 deg, even though the underlying extrapolated field is nearly symmetric. The fitted peak exceeds sigma_tt(0) by only about 9.82e-4 in absolute stress, or roughly 2.7e-4 relative, so the local five-point fit is locking onto a tiny mesh-scale ripple on a very flat physical maximum. Its small RMSE (1.28e-6) only shows that a quadratic fits those five neighboring samples well; it does not prove that this very narrow neighborhood identifies the continuum peak.

Conclusion: retain the redesigned radial boundary-stress estimator, but do not yet freeze the five-point angular peak locator. The angular regression window should be tied to the boundary mesh angular scale h_hole/R and made wide enough to average mesh-scale oscillations while shrinking consistently under refinement. A dedicated window-width sensitivity audit is the next gate.


## Step 16 result: mesh-scaled angular regression works; benchmark has two equivalent maxima

Step 16 revealed that the raw cross-mesh summary must respect the centered-hole symmetry. At Npoly=240 the global discrete maximum lies near the right-hand exact peak (phi=0), whereas at Npoly=480 tiny numerical differences make the left-hand exact peak (phi=pi) marginally larger. These are physically equivalent initiation sites for the centered benchmark; comparing +0 deg and +/-180 deg as if they were different directions creates a false ~180-deg 'spread'.

Measured relative to the nearest exact symmetry-equivalent peak (0 or pi), the fitted angular offsets for half-width factor c=[1,1.5,2,3] are approximately: Npoly=240: [+0.0381,-0.0174,+0.0329,+0.00604] deg; Npoly=480: [-0.00270,-0.03845,+0.00216,+0.00430] deg. The c=3 window is clearly the most mesh-stable in this benchmark: maximum absolute offset about 0.0061 deg and cross-mesh signed spread about 0.00174 deg. Its quadratic-fit residual remains small (relative RMSE about 1.30e-4 on Npoly=240 and 3.69e-5 on Npoly=480), and the fitted curvature remains consistent with the narrower windows.

The c=0.5 case is too narrow on Npoly=480 because it contains only three angular samples, while the audit deliberately requires at least five. Factors c=1--3 all remove the mesh-scale discrete peak hopping; c=3 gives the best centered-hole invariance. Because the angular half-width is c*h_hole/R, even c=3 shrinks to zero under mesh refinement (4.5 deg half-width at Npoly=240 and 2.25 deg at Npoly=480).

Two production-design consequences follow. First, the angular locator should be expressed by a mesh-scaled half-width rather than a fixed five-point neighborhood. Second, the initiation finder must distinguish 'where is the local peak?' from 'which of several physically equivalent global peaks should a single-crack simulation select?'. Mesh noise should not choose between degenerate sites. The corrected Step-16 driver now reports local offsets relative to the nearest exact candidate peak rather than naively comparing raw angles.

Candidate production setting from the centered-hole audit: angular_fit_halfwidth_factor = 3.0. Before freezing it, verify on an intentionally asymmetric geometry that this window preserves a genuinely nonzero/unique peak rather than oversmoothing it.


## Step 17: centered-hole right half-domain verification

The centered-hole benchmark is now being reformulated as a right-half symmetry model rather than using the full domain with two physically equivalent initiation sites. The retained domain is x in [A/2,A], y in [-B,B]. The circular hole becomes a semicircular traction-free cutout on the vertical symmetry boundary. This removes left/right initiation degeneracy while retaining both y>0 and y<0 material, so later Stage-II upward/downward kinking is not imposed by the symmetry reduction.

`cfg_centered_half_domain.m` defines the benchmark. `geom_centered_half_hole.m` builds the concave half-domain polygon with Npoly/2 straight segments on the right semicircle. `solve_hole_only.m` now supports `C.bc.anchor_mode='symmetry_half_x'`: ux=0 is imposed on all T6 nodes on x=x_sym, with one uy gauge constraint to remove rigid vertical translation. The standard remote y-tension loading and stress recovery are otherwise unchanged.

`sample_hole_boundary_stress_v2.m` now also supports a finite nonperiodic circular arc through `C.stage1.phi_range`; full-circle behavior is unchanged when that field is absent. `find_hole_initiation_point_v2.m` supports an optional mesh-scaled angular half-width. The half-domain benchmark uses phi in [-pi/2,pi/2], 721 samples (0.25-deg spacing), and the Step-16 candidate half-width factor c=3.0. This setting is verification-specific and has not yet been frozen into the general asymmetric full-domain production configuration.

`verification/sif_audit/main_step17_centered_half_domain_verification.m` compares the half-domain solution with the right semicircle of an independently meshed full-domain solution at Npoly=240 and 480. It reports sigma_tt(0), full-curve Linf/L2 differences, symmetry defects, fitted right-peak angles, initiation loads, traction-free residuals, and the symmetry-BC displacement residual. Overlay, difference, convergence, and half-domain geometry plots are produced.


## Step 17 result: centered right-half benchmark validated

Step 17 completed locally on 2026-10-01 for Npoly=240 and 480. The symmetry boundary is enforced exactly at the algebraic level (max |ux| on x=x_sym equals zero). The half-domain and independently meshed full-domain right-semicircle boundary-limit stress fields converge rapidly toward one another: relative Linf difference decreases from 1.35e-3 at Npoly=240 to 1.50e-4 at Npoly=480, while relative L2 decreases from 5.36e-4 to 4.24e-5.

The fitted half-domain initiation angle converges to the exact rightmost symmetry point: +0.0170 deg at Npoly=240 and -0.00428 deg at Npoly=480. The corresponding full-domain right-peak fits are +0.00604 deg and +0.00551 deg. On the fine mesh the peak stresses are 3.617550 (half) and 3.617529 (full), and the initiation load factors are 82.92905 and 82.92953, respectively. Thus the symmetry-reduced benchmark reproduces the full-domain right-hand solution while removing the left/right degeneracy.

The half-domain traction-free residuals also decrease strongly under refinement: max|sigma_nn|/sigma_tt-scale from 9.56e-4 to 2.25e-4 and max|sigma_nt|/sigma_tt-scale from 2.83e-4 to 6.96e-5.

One diagnostic-only bug was found in the first Step-17 print/table output: the full-domain sigma_tt(0) value was indexed with the full-circle sorted-grid index after the full solution had already been interpolated onto the half-domain grid, so the reported values -1.421 and -1.424 were actually values from the wrong half-domain array index. This did not affect the curve comparison, fitted peaks, peak stresses, or initiation loads. The driver has been corrected to use the half-domain phi=0 index for both aligned curves.

Conclusion: the centered right-half model is accepted as the preferred symmetry benchmark for Stage I. It eliminates the artificial choice between two equivalent left/right initiation sites while retaining the complete upper/lower material domain needed for a later unbiased Stage-II KII=0 / theta=0 verification.


## Step 18 prepared: centered half-domain Stage-II signed-EDI symmetry test

The Stage-II continuation of the centered right-half benchmark is now implemented as a dedicated verification problem. The crack mouth is placed at the exact rightmost circular-hole point, independent of the small residual Stage-I fitted-angle error, so Stage II can be tested in isolation. The half-domain keeps the vertical symmetry boundary ux=0, which means that when reflected back to the full plate the Stage-II benchmark corresponds to a centered hole with two opposite, mirror-related short cracks. This is intentional: it removes left/right degeneracy while leaving the entire upper/lower material domain free, so the right-hand crack may be swept through positive and negative kink angles without prescribing theta=0.

`build_domain_centered_half_pencil.m` constructs the concave half-domain outer boundary with the semicircular notch and a sharp temporary two-face appendix. `build_stage2_centered_half_cracked_mesh_for_theta.m` meshes this geometry, identifies the two appendix faces, and collapses them to the zero-thickness crack midline. `solve_cracked_LEFM.m` now supports the same `symmetry_half_x` boundary condition used in Stage I.

`verification/sif_audit/main_step18_centered_half_stage2_symmetry.m` performs the signed interaction-EDI audit at Npoly=[240,480], theta=[-2,-1,-0.5,0,0.5,1,2] deg, and r_outer/a0=[0.50,0.65,0.80]. For every theta a fresh cracked mesh is solved once; all EDI radii are evaluated on that same displacement field. The driver reports KII(0)/KI, the odd-parity defect of KII(theta), the even-parity defect of KI(theta), a linear signed-root estimate, a sign-change/interpolated root, EDI-domain spread, tip mesh scale, and the exact ux=0 symmetry residual. The theta=0 collapsed mesh is plotted for the first refinement level.

The acceptance gate is: KII(0)/KI -> 0 under refinement, KII(-theta) approximately -KII(+theta), KI(-theta) approximately KI(+theta), and the signed-EDI local-symmetry root theta* -> 0 with small dependence on the EDI outer radius. Absolute KI values should not be compared directly with the earlier full-domain single-right-crack benchmark because the symmetry-reduced Stage-II model represents the mirrored two-crack configuration.


## Step 18 result: centered half-domain Stage-II symmetry benchmark passes

Step 18 completed locally on 2026-10-01 for Npoly=240 and 480, theta=[-2,-1,-0.5,0,0.5,1,2] deg, and r_outer/a0=[0.50,0.65,0.80]. The signed FE-nodal interaction EDI shows the expected parity: negative trial angles give negative KII and positive trial angles give positive KII, while KI remains nearly even in theta.

At theta=0, KII/KI is already between about -1.59e-5 and -6.78e-6 on Npoly=240 and between about -3.09e-6 and -1.28e-6 on Npoly=480, depending on EDI outer radius. Thus the worst zero-angle modal contamination decreases by a factor of about 5.1 under the 240->480 refinement. The KI even-parity defect improves by about an order of magnitude, from roughly 2.3--2.6e-4 to about 2.0e-5, while the KII odd-parity defect remains at only O(1e-5) and also decreases.

The signed local-symmetry root is essentially zero. Linear fits give |theta*| <= 9.1e-5 deg on Npoly=240 and <=1.17e-4 deg on Npoly=480; the sign-change/interpolated roots are at most about 1.66e-3 deg and 3.23e-4 deg, respectively. The apparent lack of monotonic improvement in the tiny linear-fit root is below the numerical noise scale and is not physically meaningful; the bracket root and parity defects give the clearer convergence evidence.

EDI-domain sensitivity is negligible: at theta=0 the relative KI spread over r_outer/a0=[0.50,0.65,0.80] is 1.58e-5 for Npoly=240 and 9.06e-6 for Npoly=480. The maximum |KII/KI| at theta=0 falls from 1.59e-5 to 3.09e-6. The symmetry boundary remains exact algebraically (max |ux|=0).

The geometry-only sharp-pencil edge identifier fails for this concave outer-boundary representation, but the existing temporary-mesh fallback consistently identifies the correct upper/lower appendix edges (64/65 for Npoly=240 and 124/125 for Npoly=480). This is an efficiency/cleanliness issue rather than a correctness issue; the collapsed meshes and SIF symmetry results are stable.

Conclusion: the centered right-half Stage-II benchmark is accepted. Together with Step 17, it closes the symmetric verification chain: Stage-I boundary-stress initiation converges to phi*=0, signed EDI gives KII(0)->0, KII is odd and KI even in theta, and the local-symmetry direction converges to theta*=0. The next scientific benchmark should therefore be a genuinely asymmetric full-domain hole configuration, where a nonzero unique initiation angle and nonzero Stage-II direction correction can be tested without symmetry degeneracy.


## Step 19 prepared: asymmetric full-domain Stage-I benchmark

After closing the centered symmetry chain in Steps 17--18, the next benchmark is a genuinely asymmetric full-domain circular hole. `cfg_asymmetric_full_domain.m` uses the existing A=0.30 m, B=0.10 m, R=0.03 m plate/hole dimensions but moves the hole center to [0.17,-0.02] m. Both x- and y-reflection symmetries are therefore broken. The minimum ligament is still 0.05 m, so this is not a near-contact geometry.

`verification/sif_audit/main_step19_asymmetric_stage1_validation.m` tests Npoly=[240,480] and angular half-width factors c=[1,1.5,2,3] on the same redesigned boundary-limit stress field for each mesh. It reports the discrete and fitted global initiation angles, fitted stress and initiation load, fit residual, radial extrapolation residual, traction residuals, and cross-mesh angle/stress changes. It also ranks independent local tensile maxima on the full circular boundary, reporting the relative stress gap and angular separation between the dominant and secondary peaks.

The purpose is to validate the Step-16 candidate c=3 away from symmetry. Acceptance requires a clearly preferred local maximum, a genuinely nonzero fitted initiation angle, stability under Npoly=240->480 refinement, and agreement of c=3 with narrower mesh-scaled windows. Only after this gate passes should the asymmetric Stage-II signed-EDI direction sweep be launched.


## Step 19 result: asymmetric Stage-I angle is stable; peak-ranking diagnostic needs clustering

Step 19 completed locally on 2026-10-01 for the full-domain hole centered at [0.17,-0.02] m. The redesigned radial estimator remains well conditioned: all three query rings stay inside the FEM domain; the maximum relative radial-fit RMSE decreases from 1.78e-5 at Npoly=240 to 4.76e-6 at Npoly=480. The traction residuals also fall strongly (sigma_nn relative maximum about 1.05e-3 -> 2.28e-4; sigma_nt about 7.33e-4 -> 1.58e-4).

The asymmetric initiation direction is clearly nonzero. For the mesh-scaled c=3 angular window, phi*= -1.57371 deg at Npoly=240 and -1.56061 deg at Npoly=480, a cross-mesh change of only 0.0131 deg. The fitted peak stress changes by only about 3.50e-4 relative and the initiation load by the same amount with opposite sign. c=2 gives essentially the same cross-mesh stability; c=1 is also stable in this particular case but was less robust in the centered symmetry benchmark. c=1.5 is visibly more sensitive to which discrete mesh ripple is selected. Taken together with Step 16, c=3 remains the preferred production candidate.

The initial Step-19 'dominant versus second local peak' output must not be interpreted as two competing physical initiation sites. Their fitted angles differ by only about 0.021 deg at Npoly=240 and 0.010 deg at Npoly=480, far below the mesh angular scale (1.5 and 0.75 deg) and far below the regression half-width. They are duplicate local-ripple detections that collapse onto the same broad physical maximum. Therefore the raw top-second stress gap of O(1e-5)--O(1e-6) is not a degeneracy measure.

`main_step19b_asymmetric_peak_hierarchy.m` was added as a zero-cost postprocessor for an existing O19 result. It performs angular non-maximum suppression after quadratic refinement, clustering candidate maxima whose separation is less than c_cluster*h_hole/R (default c_cluster=3). This produces a ranking of physically separated maxima without repeating any FEM solve. The independent-peak hierarchy should be checked before launching the asymmetric Stage-II direction sweep.


## Step 19B result: asymmetric Stage-I has a stable preferred physical peak

Step 19B clustered mesh-scale duplicate detections into physically separated maxima without repeating any FEM solve. For Npoly=240, the dominant independent peak is near phi=-1.553 deg with sigma=3.71282, while the second independent peak is near phi=-177.636 deg with sigma=3.66350. For Npoly=480, the corresponding peaks are near -1.561 deg with sigma=3.71407 and -177.617 deg with sigma=3.66444.

The dominant-to-secondary stress gap is therefore about 1.33% on both meshes (1.3286% and 1.3363%), and the two physical peaks are separated by about 176.1 deg. The ratio sigma2/sigma1 is approximately 0.9867 at both refinements. Thus the primary site is not overwhelmingly stronger, but its preference is highly mesh-stable and comfortably larger than the numerical postprocessing errors measured in Steps 14--19. Under proportional loading and before crack-induced redistribution, the secondary site would reach the same tensile threshold at roughly 1.35% higher load.

Conclusion: the asymmetric Stage-I benchmark is accepted for a single-crack first-initiation example. The production initiation point continues to use the global c=3 Stage-I fit from `find_hole_initiation_point_v2`; Step 19B is only a physical-peak hierarchy diagnostic. Because the secondary peak is close in strength, later multi-crack studies should not assume it remains inactive after the first crack grows; redistribution must be solved explicitly.


## Step 20 prepared: asymmetric full-domain Stage-II signed-EDI probe

With the independent Stage-I peak hierarchy accepted, the next step is deliberately modest: a signed-EDI direction probe on the Npoly=240 asymmetric full-domain benchmark before any refined root solve or Npoly=480 Stage-II sweep. `main_step20_asymmetric_stage2_probe.m` uses the redesigned Stage-I initiation point and local material frame, then remeshes and solves fresh Stage-II cracked geometries at theta=[-4,-2,-1,0,1,2,4] deg relative to the fitted Stage-I material normal.

For each theta, the FE-nodal interaction EDI is evaluated at r_outer/a0=[0.50,0.65,0.80], with r_inner=max(0.1*r_outer,2*h_tip). The driver reports KI, signed KII, KII/KI, tip mesh scale, interpolated sign-change roots, global linear-fit roots, and EDI-domain root spread. No extra root-refinement mesh is generated in this step.

The acceptance gate is that all three EDI domains exhibit the same KII(theta) sign trend and produce a consistent sign-change root inside the probe interval. If that condition is met, the next step will refine only the root neighborhood and then repeat the confirmed direction on Npoly=480.


## Step 20 result: signed response consistent, but the asymmetric direction correction is unresolved

Step 20 completed locally on 2026-10-02 for the asymmetric full-domain hole center [0.17,-0.02] m and Npoly=240. Stage I gives phi*=-1.57371281 deg, x*=[0.19998868461,-0.020823890498] m, boundary stress 3.71277044 at unit load, and lambda_ini=80.80219463. Seven freshly meshed Stage-II angles theta=[-4,-2,-1,0,1,2,4] deg were evaluated by signed FE-nodal EDI at r_outer/a0=[0.50,0.65,0.80].

The angular dependence is extremely clean: at theta=-4 deg KII is about -0.01394 and at theta=+4 deg about +0.01392 (KI about 0.365); near theta=0 the local slope dKII/dtheta is about 0.203 per rad, and the three EDI domains agree closely. However KII(theta=0)=[+5.20e-6,+5.82e-6,+1.93e-6] for the three domains, corresponding to KII/KI=[+1.42e-5,+1.59e-5,+5.29e-6]. These are comparable to the centered Step-18 numerical symmetry floor. They do not establish a physically nonzero angle correction.

The sign-change/interpolated roots are [-0.0014628,-0.0016402,-0.0005450] deg, mean -0.0012160 deg and domain spread 0.0010952 deg. The global seven-point linear fits instead yield [+0.00610,+0.00671,+0.00685] deg, illustrating that numerical ripple/nonlinearity dominates a correction this small. Neither should be reported as a measured physical kink angle; the supported interpretation is theta* approximately zero within current numerical resolution even though the Stage-I initiation point phi* is displaced by about -1.57 deg.

This finding is physically plausible: the location of maximum boundary hoop stress is naturally close to a locally pure opening orientation for a short radial appendix. It must not be claimed as an exact theoretical identity for finite a0 without further study. Changing the geometry merely to force a nonzero correction would skip an important model verification question.

## Step 21 prepared: controlled off-peak mouth perturbation

`verification/sif_audit/main_step21_asymmetric_initiation_perturbation.m` reuses an existing O20 Stage-I solution and deliberately moves the crack mouth in hole polar angle by delta_phi=[-3,-1.5,0,+1.5,+3] deg around the Stage-I fitted maximum (one and two nominal boundary-mesh angular scales for Npoly=240). Each location gets a fresh, normal-oriented full-domain cracked FEM solve, followed by signed EDI evaluation on the same three domain sizes as Step 20. The original Stage-I solution is not recomputed. A theta-correction proxy -KII(0)/(dKII/dtheta) is obtained from the independent Step-20 +/-1 deg sweep and explicitly labeled diagnostic only; a true root at the shifted mouth would require a fresh theta sweep.

Step 21 checks two hypotheses: (1) the near-zero baseline KII is reproducible under an independent rebuilt solve, and (2) deliberately off-peak positions produce a clear signed KII response above the Step-18/20 noise floor. If both hold, the nearly normal Stage-II direction at the true Stage-I peak is a supported physical/numerical result, rather than evidence that the EDI extractor is insensitive to the actual initiation geometry. Only then should a0- or geometry-sensitivity experiments be considered.


## Step 21 result: the near-normal direction is geometry-sensitive, not an insensitive signed-EDI selector

Step 21 completed locally on 2026-10-02 on the Npoly=240 asymmetric configuration at fixed a0=0.004 m. The hole mouth was moved by [-3,-1.5,0,+1.5,+3] degrees relative to the fitted Stage-I maximum phi*=-1.57371281 degrees, and a new normal short crack was meshed and solved at every point. All three FE-nodal EDI domains [0.50,0.65,0.80]*a0 gave consistent signed responses. At the reference domain r_outer/a0=0.65, the measured KII values were [-3.95836e-3,-2.01048e-3,+5.82449e-6,+1.91112e-3,+3.94075e-3], respectively; KI remained about 0.364--0.366. Therefore KII is strongly responsive to mouth relocation, reverses sign across the Stage-I peak, and rises by about three orders of magnitude above the near-zero baseline residual for the +/-3-degree shifts.

The independent Step-20 slope dKII/dtheta is about 0.200 per rad (approximately 3.49e-3 per degree). Combining that slope with the off-peak KII measurements gives linearized correction proxies of approximately [+1.135,+0.576,-0.00167,-0.548,-1.130] degrees for the five mouth offsets. These are NOT solved local-symmetry roots at shifted mouths; they merely predict a compensating direction. The numerical results imply a sensitivity of order dKII/d(phi_mouth) ~ 1.3e-3 per degree and d(theta_correction)/d(phi_mouth) ~ -0.38 in this finite-length/geometry configuration.

The zero-mouth-offset case exactly reproduces the Step-20 baseline KII for every EDI domain. This is deterministic same-configuration reproducibility, not independent discretization validation. The boundary_stress_over_peak=0.99987 printed at zero offset arises because Step 21 linearly interpolates the discrete boundary-stress samples whereas Stage I reports the fitted angular quadratic peak; do not interpret that small difference as a physical shift of the crack mouth.

IMPORTANT precision qualification: the Step-19 c=3 initiation angle shifts by 0.0131 deg between Npoly=240 and 480. The Step-21 observed mouth sensitivity means such a small angle perturbation could change KII by about 1.7e-5, corresponding to a linearized direction change of about 0.005 deg. This is larger than the apparent 0.0005--0.0016 deg Step-20 local roots, so those roots must not be reported as physical corrections. The robust conclusion is that Stage II is almost normal at the Stage-I stress maximum to within the presently resolved angular accuracy; this is not an assertion of an exact theta*=0 identity for arbitrary finite crack lengths.

## Step 22 prepared: refined-mesh crack-length sensitivity

`verification/sif_audit/main_step22_asymmetric_length_sensitivity.m` computes Stage I independently at Npoly=480 using the c=3 mesh-scaled angular window, then fixes that initiation point and its material-side normal for four separate Stage-II short-crack lengths a0=[0.002,0.004,0.006,0.008] m. It solves each cracked geometry at theta=0, measures signed KI/KII using three FE-nodal EDI annuli r_outer/a0=[0.50,0.65,0.80], and reports tip mesh scale and inner/outer-radius ratios, within-length domain spread, and KII/KI versus a0/R. The existing Step-20 Npoly=240,a0=0.004 solution is compared with the Npoly=480,a0=0.004 result, with an explicit caveat that this comparison changes both mesh density and the Stage-I fitted mouth angle. Optional Step-21 input prints the coarse-mesh mouth-sensitivity scale for that angle change.

The next gate is whether the near-zero normal-crack KII persists on a finer mesh and over the selected crack lengths, or becomes a resolved finite-length effect. This is a controlled finite-length sensitivity study, not an asymptotic proof and not a direction-root computation. If KII becomes substantially resolved, a separate theta sweep/root search will be needed for those selected lengths.


### Step 22 geometry-control qualification

The existing `build_stage2_cracked_mesh_for_theta.m` hardcodes `nArc=160` for the retained Stage-II hole arc. Thus Step 20 and the planned Step 22 length scan use the SAME fixed 160-point Stage-II hole polygon even as Npoly drives finer Stage-I geometry and smaller overall Stage-II FEM hmin. This improves FE resolution while keeping the Stage-II hole geometry fixed. Such a comparison is controlled for finite-length effects, but it does NOT constitute full geometric convergence of the Stage-II hole boundary. A subsequent isolated nArc sweep at fixed fine FEM resolution is needed before interpreting residual KII smaller than the potential polygon-geometry floor as a physical kink.


## Step 22 result: finite-length residual requires a geometry-convergence audit

The fine-mesh asymmetric Stage-I solution (Npoly=480) places the peak at phi*=-1.56061278 deg, compared with phi*=-1.57371281 deg on Npoly=240. The change is +0.01310003 deg. From Step 21, the coarse mouth-sensitivity estimate is dKII/dphi approximately 1.3072e-3 per degree at a0=0.004 m, implying a KII change of approximately 1.7124e-5 at that length from the angle shift alone. Therefore Step-20 versus Step-22 comparison at 4 mm is NOT a fixed-mouth convergence test.

With the fine-mesh Stage-I point held constant and theta=0 for all lengths, the reference EDI annulus r_outer/a0=0.65 gives:

| a0 (mm) | a0/R | KI | KII/KI |
|---:|---:|---:|---:|
| 2 | 0.0667 | 0.28853 | +6.7963e-5 |
| 4 | 0.1333 | 0.36617 | +1.6721e-5 |
| 6 | 0.2000 | 0.40943 | +3.5461e-5 |
| 8 | 0.2667 | 0.43766 | +9.1118e-5 |

The 2-mm value is less secure because r_inner/r_outer ranges from about 0.31 to 0.49, and the relative KII/KI domain spread is 3.90e-5. The 4-mm residual is tiny, and the near-identical 240/480 ratios do not establish convergence because the mouth moved between meshes. The 6-mm and especially 8-mm residuals are larger, with 8-mm KII/KI in [8.44e-5,9.11e-5] and domain spread 6.73e-6, but domain independence alone does not exclude a systematic boundary-polygon or mouth-geometry error.

Two geometry details are deliberately NOT conflated with physical finite-length effects: (i) the full-domain Stage-II builder currently retains only nArc=160 points along the remaining hole boundary, independent of Npoly; (ii) the temporary appended-hole mouth half-shift stays fixed at 0.1 mm, so its ratio to a0 varies substantially over this length study. Either can matter when interpreting very small signed KII values. In particular, the 8-mm residual corresponds to a putative angular correction comparable with the existing initiation-angle uncertainty, but the dKII/dtheta slope must be recomputed at 8 mm before estimating any direction.

## Step 23 prepared: same-mouth polygon and mesh convergence

The Stage-II builder now accepts an optional 'NArc' name-value argument (default 160, so historical callers retain their geometry). The focused verification driver 'verification/sif_audit/main_step23_asymmetric_fixed_mouth_convergence.m' takes an existing O22, reuses its fine Stage-I mouth EXACTLY and reuses its already computed 160-point fine-mesh Stage-II solutions at a0=4 and 8 mm. It adds independent cracked FE solves with NArc=320 and 480 on the same nominal fine FE h and, for NArc=480, a separate coarse nominal FE h at the identical mouth. It evaluates the same three signed FE-nodal EDI domains on each new solution.

The arc-resolution sweep necessarily remeshes the whole polygonal domain, so it measures combined geometric/remeshing variation. The second, fixed-NArc coarse/fine comparison isolates the nominal FE-size effect to the extent that enforced polygon edge discretization permits; the reported tip mesh scales should be used to check whether the local meshes truly differ. No new Stage-I fit is done in Step 23, preventing the 0.0131-deg mouth movement from contaminating this comparison.

Interpretation gate: the 8-mm signed residual must remain stable under BOTH controls before we treat its growth with a0 as evidence of a finite-length effect. If it stabilizes, the next experiment is to recompute signed dKII/dtheta near theta=0 separately at each selected length and solve verified local-symmetry roots. If it does not stabilize, investigate the hole polygon, the fixed mouth-width parameter and local FE tip refinement first. MATLAB execution of Step 23 is pending.


## Step 23 results (2026-10-02): finite-length Mode II is not yet numerically converged

Step 23 used the same fine Stage-I initiation mouth at phi*=-1.5606127781 deg, position [0.19998887220,-0.020817033905] m, and studied a0=4 and 8 mm. The existing Step-22 results with NArc=160 and nominal fine FEM h=0.00039269908 m were reused; additional NArc=320/480 full-domain cracked solves were computed on the same nominal fine mesh, and NArc=480 was also computed with nominal coarse h=0.00078539816 m. Signed FE-nodal EDI was evaluated at r_outer/a0=[0.50,0.65,0.80].

At the reference EDI radius r_outer/a0=0.65:
| a0 | fine NArc=160 | fine NArc=320 | fine NArc=480 | coarse NArc=480 |
|---:|---:|---:|---:|---:|
| 4 mm | +1.6721e-5 | -6.7491e-6 | -1.8726e-5 | +3.8870e-5 |
| 8 mm | +9.1118e-5 | +1.26555e-4 | +7.80814e-5 | -8.49684e-4 |

Entries are signed KII/KI. At 4 mm, polygon changes reverse the sign of the tiny residual; thus no physical nonzero kink is resolved. At 8 mm, the fine-mesh residual is positive at all three polygon resolutions but NONMONOTONIC (range ~4.85e-5 at the reference EDI radius), and fixed-NArc coarse versus fine FE changes KII/KI by ~9.28e-4 including sign reversal. This is much larger than the putative finite-length signal. Although the fine NArc=480, 8-mm EDI-domain spread is only ~1.29e-6 (r_out/a0 range [0.50,0.80]), that establishes only integration-domain consistency on ONE displacement field; systematic mesh/geometry errors can be common to all radii.

The matched NArc=480 tip-edge medians are 0.0004215 m coarse and 0.0002146 m fine for 8 mm, so the FE discretizations genuinely differ locally. However, Step 23 also changed global Hmax together with nominal Hmin; its coarse/fine comparison was not isolated tip refinement. K_I differs by about 0.3% while K_II changes much more; signed K_II is especially sensitive near zero.

The command completed every FE solve and printed both result tables, but a MATLAB fprintf multiline text concatenation without continuation ellipses raised a vertcat error after the computations and BEFORE Out was assigned. Consequently no O23 exists in the workspace despite the completed computation. The fprintf source is fixed in commit c814391b; the error was cosmetic and did not affect the printed numerical data. Do not demand a full rerun solely for that formatting problem.

## Step 24 prepared: focused fixed-polygon, fixed-global-mesh 8-mm convergence

The dedicated script verification/sif_audit/main_step24_fixed_geometry_FE_refinement.m accepts O22 (not O23) and defaults to a0=8 mm, NArc=480, RefineFactors=[1,2]. It keeps the fine Step-22 Stage-I mouth/normal, a0, retained-hole polygon, temporary appendix mouth half-shift w=0.0001 m, global Hmax, and Hgrad unchanged, then changes only nominal Hmin/Hhole/Hcrack to half. Both solved meshes use signed FE-nodal EDI at the existing three radii, with measured actual tip-edge medians and EDI-domain spreads. The factor-1 run should independently reproduce the printed Step-23 fine NArc=480, 8-mm KII/KI ~+7.808e-5 at reference r_outer/a0=0.65. One further factor-2 solve tests whether the signal remains stable under controlled local refinement. If needed, rerun with RefineFactors=[1,2,4] to add finer refinement, or isolate the appendix mouth-width dependence after this gate.

IMPORTANT: even if the 8-mm fixed-polygon FE ratios converge, the NArc nonmonotonicity remains unresolved; geometric and mouth-width sensitivity must be audited independently before interpreting a finite-length kink as physical.


## Step 24 result: refining only the local FE scale causes severe EDI-radius sensitivity

Step 24 completed locally with a0=0.008 m, the same fine Stage-I mouth at phi*=-1.560612778 deg, the same retained-hole polygon NArc=480, the same appendix mouth half-shift 0.0001 m, and the same global Hmax=0.007853982 m. Only nominal Hmin/hhole/hcrack changed by a factor of two. Measured median tip-adjacent T3 edge lengths decreased from 2.14597e-4 m to 1.06855e-4 m; the global T6 node count changed from 16491 to 17317 (only ~5% growth). Thus this test substantially refines only the immediate tip region, while much of the EDI annulus may remain at the former mesh scale.

For r_outer/a0=[0.50,0.65,0.80], the unrefined and refined KII/KI values were:
| r_outer/a0 | htip=0.2146 mm | htip=0.1069 mm |
|---:|---:|---:|
| 0.50 | +7.83745e-5 | -8.42691e-7 |
| 0.65 | +7.80814e-5 | +1.12955e-5 |
| 0.80 | +7.70819e-5 | +1.26501e-4 |

The EDI-domain KII/KI spread INCREASED from 1.2926e-6 to 1.2734e-4 (~98.5-fold). The refined K_I values also become radius-sensitive: [0.4373923,0.4374202,0.4374963], about 2.4e-4 relative across radii, whereas the coarser local mesh gave K_I approximately 0.43761 at all three radii. Thus the fine-mesh near-zero signed Mode II is not path-independent, and no physical 8-mm crack kink may be inferred from it.

A possible explanation is that the local Hmin change refined only the immediate tip, while the outer EDI annuli encountered an insufficiently refined or differently remeshed displacement field. A second possibility is sensitivity to FE-interpolated radial q and the varying r_inner rule. Neither explanation has been established; one should NOT ascribe the discrepancy exclusively to the extractor before testing its weighting implementation on the SAME solved field. The FE-nodal versus analytic-radial comparison is diagnostic (the analytic discontinuous gradient may itself have quadrature error), not an automatic ranking of implementations.

## Step 25 prepared: same-field q-gradient and radial-domain audit

The next minimal experiment is verification/sif_audit/main_step25_EDI_weight_audit.m. It independently solves ONE full-domain Stage-II mesh at the exact Step-24 fine local factor=2, a0=8 mm, NArc=480 and fixed global Hmax. It then runs the *same* displacement solution through the EDI routine with both FE-nodal q and analytical radial q-gradient, each under (A) the adaptive inner radius max(0.1 r_outer,2h_tip) and (B) a common physical inner radius 0.0008 m for all three EDI outer radii. For each setting it reports K_I, signed K_II, their ratio, Gauss-point/element participation and auxiliary-field strain consistency, and summarizes the domain spreads. It stores the solved mesh, displacements and material for later postprocessing without further FEM solves.

Step-25 interpretation: if q-method differences explain the large outer-radius spread, audit the q construction, shape gradients and discontinuous radial integration. If both q variants agree but retain the anomaly, prioritize actual FE displacement/equilibrium and spatial annulus-resolution checks. If a common r_inner resolves the instability, the adaptive annulus was a confounder. Run the code in MATLAB before drawing any such causal conclusion.


## Step 25 result: q implementation matters, but r_inner variation is not the cause

Step 25 completed locally for a0=8 mm, NArc=480, fixed fine Stage-I mouth phi*=-1.560612778 deg, unchanged global Hmax, and Step-24 local-refinement factor=2. The new solved mesh reproduced Step 24 **exactly**, with 17,317 T6 nodes, 8,317 T3 elements and median tip-edge length 1.06855353e-4 m. The same displacement field was postprocessed under all four EDI settings: FE-nodal versus analytic radial q-gradient, and adaptive versus fixed physical inner radius r_inner=0.0008 m.

The primary FE-nodal q results at r_outer/a0=[0.50,0.65,0.80] were [-8.42691e-7,+1.12955e-5,+1.26501e-4] for KII/KI and [0.437392345,0.437420227,0.437496280] for KI. With common inner radius, KII/KI changed only to [+2.53545e-6,+1.37993e-5,+1.30252e-4], with virtually unchanged domain spread (1.2734e-4 versus 1.2772e-4). Therefore the variable inner-radius rule is NOT the main explanation.

The analytic-radial gradient with abrupt GP inclusion at r_inner/r_outer performed much more poorly under the present seven-point T6 quadrature: adaptive KII/KI=[-4.57462e-4,+1.34966e-3,-9.19102e-3], with domain spread ~1.0541e-2; common-inner results were similarly unstable. These numbers do NOT justify rejecting the mathematics of the radial weight function. The current analytic-radial code cuts off integration at GP locations without an element-conforming integration submesh, making integration error across sharp circular q-gradient cutoffs a plausible explanation. Treat the FE-nodal q route as the working reference while auditing it.

Internal auxiliary stress/strain consistency is strong: median relative mismatch of order 1e-9 or smaller for both auxiliary modes in both q settings. The FE-nodal maximal mismatch reached order 2e-5 at some near-branch-cut GP locations. This tests local auxiliary field implementation, not complete EDI path independence. The most important unanswered question remains whether the FE-nodal residual arises in EDI interpolation/quadrature on this nonuniform T6 mesh or in the actual FEM displacement/equilibrium field.

## Step 26 prepared: exact Williams field replay on the identical Step-25 mesh

The driver verification/sif_audit/main_step26_EDI_exact_field_replay.m reuses O25.mesh, O25.mat, O25.crack and O25.U and performs **no additional FEM solve**. It constructs exact leading-order unit mode-I and mode-II Williams displacement fields at the same T6 nodes, with duplicate upper/lower crack-face T3 nodes and their T6 midsides explicitly classified by the original collapsed mesh topology, including tip-adjacent face midsides. The exact local displacements are rotated into the actual crack's global reference frame, interleaved into the same DOF convention, and passed to the *unchanged* signed FE-nodal EDI extractor.

For each of the three EDI outer radii, it uses the same fixed r_inner=0.0008 m and records the 2x2 recovery matrix M. Ideally its first column is [1,0] for unit pure mode I and second [0,1] for unit pure mode II; off-diagonal entries quantify artificial modal leakage on this exact T6 mesh. It separately recovers K from the actual O25.U at identical domains. A diagnostic M-inverse(actual_K) is printed only when M is well-conditioned and is NOT automatically an acceptable physical correction.

If the recovered exact-field matrix is nearly identity and domain-stable, the EDI machinery can reproduce exact singular fields on this mesh, shifting attention to FEM solution accuracy, tip-zone/collapsed-geometry artifacts, FE equilibrium, or missing nonsingular terms. If the exact-field recovery matrix is inaccurate or radius-sensitive, EDI quadrature/interpolation and crack-face labeling require correction before any further physical interpretation of residual KII. This is more informative than another expensive geometric or crack-angle sweep. Step 26 MATLAB execution is pending.


## Step 26 result: exact-field recovery is good in relative terms but not accurate enough to resolve near-zero physical KII

The same-mesh exact Williams replay ran successfully on the Step-25 NArc=480, locally refined a0=8 mm T6 mesh (17,317 nodes). Topologically classified crack-face T6 nodes: 76 upper, 76 lower, zero unclassified on the negative-axis crack line. With FE-nodal q and fixed physical r_inner=0.0008 m, recovery of the exact [KI=1,KII=0] and [KI=0,KII=1] unit fields at r_outer/a0=[0.50,0.65,0.80] produced identity-matrix errors [4.073e-4,6.681e-4,4.460e-4] and well-conditioned matrices (cond ~1.0003).

Crucially, exact Mode-I artificial KII leakages were [+1.23857e-4,+2.20224e-5,+2.76239e-4] at these radii. The actual KII measured on precisely the same FEM mesh was [+1.10899e-6,+6.03616e-6,+5.69853e-5], with KI ~0.4374. Therefore the artificial mixed-mode contribution of a large KI is comparable to or larger than the actual KII signal. The diagnostic M-inverse applied to actual [KI,KII] gives KII=[-5.30917e-5,-3.60204e-6,-6.39165e-5], which is also strongly integration-domain dependent. This calibration is NOT accepted as a physical correction because the FE solution contains nonsingular/higher-order terms absent from the synthetic leading-order Williams basis, and a field-dependent error need not follow the same recovery matrix.

The exact-field consistency test thus validates gross EDI normalization and crack-face sign bookkeeping on the given T6 mesh but does NOT reach the precision needed to infer whether theta* departs from the normal by a few thousandths of a degree. Moreover, an exact field sampled at mesh nodes incurs T6 interpolation error near the singularity; the observed leakage may be caused by interpolation, Gauss quadrature, mesh asymmetry, or combinations of these. No single cause is isolated yet.

## Step 27 prepared: compare Dunavant quadrature order on precisely the same T6 fields

The signed SIF_LEFM_interaction_EDI.m now accepts optional 'QuadratureRule' 7 (unchanged historical default), 12 (degree-6 Dunavant) or 16 (degree-8 Dunavant). The 12/16 rules were checked against all reference-triangle polynomial moments up to their nominal degree, with maximum coefficient-level residuals near machine precision. Historical production calls are unchanged unless they explicitly select a different rule.

The new verification/sif_audit/main_step27_EDI_quadrature_replay.m accepts the existing O25 and O26: no new geometry, mesh, solve, face reconstruction, or Stage-I optimization is needed. The default [7,16] quadrature comparison uses the exact Mode-I and Mode-II sampled fields from Step 26 and the same actual O25 displacement field. It additionally creates a manufactured EXACT AFFINE T-stress field with local stress [1,0,0], which is traction-free on ideal crack faces and representable exactly by T6 shape functions. The affine unit stress is an analytic diagnostic, NOT a measured T-stress from the actual boundary-value problem. For all three fixed-r_inner=0.0008 m EDI outer domains, the driver reports exact modal recovery matrices, actual KI/KII, artificial mode-I leakage, and spurious SIFs under the affine T-stress field.

The historical 7-point exact and actual results can be reused from O26 where settings match; the T-stress field and all 16-point evaluations are new postprocessing computations on O25.mesh/O25.U. If 16-point integration significantly improves both manufactured-field recovery and the real-field integration-radius consistency, quadrature was an important cause. If synthetic Williams leakage improves without improving real-field domain spread, the real computed displacement field, nonsingular terms, or the FE-nodal q discretization remains suspect. If neither improves, investigate nodal interpolation of singular fields, q-gradients and full-annulus spatial FE resolution rather than spending additional runs solely on integration points.


## Step 27 result: smooth-field quadrature improves, Williams leakage barely changes

The existing Step-25 8-mm, NArc=480, locally factor-2 refined T6 mesh (17,317 nodes) was replayed using historical seven-point and degree-8 sixteen-point Dunavant EDI quadrature. The same FE-nodal q and fixed physical r_inner=0.0008 m were retained. No new FEM solve was performed.

| Metric | 7 points | 16 points |
|---|---:|---:|
| Maximum unit-Williams 2x2 recovery matrix error ||M-I||F | 6.6813e-4 | 6.6932e-4 |
| Maximum |unit Mode-I leakage into KII/KI| | 2.7632e-4 | 2.6657e-4 |
| Exact Mode-I leakage domain spread | 2.5429e-4 | 2.4924e-4 |
| Actual FEM KII/KI domain spread | 1.2772e-4 | 1.1970e-4 |
| Max spurious KII per unit manufactured affine T stress | 6.3688e-7 | 9.5543e-9 |

At reference r_outer/a0=0.65, the actual ratio changes from +1.37993e-5 to +9.75221e-6, but the marked radius dependence persists. For the manufactured nonsingular T stress, stronger quadrature reduces the spurious KII by about a factor of 67, confirming that the new integration order has a real numerical effect. However, the singular Williams recovery error and artificial modal leakage remain comparable: the principal near-zero Mode-II problem cannot be resolved merely by increasing Gauss-point order on the present mesh.

Interpretation limit: the nodally sampled exact Williams fields are interpolated by T6 elements. Their computed gradients are not the exact analytical gradients. Therefore nearly quadrature-insensitive artificial modal leakage could originate in nodal interpolation of the singular fields; it could also originate in the FE-nodal q construction or other EDI errors. The test has NOT yet isolated the cause, and no numerical value from these experiments establishes a finite physical kink angle.

## Step 28 prepared: exact Williams gradients at the same Gauss points

SIF_LEFM_interaction_EDI.m now has an opt-in verification-only 'AnalyticActualK'=[KI,KII] argument, default empty. With empty input, all production FE displacement-gradient extraction behaves exactly as before. When enabled, the 'actual' field in the EDI density is evaluated directly from the exact Williams stress/strain and derivative formulas at each Gauss point, while retaining the SAME T6 mesh coordinates, FE-nodal q-gradient, Gauss quadrature, crack-tip frame, and EDI annulus. This bypasses only T6 nodal interpolation of the manufactured singular field and must NEVER be used to postprocess an actual numerical FEM field.

The verification driver 'verification/sif_audit/main_step28_EDI_interpolation_audit.m' accepts O25, O26 and O27. It reuses O27's existing nodally sampled exact-field recovery matrices and actual-field results at the historical 7- and 16-point rules, and performs 12 new extraction-only operations (two unit Williams modes, three radii, two quadrature orders) with exact gradients at Gauss points. It reports both recovery matrices, errors from identity and artificial Mode-I leakage. No new FE mesh or solve is performed. If exact-Gauss-point recovery improves strongly compared with nodally sampled Williams recovery, the test isolates substantial T6 interpolation error in the MANUFACTURED singular field. If exact-Gauss-point recovery still shows significant leakage, investigate q interpolation/integration or EDI formulation. Either result by itself does NOT establish the cause of error in the actual FEM solution, which can include mesh and equilibrium errors.


## Step 28 result: T6 interpolation of nodal Williams fields is the dominant synthetic-field recovery error

Step 28 successfully replayed both EXACT unit Williams fields two ways on the same Step-25 17,317-node T6 mesh, FE-nodal q, r_inner=0.0008 m, and three outer annuli r_outer/a0=[0.50,0.65,0.80]. Its first path reused nodally sampled unit modes from Step 26/27; its second bypassed ONLY the manufactured field's T6 interpolation by evaluating exact singular stress, strain and displacement gradient directly at each Gauss point (new verification-only `AnalyticActualK` EDI option). No FEM mesh/solve was performed.

| Quadrature | Max nodal Williams recovery ||M-I||F | Max exact-GP recovery ||M-I||F | Max nodal unit Mode-I artificial KII/KI | Max exact-GP unit Mode-I artificial KII/KI |
|---:|---:|---:|---:|---:|
| 7 points | 6.6813e-4 | 3.5607e-5 | 2.7632e-4 | 2.4520e-5 |
| 16 points | 6.6932e-4 | 7.6912e-7 | 2.6657e-4 | 4.5720e-7 |

With 16-point quadrature, direct Gauss-point evaluation reduces the worst 2x2 recovery-matrix error by about 870-fold, while reducing the worst artificial Mode-I KII leakage by about 583-fold. The exact-Gauss-point modal leakage changes from +5.0772e-8 to -4.5720e-7 across the three annuli; the nodal-Williams leakage changes from +1.2339e-4 to +2.6657e-4. Thus interpolation of the SYNTHETIC singular field by quadratic shape functions, not the fundamental interaction-integral normalization, is the leading source of the observed manufactured-field error at the tested mesh sizes. FE-nodal q/integration residuals are small when exact fields are supplied at 16-point Gauss locations.

The ACTUAL computed FEM displacement field is still domain-sensitive under the same FE-nodal q and 16-point rule: KII/KI=[+2.19786e-6,+9.75221e-6,+1.21895e-4], domain spread 1.1970e-4. One cannot infer automatically that its error is solely T6 interpolation: the FE displacement solution and its near-tip gradient may additionally contain errors from mesh topology, stress equilibrium, and collapsed crack-face geometry, as well as valid higher-order physical terms. The tiny nonzero signed KII remains unresolved physically.

## Step 29 prepared: independent crack-face displacement jump on the SAME FEM solution

The new driver `verification/sif_audit/main_step29_crack_face_COD_audit.m` accepts the existing O25, O26 and O27 and **performs no FEM solves**. It uses the authoritative O26 upper/lower T6 crack-face classification (including mid-edge nodes), extracts the actual O25 displacement field and the two nodally sampled exact Williams fields on each side, and separately interpolates upper and lower local displacement components against the exact physical distance r BEHIND the crack tip. For a straight traction-free crack, the asymptotic jump formulas give

- `KI_app(r) = mu/(kappa+1) * sqrt(2*pi/r) * [u_normal^upper(r) - u_normal^lower(r)]`,
- `KII_app(r) = mu/(kappa+1) * sqrt(2*pi/r) * [u_tangent^upper(r) - u_tangent^lower(r)]`,

using signs consistent with the existing EDI auxiliary displacements and checked independently with O26's exact unit modes. The script samples r/a0=[0.04,0.40] at 61 points and fits `K_app(r)=K0+b*r/a0` over four windows extending to r/a0=0.12,0.20,0.30,0.40. The fit intercepts provide near-tip COD SIF estimates, reported alongside each fit's exact-unit-mode 2x2 synthetic recovery matrix and the reference 16-point EDI values.

This is an independent EXTRACTION METHOD, not an independent FEM solution. COD's own interpolation and finite-window extrapolation uncertainties must be assessed through the exact-unit synthetic matrix and window convergence. If real-FEM COD KII is stable across windows and consistent with a small near-normal direction while EDI remains strongly radius-dependent, revisit EDI interactions with higher-order FE fields and annulus-region residuals. If COD KII is also window-sensitive, the actual FE crack-tip displacement field is not resolved at the small signed KII scale, and more physically directed mesh/geometry validation is required. No calibration inverse should be presented as a production physical KII without separate justification. MATLAB runtime testing of Step 29 remains pending.
