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
