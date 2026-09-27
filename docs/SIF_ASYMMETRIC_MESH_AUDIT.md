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
