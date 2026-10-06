# Symmetric crack-path stability experiment

## Purpose

This experiment separates two questions that are identical only for an
unperturbed symmetric crack:

1. Is the straight crack on the symmetry axis an exact mode-I solution?
2. Is that straight solution locally stable when the crack itself is allowed
   to break the symmetry?

The external problem is symmetric, but a kinked crack is not. Therefore the
stability experiment must use the **full plate**. A half-domain model is valid
for the straight control only and would constrain the nonzero angular probes.

## Geometry and loading

The plate is

\[
0\le x\le A,\qquad -B\le y\le B,
\]

with

\[
A=0.300\ \mathrm{m},\qquad B=0.100\ \mathrm{m}.
\]

The circular hole is centered at

\[
(x_c,y_c)=(0.150,0)\ \mathrm{m},\qquad R=0.030\ \mathrm{m}.
\]

The right-hand initiation point and local frame are prescribed exactly by
symmetry:

\[
P_0=(0.180,0)\ \mathrm{m},\qquad
\mathbf n=(1,0),\qquad
\mathbf t=(0,1).
\]

Remote vertical unit tension and the existing minimal anchoring are retained.
The crack increment is the already-audited value

\[
\Delta a=4\ \mathrm{mm}.
\]

No Stage-I FEM solve is needed to locate \(P_0\).

## Perturbation experiment

The first segment is fixed:

\[
P_1=P_0+\Delta a\,\mathbf n,
\qquad \theta_1=0.
\]

For each prescribed second-segment perturbation

\[
\theta_2\in
\{-0.1^\circ,-0.05^\circ,0,
  +0.05^\circ,+0.1^\circ\},
\]

construct

\[
P_2(\theta_2)
=
P_1+\Delta a
\left(
\cos\theta_2\,\mathbf n+
\sin\theta_2\,\mathbf t
\right).
\]

Each two-segment path is independently qualified with the existing structured
new-tip core, native COD sampling, interaction EDI, mesh gates, and prescribed
Williams replay. Exactly one physical solve is then performed at \(P_2\).
The current geometry is never searched over. MTS predicts the next absolute
direction \(\theta_3\), but no third segment is generated.

The exactly straight two-segment control requires the explicit
`AllowStraightPath` qualifier option. Its default is false, so the production
incremental workflow retains its original kink requirement.

## Symmetry identities

For exact reflection symmetry about \(y=0\),

\[
K_I(+\theta)=K_I(-\theta),
\qquad
K_{II}(+\theta)=-K_{II}(-\theta).
\]

The MTS response should therefore satisfy

\[
\Delta\theta_3(+\theta)
=
-\Delta\theta_3(-\theta),
\qquad
\theta_3(+\theta)
=
-\theta_3(-\theta).
\]

The straight control should give

\[
K_{II}(0)=0,
\qquad
\Delta\theta_3(0)=0,
\qquad
\theta_3(0)=0
\]

up to numerical symmetry/extraction error.

Because the full-domain mesh is generated numerically rather than constructed
by exact reflection, the experiment measures any residual numerical
symmetry-breaking instead of assuming it away.

## Stability map

Define the one-step direction map

\[
F(\theta_2)=\theta_3.
\]

The straight path is a fixed point if \(F(0)=0\). Near zero,

\[
F(\theta)\simeq m\theta,
\qquad
m=F'(0).
\]

The driver estimates \(m\) from the symmetric finite probes, both pairwise and
by a least-squares fit through the origin.

Interpretation:

- \(|m|<1\): a small angular perturbation decays;
- \(m<0\) and \(|m|<1\): decay with alternating sign;
- \(|m|>1\): the perturbation is amplified and the straight path is locally
  symmetry-breaking;
- \(|m|\approx1\): the finite probe set is insufficient to distinguish
  stability from neutrality.

This classification is a finite-perturbation numerical estimate, not an
analytical bifurcation proof.

## Running

The physical solves are guarded. From the repository root:

```matlab
R = main_symmetric_path_stability( ...
    'AllowPhysicalSolves',true);
```

The default test performs five physical solves. To test a smaller amplitude
after inspecting the first result, for example:

```matlab
Rfine = main_symmetric_path_stability( ...
    'ThetaProbeDeg',[-0.05 -0.025 0 0.025 0.05], ...
    'AllowPhysicalSolves',true, ...
    'OutputDir','verification/crack_path/symmetric_stability_fine');
```

Use a separate output directory for a different probe set.

## Outputs

Default directory:

```text
verification/crack_path/symmetric_stability/
```

The driver writes:

- one candidate, qualification summary, physical checkpoint, and compact
  physical result per probe angle;
- `symmetric_stability_probes.csv`;
- `symmetric_stability_parity.csv`;
- `symmetric_stability_straight_control.csv`;
- `symmetric_stability_linearization.csv`;
- `symmetric_stability_result.mat`;
- four EPS/PNG figures for mode mixity, MTS correction, the one-step stability
  map, and mode-I parity.

The first run should be interpreted before reducing the perturbation amplitude
or propagating either symmetry-broken branch.
