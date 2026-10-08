# Centered-hole symmetric path-stability control

## Purpose

This control checks the sign convention and local stability of the same
EDI--MTS directional recurrence used for the asymmetric circular-hole
trajectory.

The experiment is deliberately local: it tests the one-step map
`theta_2 -> theta_3` near the exactly symmetric straight configuration.
It is not an increment-size convergence study and is not used as a second
main physical example.

## Geometry and recurrence

- plate: `A = 0.300 m`, `B = 0.100 m`;
- centered circular hole: center `(0.150,0) m`, radius `R = 0.030 m`;
- crack mouth: `P0 = (0.180,0) m`;
- prescribed first tip: `P1 = (0.184,0) m`;
- crack increment: `Delta a = 4 mm`;
- probe directions:
  `theta_2 = {-0.10,-0.05,0,+0.05,+0.10} deg`.

Each two-segment probe was rebuilt and physically solved independently with
the same paired-core, EDI, linear-solver, and MTS machinery used for the
asymmetric calculation.

Define

```text
F(theta_2) = theta_3
           = theta_2 + alpha_MTS(KI(theta_2),KII(theta_2)).
```

## Straight symmetric state

At `theta_2 = 0`:

- `KI = 0.425410776939 MPa*sqrt(m)`;
- `KII/KI = 6.224e-9`;
- `theta_3 = -7.13e-7 deg`.

Thus the straight symmetry line is recovered as a numerical mode-I fixed
path to the accuracy relevant for the trajectory calculation.

## Five-probe map

| theta_2 (deg) | theta_3 = F(theta_2) (deg) |
|---:|---:|
| -0.10 | +0.011415 |
| -0.05 | +0.005708 |
| 0 | -0.000000713 |
| +0.05 | -0.005708 |
| +0.10 | -0.011415 |

The recovered `KI` is even with respect to the imposed perturbation, while
`KII/KI` and `theta_3` are odd to numerical mismatches of order `1e-9`.

Pairwise centered slopes are:

- from the `+/-0.05 deg` probes: `-0.1141583174`;
- from the `+/-0.10 deg` probes: `-0.1141509209`.

The least-squares origin fit gives

```text
F'(0) = -0.1141524002.
```

## Interpretation

Because `abs(F'(0)) < 1`, the straight symmetric path is locally attracting
for this one-step discrete recurrence.  The negative multiplier means that a
small angular perturbation changes sign on the next update while its
magnitude is reduced by approximately a factor of 8.8.

This control supports two conclusions used in the manuscript:

1. the signed EDI/MTS implementation respects the expected reflection
   symmetry and mode-I fixed path;
2. the curved asymmetric trajectory is not generic numerical drift away from
   a straight crack path.

The result is only a local stability statement for the chosen `4 mm`
increment.  It does not establish global stability or crack-increment
convergence.
