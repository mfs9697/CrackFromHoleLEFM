# MTS signed angle: admissible stationary-root safeguard

This narrow change prepares the established maximum tangential stress
(MTS) helper for later incremental crack-trajectory work. It is
**not** an analytical Williams-field EDI benchmark, does not use or
alter physical FEM displacements or saved source geometry, and does
not simulate a crack trajectory.

## Correct leading-order theory

From
`KI sin(theta)+KII(3 cos(theta)-1)=0`,
setting `t=tan(theta/2)` gives
`2*KII*t^2-KI*t-KII=0`.
The near-forward tensile root for `KI>0` is

`theta=2*atan(-2*KII/(KI+sqrt(KI^2+8*KII^2)))`.

In particular, the numerator is **-2*KII**, not **-4*KII**.
For the qualified original edge-crack L2 signed pair
`KI=0.659352421717419`, `KII=+0.0183612029902817`
MPa sqrt(m), the angle is **-3.185314152245 deg**.

## Legacy implementation problem

The former `kink_angle_LEFM_MTS.m` numerically found one
stationary root, but if it produced negative leading
circumferential stress it tried `-theta`.
Except for special symmetric cases, `-theta` is not
another stationary root for the same signed mixed-mode pair.
Also, returning `theta=0` unconditionally for
`KII≈0` would prescribe zero turn for pure
compressive Mode I.

## Narrow remedy

The function now analytically forms both roots of the
stationarity quadratic. Candidates are retained only if:
- finite and within the original forward ±80-degree sector;
- original stationarity residual within normalized tolerance;
- positive leading hoop stress (tensile);
- negative exact local angular curvature (strict local maximum).

The admissible candidate with largest tensile stress is
returned. If no candidate exists (including pure compression,
zero loading and certain compressive mixed-mode cases),
the function returns a **NaN triplet** and emits a warning
with identifier `kink_angle_LEFM_MTS:NoTensileStationaryRoot`
rather than inventing a crack-growth direction. Calling
code must handle that status and separately enforce
a fracture initiation/growth condition.

The original three-output public MATLAB function signature
is preserved. The existing `kinking_LEFM_1leg.m` call site is also
guarded: if the MTS helper reports NaN or a non-tensile stationary
state, it raises `kinking_LEFM_1leg:NoAdmissibleTensileMTS`
instead of passing an undefined angle downstream.
**No SIF calculation, material model, FEM solve or mesh routine changes.**

## Pure theory regression

`test_kink_angle_LEFM_MTS_roots_only.m`

Checks 13 signed algebraic gates:
true original L2 baseline, pure tensile I,
positive/negative mixed modes, pure positive/negative II
at ±70.5288°, inadmissible pure compressive I,
weak- and strong-shear compressive examples,
true stationarity of every retained direction,
positive tensile hoop maximum/negative curvature,
sector bounds, and zero physical computations.

From repo root in MATLAB R2023a:

```matlab
addpath(genpath(pwd));
T = test_kink_angle_LEFM_MTS_roots_only();
assert(T.pass);
```

**Status:** the source and mathematical branch behavior have been
reviewed, but the new code has not been executed in the author's
MATLAB installation. No new physical results, crack trajectory,
fracture-resistance threshold, or analytical Williams-field benchmark
are claimed. Do not merge until MATLAB tests and any relevant call
sites are checked for the new NaN-on-inadmissible behavior.
