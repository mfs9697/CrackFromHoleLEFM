# Centered-hole symmetry benchmark provenance

Physical run date: 2026-10-07.

- Computational source branch: `symmetric-path-stability`
- Computational source commit: `6945371dfdfac441b90970690b7e1728c7f06259`
- MATLAB command: `R = main_symmetric_path_stability('AllowPhysicalSolves',true);`
- Probe set: `[-0.10 -0.05 0 0.05 0.10]` degrees
- First segment: 4 mm, exactly along the centered-hole symmetry axis
- Second segment: 4 mm, prescribed at each probe angle
- One physical solve per probe; no third segment generated
- Uploaded complete console log SHA256:
  `f4dd0397a86fb1ad395270a496b1eac209effe450dddc8dc3603ac1852695301`
- Readable numerical export: `paper/data/symmetric_stability_benchmark.csv`

All five qualification and physical-solve gates passed. The straight
control gives `KII/KI = 6.22351746186701e-09` and
`theta3 = -7.1316256897125e-07 deg`. Pairwise map multipliers are
`-0.114158317407462` for 0.05 deg and `-0.114150920881333` for
0.10 deg. The four nonzero probes give the origin-constrained fit
`F'(0) = -0.114152400186559`.

This benchmark is separate from the asymmetric production trajectory. It
tests symmetry preservation, sign conventions, SIF extraction, and local
behavior of the implemented EDI/MTS direction update. It does not establish
4/2/1-mm trajectory convergence and does not isolate the free boundary as
the cause of the late mode-mixity reversal.
