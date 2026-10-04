# SIF audit synthesis and closure

**Status:** scientifically closed; repository integration is recorded below.

**Scope:** validation of the very small positive Mode-II component in the asymmetric cracked-hole LEFM problem, with special attention to mesh-induced parity leakage, extraction bias, linear-solver effects, and mesh convergence.

## Executive conclusion

The audit supports a physical mixed-mode ratio

\[
\boxed{K_{II}/K_I \approx 1.066\times10^{-4}}
\]

for the audited asymmetric cracked-hole problem.

The strongest physical evidence is the three-scale structured C03 family:

| Scale \(s\) | Target tip size (mm) | T3 | T6 nodes | Matched EDI \(K_{II}/K_I\) |
| ---: | ---: | ---: | ---: | ---: |
| 1 | 0.05402465 | 32,980 | 66,854 | \(1.0631\times10^{-4}\) |
| \(1/\sqrt2\) | 0.03820120 | 62,570 | 126,112 | \(1.0646\times10^{-4}\) |
| 1/2 | 0.02701233 | 122,691 | 246,701 | \(1.0652\times10^{-4}\) |

The matched EDI sequence is monotone. The relative increments decrease from approximately +0.1378% to +0.0596%.

Using the recorded audit endpoints gives an **approximate** observed order

\[
p\approx 2.41
\]

and an **approximate** Richardson limit

\[
(K_{II}/K_I)_\infty \approx 1.0657\times10^{-4}.
\]

The fine-grid value differs from this approximate extrapolated value by about 0.05%. This is a numerical truncation indicator, **not** a probabilistic uncertainty bound.

All eight COD estimators are also monotone over the same three geometrically spaced scales, and every second refinement increment is smaller than the first. Their observed orders are approximately 1.37--2.10. The quadratic COD extrapolations cluster around \(1.064\times10^{-4}\) to \(1.065\times10^{-4}\), close to the EDI estimate.

## What the audit was trying to decide

Historical asymmetric calculations repeatedly produced a very small positive ratio of order \(10^{-4}\). The central question was not whether two numerical runs happened to agree, but whether this signal could instead be created by:

1. mesh asymmetry;
2. interpolation of the singular field on ordinary T6 elements;
3. EDI mode leakage;
4. COD fitting-window effects;
5. an inadequately controlled refinement path;
6. or linear-solver artifacts.

The audit therefore isolated these mechanisms one at a time rather than immediately refining the physical problem.

## Evidence chain

| Audit stage | Question | Main result | Evidential role |
| --- | --- | --- | --- |
| Historical Steps 34/38 | Does the original asymmetric calculation repeat under local refinement? | \(K_{II}/K_I\approx1.0754\times10^{-4}\) and \(1.0796\times10^{-4}\), about 0.39% apart | Motivation only; meshes were not a controlled family |
| Step44 | Can COD/EDI recover prescribed Williams fields? | COD essentially exact; EDI pure-I to II leakage about 0.0235% of the tiny mixed signal | Extractor qualification, not physical validation |
| Steps45--51 | Why does a symmetric FEM problem show nonzero/sign-changing \(K_{II}\)? | Exact nodal pure-I leakage changes sign with mesh; direct Gauss-point replay reduces leakage by orders of magnitude | Identified singular-field interpolation/mesh representation as a contamination mechanism |
| Step49 | Are the historical symmetric meshes reflection paired? | No | Explains why symmetry alone did not guarantee zero discrete Mode II |
| Steps52--56 | Can an exactly reflection-paired mesh be constructed? | Yes | Controlled mesh design |
| Steps57--59 | Does a symmetric physical FEM solution on that paired mesh still leak Mode II? | Actual EDI \(K_{II}/K_I\approx-4.69\times10^{-14}\); exact-nodal pure-I about \(-1.09\times10^{-15}\) | Strong validation that exact reflection pairing removes the spurious symmetric Mode-II mechanism |
| Step60 | Is the old asymmetric Step38 mesh a clean refinement family? | No; support intersects an irregular patchwork mesh | Justifies replacing ad hoc refinement |
| Step62 | Can a deterministic structured paired-core family be built? | Yes; exact pairing and synthetic mixed recovery | Family construction |
| Step62B | Can the exterior be calibrated without using physical SIFs? | C03 selected with max adjacent-size ratio 1.79678451 | Avoids result-driven mesh tuning |
| Step63/63R | Is the Level-0 physical COD signal reproducible? | Replacement solve reproduced the first Step63 COD fingerprint essentially exactly | Physical field reproducibility |
| Step64 | Do COD and EDI agree on the same Level-0 field? | EDI \(1.0631\times10^{-4}\); closest COD fit within about 0.46% | Cross-extractor evidence |
| Step65 | Is Level 1 a clean \(h/2\) family member? | Yes; exact scale halving, all mesh/synthetic gates passed | Controlled refinement |
| Step66 | Is the direct solver safe for Level 1? | Factorization risk high; assembly itself not limiting | Numerical linear-algebra audit |
| Step67 | Does fixed ICT preconditioning work? | No; nonpositive pivot before PCG | Negative result preserved; no physical solve consumed |
| Step67A | Does SGS-PCG reproduce the direct Level-0 solution? | COD and EDI reproduced far inside 0.05% gates | Solver qualification |
| Step68 | What changes under \(s:1\to1/2\)? | EDI ratio +0.201%; eight COD ratios +0.32% to +0.52% | Two-level physical convergence evidence |
| Step69 | Does an intermediate \(s=1/\sqrt2\) lie on a smooth trend? | Yes; all COD and EDI sequences monotone with decreasing increments | Three-level convergence evidence |

## Why the final Mode-II signal is accepted as physical

### 1. The known spurious symmetry mechanism was isolated and removed

On non-paired symmetric meshes, prescribed exact nodal pure-I fields produced mesh-dependent artificial \(K_{II}\), including sign changes. Evaluating the same analytical field directly at Gauss points reduced this leakage drastically.

When the mesh was made exactly reflection paired, both the actual symmetric FEM field and exact-nodal pure-I replay produced \(K_{II}/K_I\) at roundoff scale.

Therefore the historical symmetric leakage mechanism is real, understood, and absent from the final C03 paired-core design.

### 2. Prescribed mixed fields recover the \(10^{-4}\) signal

At every qualified C03 scale, pure-I, pure-II, and tiny-mixed prescribed-field controls pass. The synthetic \(K_{II}=10^{-4}\) component is recovered to essentially machine precision relative to the physical scale of interest.

These controls establish extractor/topology capability. They do not by themselves validate the physical FEM solution.

### 3. Independent physical extractors agree

At Level 0 the matched EDI ratio is \(1.0631\times10^{-4}\). The COD fits form a narrow positive band of the same order, with quadratic fits closest to EDI.

At Level 1 the matched EDI ratio is \(1.0652\times10^{-4}\). COD mean and median remain within about 0.82% and 0.55% of EDI, respectively.

The two extraction methods estimate different finite-distance quantities before extrapolation, so exact equality is not expected.

### 4. The physical signal converges in a controlled mesh family

The family changes only the global size scale

\[
h(r;s)=s\,[h_{\rm base}+0.028\,r],
\]

while keeping the physical geometry, 6-mm reflection-paired core, C03 exterior grading, support, solver protocol, and EDI annulus fixed.

The sequence \(s=1,\ 1/\sqrt2,\ 1/2\) is geometrically spaced with constant refinement ratio \(\sqrt2\). The physical EDI and all eight COD estimates are monotone with shrinking increments.

This is the strongest evidence in the audit.

## Final numerical recommendation

For reporting the audited physical result, use

\[
\boxed{K_{II}/K_I \approx 1.066\times10^{-4}}.
\]

If an extrapolated value is useful, state it only as

\[
(K_{II}/K_I)_\infty \approx 1.0657\times10^{-4},
\qquad p\approx2.4,
\]

with the qualifier **approximate three-level estimate**.

Do not report the extra digits as formal precision because the exact Step67A/Step68 compact endpoint files did not survive branch switching; Step69 therefore used recorded rounded/reconstructed audit fingerprints for at least part of the endpoint data.

A conservative manuscript-level formulation is:

> The matched interaction-integral estimate converges monotonically over three geometrically spaced structured mesh levels to \(K_{II}/K_I\approx1.066\times10^{-4}\). The refinement increments decrease from about 0.14% to 0.06%; an approximate three-level extrapolation gives \(1.0657\times10^{-4}\). Independent COD estimates exhibit the same monotone trend.

## Recommended production mesh policy

1. Use a reflection-paired structured crack-tip/core mesh whenever symmetry leakage at the \(10^{-4}\) level matters.
2. Keep the entire primary EDI support inside the paired region. In the audited problem the fixed annulus is \(r_{\rm inner}=0.8\) mm and \(r_{\rm outer}=5.2\) mm, inside the 6-mm paired core.
3. Use the C03 exterior calibration as the reference design:
   - transition length 8 mm;
   - far-field slope 0.10;
   - boundary-metric growth 0.25;
   - neighbor-ratio target 1.8.
4. Use the continuous global scale \(s\) as the primary convergence parameter. Do not claim convergence from unrelated changes in grading, paired radius, or mesher settings.
5. Treat paired radius, structured slope, transition length, and other grading parameters as robustness/sensitivity variables, not as substitutes for \(h\)-convergence.
6. Run pure-I, pure-II, and a tiny mixed prescribed-field extraction control on each new mesh family member before any physical interpretation.
7. Prefer the free-DOF SPD formulation with symamd plus parameter-free SGS-PCG for large meshes. It was qualified against the direct Level-0 solution and avoided the direct-factorization memory risk.
8. Checkpoint the physical displacement field immediately after the linear solve, before COD/EDI postprocessing.
9. Preserve compact numerical results outside branch-fragile generated paths or commit a small text/CSV summary. Do not rely on untracked MAT files as the sole record of endpoint values.

## Recommended extraction and reporting policy

Use the interaction EDI as the primary SIF estimator with the audited fixed settings:

- FE-nodal q;
- 16-point quadrature;
- \(r_{\rm inner}=0.8\) mm;
- \(r_{\rm outer}=5.2\) mm.

Use native crack-face COD as an independent secondary diagnostic with the fixed four windows:

- \(0.04\le r/a_0\le0.20\);
- \(0.04\le r/a_0\le0.30\);
- \(0.08\le r/a_0\le0.30\);
- \(0.12\le r/a_0\le0.30\);

and polynomial degrees 1 and 2.

Do not select a COD window after looking for the value closest to EDI. Report the predeclared family of fits. Quadratic fits currently cluster most closely around the EDI extrapolation, but that is an observed property, not a tuning rule.

## What is closed

The following questions are considered closed for the audited configuration:

- whether the \(O(10^{-4})\) positive physical Mode-II signal is merely symmetric-mesh parity leakage: **no**;
- whether the EDI extractor can represent a \(10^{-4}\) mixed component on the final mesh family: **yes**;
- whether the Level-0 physical result is reproducible: **yes**;
- whether the tiny signal survives controlled paired-mesh refinement: **yes**;
- whether the physical sequence shows a smooth three-scale convergence trend: **yes**;
- whether the memory-safe SGS-PCG formulation reproduces the direct Level-0 solution: **yes**.

## What remains outside the claim

The audit does **not** establish:

- a probabilistic uncertainty interval;
- a formally high-precision Richardson limit;
- a universal error bound for all asymmetric crack geometries;
- that every exterior-mesh parameter is irrelevant;
- superiority of quadratic COD fits in every problem;
- or that the historical unpaired meshes can be corrected a posteriori.

A further physical mesh is not required to support the present conclusion. It would only be justified if a higher-confidence asymptotic-order study or formal error estimator becomes a separate research objective.

## Closure decision

**Scientific verdict: PASS.**

The physical positive Mode-II component is accepted as numerically resolved for the audited problem, at

\[
K_{II}/K_I\approx1.066\times10^{-4}.
\]

The SIF forensic audit should now be treated as closed, subject only to repository integration and preservation of the audit trail.

## Repository integration record

The audit was developed as a stacked sequence, but final integration is performed through the cumulative closure PR.

Target:

**sif-asymmetric-mesh-audit**

PR #45 is retargeted from its stacked parent to \`sif-asymmetric-mesh-audit\`. Because that target is an ancestor of the Step70 head, the cumulative PR contains the complete Step62--Step70 audit history without rewriting or dropping the intermediate commits.

The earlier stacked PRs #34--#44 remain useful provenance for the individual audit stages, including PR #41's failed ICT experiment. That negative result documents why the parameter-free SGS route was adopted without post-hoc preconditioner tuning.

After PR #45 is merged, perform one final static/documentation regression review on \`sif-asymmetric-mesh-audit\`. Only then decide separately whether any audited production changes should be proposed toward \`main\`.

