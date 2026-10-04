# Step 66: solver-memory preflight for the C03 Level-1 mesh

## Purpose

Step66 answers one question before any Level-1 physical solve is considered:

> Is the present direct linear-algebra path likely to fit safely in memory, and if not, what solver qualification should come next?

Step66 is **symbolic only**.

It performs:

- no physical FEM solve;
- no physical displacement load;
- no stiffness-value assembly;
- no load-vector assembly;
- no `K\F`;
- no PCG/GMRES/MINRES;
- no physical COD or EDI.

## Mesh reconstruction

Step66 deterministically reconstructs Level 0 and Level 1 directly from the committed archived Step62 baseline using the same C03 calibration.

Level 0 must reproduce:

- 32,980 T3 triangles;
- 66,854 T6 nodes;
- max adjacent-size ratio 1.79678451 within tolerance.

Level 1 must preserve:

- scale 0.5;
- paired radius 6 mm;
- max adjacent-size ratio <=1.8;
- complete T3/T6 pairing;
- full structural pass.

The shared mesh builder gains a default-off `ReturnCandidate` option so Step66 can receive connectivity in memory without writing another mesh package. Existing callers are unchanged.

## Exact assembly-triplet cost

`stif_assem.m` currently allocates three double triplet arrays `rw`, `cl`, and `st` with exactly 144 entries per T6 element.

Therefore Step66 reports the exact storage of those three arrays:

`3 * 8 bytes * 144 * nElem`.

This is not the full assembly peak because MATLAB also needs the resulting sparse matrix and temporary workspace during duplicate summation.

## Structural K sparsity

For each T6 mesh, Step66 constructs only the six-node element adjacency graph.

A structurally connected pair of T6 nodes corresponds to a potentially dense 2x2 displacement block. If `A_node` is the symmetric node graph, the full scalar stiffness-pattern estimate is therefore

`nnz(K_pattern) = 4 * nnz(A_node)`.

No element stiffness values are evaluated.

## Symbolic direct-factor estimate

Step66 applies `symamd` to the node graph and `symbfact` to the permuted graph.

If the symbolic lower factor contains `nnzL_node` node blocks, the corresponding 2-DOF scalar lower-factor estimate is

`nnzL_scalar = 4*nnzL_node - nNode`.

The subtraction accounts for each diagonal 2x2 lower block having three rather than four scalar entries.

This is a **block-symbolic SPD Cholesky estimate**, not a measurement of MATLAB's actual `backslash` peak memory.

## Important production-solver observation

The current `stif_assem.m` enforces homogeneous Dirichlet conditions by:

1. zeroing constrained rows;
2. setting the constrained diagonal entries to one;
3. leaving the corresponding columns unchanged.

That row-only clamping makes the stored matrix nonsymmetric, even though the unconstrained elasticity operator with homogeneous essential conditions is SPD.

Consequently, the current `K\F` path is not guaranteed to use sparse Cholesky and can require more memory than the SPD symbolic estimate.

Step66 does **not** change this behavior. It records it as a reason to qualify any memory-saving solver formulation on Level 0 before using it for Level 1.

## Memory quantities

Step66 reports for each level:

- T3 element count;
- T6 node count;
- displacement DOF count;
- exact triplet-array storage;
- node-pattern nnz;
- scalar K-pattern nnz;
- approximate MATLAB sparse-K storage;
- symbolic factor nnz and fill ratio;
- approximate SPD factor storage.

It also prints three planning proxies:

- `assemblyPeakProxyGiB = 2*triplets + Kpattern`;
- `spdPlanningLowGiB = triplets + Kpattern + 2*factor`;
- `spdPlanningHighGiB = 2*triplets + Kpattern + 4*factor`.

These are deliberately labeled **heuristics**. They provide a conservative planning scale; they are not promises about MATLAB peak usage.

## Machine memory

On Windows, Step66 attempts to read currently available physical memory through MATLAB's `memory` function.

An explicit override can be supplied:

    R66 = main_step66_solver_memory_preflight( ...
        'AvailableMemoryGiB', 15.0);

Use the override only when you know the amount of memory actually available to MATLAB.

## Decision logic

Step66 returns a risk class and a recommended next qualification.

Possible recommendations are intentionally limited to Level-0 validation steps, for example:

- qualify a symmetric-Dirichlet direct formulation on Level 0;
- qualify an iterative solver on Level 0;
- qualify both memory-saving assembly and an iterative solver if assembly itself is at risk.

Step66 never authorizes or performs the Level-1 physical solve.

## Local run

On branch:

`audit/step66-solver-memory-preflight`

run:

    addpath(genpath(pwd));

    R66 = main_step66_solver_memory_preflight();

    disp(R66.MemoryTable);
    disp(R66.ScalingTable);
    disp(R66.MachineMemory);
    disp(R66.Decision);

Return the complete console output. The key quantities are the Level-1 triplet memory, K-pattern storage, symbolic factor storage, high planning envelope, current available memory, and the recommended next qualification.
