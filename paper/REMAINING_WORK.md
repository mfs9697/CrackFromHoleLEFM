# Remaining work before journal submission

1. **Establish trajectory convergence.** The current path uses only 4-mm
   increments. Compare 4/2/1-mm sequences at common physical crack lengths,
   controlling mesh and extraction resolution independently. Report changes
   in tip coordinates, absolute directions, SIFs, and the sign-change bracket.
   Historical straight 8-mm mesh convergence does not answer this question.

2. **Check the late-tip extraction and mesh.** At least P17, P21, P22, and
   P23 need physically matched mesh refinement and EDI-domain sensitivity.
   The COD differences are real: the largest turning difference is about
   0.01024 degrees at P23; near P21 a rear-window linear COD ratio is 26.3%
   above EDI despite a small absolute difference. All eight fits agree on
   the sign-change bracket, but that is not a precision bound on its location.

3. **Test the boundary interpretation.** Move the right boundary farther
   away in a controlled comparison, or compare matched fixed crack
   geometries, while stating how initiation/frame/geometry are held or
   recomputed. The present sequence changes curvature and boundary distance
   together and cannot establish a uniquely boundary-caused reversal.

4. **Assemble verified literature.** No complete bibliography entries were
   found in the inspected repository. Fill the MTS, hole/notch trajectory,
   incremental LEFM, interaction-integral, COD, and quadrature TODOs using
   verified sources. Establish what is new before claiming novelty.

5. **Complete provenance of the initial and terminal states.** An exact
   investigator-saved R0 was recovered and copied without reconstruction to
   `data/accepted_stage1_source.mat`; its configuration and frame match the
   clean run. The standard historical Stage-I path remains absent, and the
   originating P1 physical checkpoint is not in the clean compact archive.
   Preserve that original P1 evidence or its full provenance. Retain the
   original P24 residual history/failure log: stagnation is reported in the
   task description but cannot be independently quantified from the current
   compact archive. Do not retrospectively lower the acceptance gate.

6. **Keep the physical scope explicit.** This is direction prediction for
   imposed finite advances in quasi-static LEFM. No propagation threshold,
   plastic-zone assessment, dynamic/branching model, arrest criterion, or
   boundary-intersection prediction is established. The Stage-I record also
   has a secondary tensile maximum about 1.34% below the primary; justify
   the single-active-crack restriction without assuming no subsequent
   initiation can occur after stress redistribution.

7. **Finalize article presentation only after the science is strengthened.**
   Review equation-level conventions with the authors, add verified
   citations, and check figure labels/legend size at the chosen publication
   width. Keep the 84.15-mm interpolation clearly distinct from solved
   states. Add author details and data/code availability statements only
   when supplied; do not optimize for a publisher template yet.
