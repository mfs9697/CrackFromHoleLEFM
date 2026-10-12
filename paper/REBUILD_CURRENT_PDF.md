# Corrective PDF/source verification — pending final PDF

This corrective branch begins at \`paper-scientific-batch2\` following PR #74
(merge \`96c311acc356a056c2de29bfb682fedc66c0691f\`).
It addresses the two reviewer findings on PR #81, without modifying
\`paper/main.tex\`, any numerical result, any figure, or \`references.bib\`.

## Why a new PDF is necessary

The tracked \`paper/main.pdf\` was produced before the final academic-language
and PCG revisions. It is **not a valid proof of the current manuscript**.
The existing old PDF and old visual-review JSON remain historical until a
rebuild is completed. Do not label either one as final.

## Rebuild from current source (Windows / TeXstudio)

1. Preserve any modified local \`paper/main.pdf\` before changing worktrees.
2. Check out this PR's branch \`paper-corrective-pdf-source-sync\`. Do not
   discard changes in another worktree.
3. Compile \`paper/main.tex\` with the current TeXstudio LaTeX/BibTeX toolchain,
   retaining all 31 active \`paper/figures\` PDFs. The project currently uses
   the ordinary \`article\` class, not yet the Archive of Applied Mechanics
   template.
4. From the repository root, run:

   \`\`\`text
   python paper/verify_pdf_source_match.py --record
   python paper/verify_pdf_source_match.py
   python paper/verify_final_typesetting.py --pdf paper/main.pdf --self-test
   \`\`\`

   \`--record\` only writes \`paper/reviews/final_typesetting/pdf_source_manifest.json\`
   if all source/PDF textual correspondence checks already pass. The second
   command fails if the PDF, TeX, bibliography, or any of the 31 active figure
   files subsequently change. The third runs all existing content safeguards,
   plus the new fail-closed correspondence check.

5. Commit the **freshly compiled** \`paper/main.pdf\` and recorded JSON manifest
   to this corrective branch. Rebuild the affected-page visual comparison and
   review manifest only if those proofs are to be advertised as current;
   alternatively, clearly label them as the archived pre-language-edits proof.
6. Visually inspect the new PDF, especially Sections 3.3 and 4.6,
   Figures 7–9 and the bibliography. Record the inspection before merge.

## Verified scope of the change

This PR must remain a **draft** until the compiled PDF and manifest have
been committed and the checks pass. GitHub PR #73 must remain unmerged until
this corrective PR has been approved and merged into
\`paper-scientific-batch2\`.

The guards deliberately require the following current-text distinctions:

- \`a_k = k Delta a\` is the accumulated crack-path length from \`P_0\`.
- The PDF contains the revised software-neutral PCG stagnation/restart wording.
- The abstract does not imply that coupled refinement isolates the increment.
- Figure 7 uses \`h_0\`/\`2h_0\`; Figure 8 uses \`M0\`/\`M1\`.
- No instance of "MATLAB" remains in the final compiled manuscript.
- The compiled PDF is linked to the exact TeX/bibliography/graphics bytes.

This is a **build and provenance** repair, not a numerical or theoretical
revision.
