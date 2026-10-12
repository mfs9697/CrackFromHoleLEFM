# Archive of Applied Mechanics — submission preparation

**Journal:** Archive of Applied Mechanics (Springer Nature)  
**Article:** *Incremental LEFM Prediction of a Crack Trajectory from a Circular Hole*  
**Author and corresponding author:** Danylo Selivanov, `fract@inmech.kyiv.ua`  
**Source baseline:** `paper-crack-trajectory` at `1b8dd5f43dcc5e43205d34cc526523832fc852e7`

This is the **submission-preparation branch**. No manuscript has been sent to
the journal, and neither `paper/main.tex` nor the verified `paper/main.pdf`
should be modified while preparing journal-specific files.

## Official requirements checked on 2026-10-12

- Journal instructions: https://link.springer.com/journal/419/submission-guidelines
- Springer Nature LaTeX author support:
  https://www.springernature.com/gp/authors/campaigns/latex-author-support
- Publisher template ZIP: use the journal-author download link on the latter page.
- The journal accepts LaTeX, asks for editable sources, original figures and
  a compiled manuscript PDF. A graphical/video abstract is **not listed as a
  mandatory component**.
- Abstract: 150–250 words. Keywords: 4–6. Numeric references with full DOI
  links if available.
- Required: publication declarations, data-availability statement and
  author-approved disclosure of competing interests; author-contribution and competing
  interest information are also entered in the submission interface.
- Appendix figures retained *in the main paper* must continue ordinary
  figure numbering. The authoring candidate therefore turns prior Fig. S1
  into Fig. 10. Alternatively we can move it to a separately submitted
  `ESM_1.pdf`, restore a separately numbered supplementary figure, and cite
  Online Resource 1. Do not do both.
- The official journal page **conflicts internally**: one paragraph recommends
  the template option `[iicol]`, another says single-column. This authoring
  candidate starts as single-column. Compare both after compilation and seek
  clarification if the submission interface imposes a layout.
- There is no need for a graphical or video abstract based on these instructions.

## Draft candidate and scientific freeze

`manuscript.tex` is an isolated Springer-class candidate. It preserves:
- exact abstract text and all audited scientific macro definitions;
- original narrative, equations and numerical tables, in sequence;
- all 31 exact `\\includegraphics` paths;
- 19 citation-key entries referenced through the unchanged
  `paper/references.bib`.

Permitted differences in the submission candidate:
- Springer front matter and the class-selected numbered mathematical-physics bibliography style;
- six indexing keywords (COD omitted as keyword, **not** from the paper);
- figure numbering in an inline appendix (S1 becomes Fig. 10);
- original scientific sources remain untouched.

For development, `\\graphicspath{{../../}}` resolves the existing vector
figures and `\\bibliography{../../references}` uses the audited bibliography
when TeXstudio compiles from this directory. **This layout is not upload-ready**:
the final ZIP must include all source and style files and all referenced
graphic assets, preferably in a flat upload directory to avoid path issues.

## First local typesetting test

1. Download the **official** Springer Nature LaTeX template package. Do not
   use a third-party modified `sn-jnl.cls`. The class/style files are
   intentionally not copied into this Git repository.
2. In a local development copy of this directory, place `sn-jnl.cls` and
   `sn-mathphys-num.bst` from the official December 2024 (v3.1) package beside `manuscript.tex`. The style file is under `sn-article-template/bst/`; the `sn-mathphys-num` class option chooses the corresponding bibliography style automatically, so **do not add a separate bibliography-style command**.
3. Open `paper/submissions/aam/manuscript.tex` in TeXstudio and set the
   working directory to `paper/submissions/aam/`.
4. Compile using pdfLaTeX, BibTeX, pdfLaTeX, pdfLaTeX. Review the log and
   all pages; do **not** treat the earlier 23-page `paper/main.pdf` as the
   Springer proof.
5. Review wide figures, subcaptions, numbering, bibliography links, long tables,
   all mathematics, orphan paragraphs, and metadata. Do not adjust scientific
   values, digit precision, curves or plot labels to repair layout.
6. Adapt the layout/asset packaging in a **separate reviewed commit**, and run
   the scientific-freeze check after each edit.

If the new TeX class conflicts with `subcaption`, `placeins`, or another
package, capture the log and fix the compatibility issue rather than
silently altering scientific material.

## Author decisions still required — do not infer

- **Funding confirmed:** regular institutional funding of the S.P. Timoshenko Institute of Mechanics, National Academy of Sciences of Ukraine; no specific grant. The candidate now includes this Funding declaration.
- Competing interests: full financial/non-financial disclosure for Danylo.
- Data/code availability: precise repository URL and version/commit; accurately
  describe retained evidence and unavailable raw histories.
- Author's ORCID if any, and confirmation of the official English affiliation.
- AI assistance: make an accurate distinction between copyediting and
  substantive generation/software assistance, in accordance with Springer
  policy; author verifies wording.
- Choose whether Figure S1 remains as main-paper Fig. 10 or becomes
  independently submitted online supplementary material.
- Optional acknowledgments: only with permission from named persons.
- Cover letter and novelty/overlap review; reviewer suggestions optional.

`cover_letter_draft.md` is editable and **not yet approved for submission**.

**Gate:** Do not upload to Springer Nature's submission portal before the
compiled Springer PDF, all declarations, package integrity, author approval,
and final source-to-PDF consistency checks are complete.
