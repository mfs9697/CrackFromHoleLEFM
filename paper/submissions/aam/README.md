# Archive of Applied Mechanics — submission preparation

**Journal:** Archive of Applied Mechanics (Springer Nature)  
**Article:** *Incremental LEFM Crack-Path Prediction from an Off-Center Circular Hole: Mode-Mixity Reversal and Numerical Sensitivity*  
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

1. Use the author's attached **official December 2024, version 3.1**
   template ZIP. Do not substitute a third-party class file.
2. From the repository root, run the verified local installer with the
   **actual path** where you saved the ZIP:

   ```powershell
   python paper/submissions/aam/install_official_template.py "C:\path\to\springer-template.zip"
   ```

   It extracts only `sn-jnl.cls` and `sn-mathphys-num.bst` beside
   `manuscript.tex`. Files are ignored by Git and never overwritten
   silently. The class option `sn-mathphys-num` selects the bibliography
   style automatically; do not add a manual `\\bibliographystyle` command.
3. Before typesetting, run the scientific preservation check:

   ```powershell
   python paper/submissions/aam/check_scientific_freeze.py
   ```

4. Open `paper/submissions/aam/manuscript.tex` in TeXstudio, set the
   working directory to `paper/submissions/aam/`, and compile with
   pdfLaTeX, BibTeX, pdfLaTeX, pdfLaTeX. Review the log and
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

## First full-manuscript compilation issue (2026-10-12)

The official December 2024 `sn-jnl.cls` redefines `table` with an internal
`center` and `threeparttable` environment. Its wrapper fails when the
`float` package's forced `[H]` placement is used. In a minimal test with the
same class and packages, a table using `[H]` failed at `\\end{table}`; the
same table with `[htbp]` compiled successfully.

Two original tables used `[H]`: the centered-hole one-step directional map
and the incremental-path pairwise comparison. **Only in the journal candidate**
these two placements are changed to `[htbp]`. The freeze checker reverses
those exact two journal-only changes before byte-comparing all scientific
content to `paper/main.tex`. No entries or numerical values were changed.

The `Extra }` / `threeparttable` errors from the first table cause LaTeX to
remain inside an environment, producing misleading downstream `Not in outer
par mode` errors and a terminal `center` mismatch. Fixing the first
environments should eliminate the cascade; assess any remaining errors in a
**fresh** log, not the previous failed log.

After a clean pdfLaTeX run, process citations with BibTeX using **the same
build basename** (`manuscript`) and then run pdfLaTeX twice. Undefined
references and `No file manuscript.bbl` from the *failed* LaTeX run do not
by themselves imply a bibliography defect. If the post-BibTeX run still
reports errors, send the first error and the relevant `.blg` lines.



## Typesetting audit after the first full PDF (2026-10-12)

The author-supplied 29-page `manuscript.pdf` was visually checked. A follow-up
coordinate audit identified two table overhangs into the right text margin
(on PDF pp. 23 and 25), in addition to the overly narrow Table 1 caption
(on PDF p. 13). No clipping at the physical paper edge was observed.

The active AAM candidate now contains **journal-only, reversible** changes:
- **Table 1:** a 0.64\\linewidth `tabular*` gives the long caption a
  readable width. The numerical table cells are unchanged.
- **Table 5:** the units under the mode-mixity ratio and turning-density
  headers are wrapped onto second lines, preserving both formulas and units.
- **Table 6:** the angular-difference unit is wrapped onto a second line.
- **PDF metadata:** explicit document title and Danylo Selivanov author.

A standalone threeparttable/LaTeX layout prototype was compiled to check
Table 1 caption behavior and the proposed Table 5/6 header widths.
**The exact Springer 3.1 manuscript has NOT been recompiled after these
changes; final line breaks and margin conformity still require author
verification under the official class.** The scientific-freeze checker
reverses only the listed formatting changes.

### Author's next local verification

1. Pull branch `paper-aam-template-compatibility`, then from
   `paper/submissions/aam/` run `python check_scientific_freeze.py`.
2. In TeXstudio, rebuild `manuscript.tex` using
   pdfLaTeX → BibTeX → pdfLaTeX → pdfLaTeX.
3. Inspect the **newly rebuilt** PDF: Table 1 on/near p. 13; Tables 5
   and 6 on/near pp. 23 and 25 (pagination can change). All headings,
   numerical cells, and captions must be legible and inside the text margins.
4. Inspect the first genuine error/warning from `manuscript.log` and check
   for `Overfull \\hbox` or unresolved citations/references; report any
   overfull lines with their associated line numbers.
5. Verify PDF document properties report a nonempty title and author.
   Upload the new PDF and LaTeX log for final visual sign-off.

Do not merge the draft PR or assemble the submission ZIP until the actual
publisher-template build and full proof audit pass.

## Author decisions still required — do not infer

- **Funding confirmed:** regular institutional funding of the S.P. Timoshenko Institute of Mechanics, National Academy of Sciences of Ukraine; no specific grant. The candidate now includes this Funding declaration.
- **Competing interests confirmed:** Danylo states that he has no financial or non-financial competing interests. The AAM candidate now includes this declaration.
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
