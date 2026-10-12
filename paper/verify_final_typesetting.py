"""Read-only safeguards for the typesetting-only change based on PR #80.

Run: python paper/verify_final_typesetting.py
Requires pypdf; no numerical calculation, plotting, or local MAT archive.
"""
from pathlib import Path
import argparse
import hashlib
import json
import re
import subprocess
from pypdf import PdfReader

PAPER = Path(__file__).resolve().parent
ROOT = PAPER.parent
BASE = '799f2e190fa90dd41b7fee178763d66ffcaada29'


def git(*args):
    return subprocess.check_output(['git', *args], cwd=ROOT)


def digest(data):
    return hashlib.sha256(data).hexdigest()


def normalized(text):
    return re.sub(r'\s+', '', text)


def without_comments(text):
    return re.sub(r'(?<!\\)%[^\n]*', '', text)


def environments(text, name):
    return re.findall(r'\\begin\{' + name + r'\}.*?\\end\{' + name + r'\}', text, re.S)


def strip_figures(text):
    return re.sub(r'\\begin\{figure\}.*?\\end\{figure\}', '', text, flags=re.S)


def caption_free(text):
    # Balanced TeX arguments, including nested mathematical subscripts.
    while (m := re.search(r'\\caption\{', text)):
        end, depth = m.end(), 1
        while depth:
            depth += (text[end] == '{') - (text[end] == '}')
            end += 1
        text = text[:m.start()] + text[end:]
    return text


def check(pdf_path=None):
    old = git('show', f'{BASE}:paper/main.tex').decode()
    new = (PAPER / 'main.tex').read_text()
    results = {}

    def require(name, condition, details=True):
        assert condition, name
        results[name] = details

    # Compare all display equations, including unnumbered displays and aligned systems.
    math = lambda s: [normalized(x) for x in re.findall(
        r'\\begin\{(?:equation|align)\}.*?\\end\{(?:equation|align)\}|(?<!\\)\\\[.*?(?<!\\)\\\]', s, re.S)]
    require('displayed_equations_unchanged', math(old) == math(new), len(math(old)))
    table = lambda s: [normalized(re.sub(r'\\begin\{table\}\[[^\]]+\]',
        r'\\begin{table}', x)) for x in environments(s, 'table')]
    require('complete_tables_unchanged_except_placement', table(old) == table(new), len(table(old)))
    audited = lambda s: normalized(s[s.index('% BEGIN AUDITED DATA'):s.index('% END AUDITED DATA')])
    require('audited_numerical_macros_unchanged', audited(old) == audited(new))
    cite = lambda s: re.findall(r'\\cite\w*\{([^}]+)\}', s)
    require('citation_keys_and_order_unchanged', cite(old) == cite(new), cite(new))
    labels = lambda s: re.findall(r'\\label\{([^}]+)\}', s)
    require('labels_and_order_unchanged', labels(old) == labels(new), len(labels(new)))
    paths = lambda s: re.findall(r'\\includegraphics\[[^\]]*\]\{([^}]+)\}', s)
    require('figure_paths_and_order_unchanged', paths(old) == paths(new), paths(new))
    layout = lambda s: re.findall(r'\\begin\{(?:minipage|subfigure)\}(?:\[[^\]]*\])?\{([^}]+)\}', s)
    require('all_panel_arrangements_and_sizes_unchanged', layout(old) == layout(new), len(layout(new)))
    titles = lambda s: re.findall(r'\\(?:sub)?section\{([^}]+)\}', s)
    require('section_titles_and_order_unchanged', titles(old) == titles(new), len(titles(new)))
    # Every numerical literal outside figure captions is retained, in the same order.
    # Layout lengths are excluded explicitly, not scientific quantities.
    def numerical_content(s):
        s = caption_free(without_comments(strip_figures(s)))
        s = re.sub(r'\\Needspace\{[^}]+\}', '', s)
        s = re.sub(r'\\setlength\{\\emergencystretch\}\{[^}]+\}', '', s)
        return re.findall(r'(?<![A-Za-z])[-+]?\d+(?:\.\d+)?', s)
    require('noncaption_numerical_literals_and_order_unchanged', numerical_content(old) == numerical_content(new), len(numerical_content(old)))
    def prose(s):
        s = strip_figures(without_comments(s[s.index('\\begin{document}'):]))
        s = s.replace('$(\\bm\\sigma,\\bm\\varepsilon,\\bm u)$, the integral is evaluated over\n'
                      'the material domain $\\Omega_k$ defined in Section~\\ref{sec:problem}.\n'
                      'The gradient of the radial weight is nonzero only in the interaction\n'
                      'annulus. Here $\\dd\\Omega$ is the area element,',
                      '$(\\bm\\sigma,\\bm\\varepsilon,\\bm u)$, let $\\Omega_k$ denote the integration\n'
                      'domain around the current tip. Here $\\dd\\Omega$ is the area element,')
        s = re.sub(r'\\(?:Needspace|bibliographystyle|bibliography)\{[^}]+\}', '', s)
        s = re.sub(r'\\setlength\{\\emergencystretch\}\{[^}]+\}', '', s)
        s = re.sub(r'\\(?:begin|end)\{samepage\}', '', s)
        s = re.sub(r'\\(?:FloatBarrier|clearpage|begingroup|endgroup|par)\b', '', s)
        s = re.sub(r'\\begin\{table\}\[[^\]]+\]', r'\\begin{table}', s)
        return normalized(s)
    require('all_nonfigure_prose_unchanged_except_authorized_Omega_correction', prose(old) == prose(new))
    author = lambda s: normalized(re.search(r'\\author\{.*?\\date\{\}', s, re.S)[0])
    require('author_affiliation_and_email_unchanged', author(old) == author(new))
    # Only main.tex/main.pdf and scratch-ignore policy may differ from the baseline.
    changed = git('diff', '--name-only', BASE, '--').decode().splitlines()
    allowed = {'paper/main.tex', 'paper/main.pdf', '.gitignore',
               'paper/verify_final_typesetting.py', 'paper/FINAL_TYPESETTING_REVIEW.md',
               'paper/review_final_typesetting.py'}
    require('no_unapproved_tracked_changes', all(p in allowed or p.startswith('paper/reviews/final_typesetting/') for p in changed), changed)
    protected = git('ls-tree', '-r', '--name-only', BASE, 'paper/data', 'paper/figures').decode().splitlines()
    hashes = {}
    for path in protected:
        raw = (ROOT / path).read_bytes()
        require('protected:' + path, raw == git('show', f'{BASE}:{path}'))
        hashes[path] = digest(raw)
    require('portable_data_and_vector_evidence_byte_identical', True, hashes)
    approved = PAPER / 'reviews/issue78/after'
    for path in paths(new):
        require('approved_candidate:' + path, (PAPER / path).read_bytes() == (approved / path).read_bytes())
    geometry = {p['asset']: p for p in json.loads((PAPER / 'reviews/issue78/geometry_report.json').read_text())['panels']}
    for path in paths(new):
        if '/m1_mesh_sensitivity/' in path:
            text = PdfReader(PAPER / path).pages[0].extract_text()
            if Path(path).name in {'trajectory.pdf', 'mode_mixity.pdf'}:
                require('active_M0_M1:' + path, 'M0' in text and 'M1' in text and 'reference' not in text.lower())
            else:
                # MATLAB outlined axis glyphs have no extractable text. Bind their
                # verified label record to the active candidate's exact PDF hash.
                g = geometry[path]
                require('active_M0_M1:' + path, g['ylabelReplacement']['label'] == 'M1 minus M0'
                        and g['output_sha256'] == digest((PAPER / path).read_bytes()))
    require('Omega_material_domain_not_redefined', 'let $\\Omega_k$ denote the integration' not in new
            and 'the material domain $\\Omega_k$ defined in Section' in new)
    require('appendix_before_bibliography', new.index('\\appendix') < new.index('\\bibliography{references}'))
    require('S1_designation_and_anchor_preserved', '\\renewcommand{\\thefigure}{S\\arabic{figure}}' in new
            and '\\renewcommand{\\theHfigure}{S\\arabic{figure}}' in new)
    require('caption_policy_local_to_three_figures', new.count('justification=raggedright,singlelinecheck=false') == 3)
    require('no_clearpage_or_graphics_trim', '\\clearpage' not in new and not re.search(r'\\includegraphics\[[^\]]*trim', new))
    # The accepted/unsolved distinction is protected by unchanged prose outside captions.
    require('P24_unsolved_distinction_retained', 'has no accepted physical' in new
            and 'No failed $P_{24}$ solve is plotted.' in new
            and 'qualified but the elastic problem at' in new)
    pdf_path = pdf_path or PAPER / 'main.pdf'
    pdf = PdfReader(pdf_path)
    require('sole_author_metadata', pdf.metadata.get('/Author') == 'Danylo Selivanov')
    require('title_metadata', pdf.metadata.get('/Title') == 'Incremental LEFM Prediction of a Crack Trajectory from a Circular Hole')
    text = [p.extract_text() for p in pdf.pages]
    require('resolved_PDF_text', all('??' not in t for t in text))
    require('S1_in_compiled_PDF', any('Figure S1:' in t for t in text))
    require('references_conclude_PDF', '[19]' in text[-1] and 'Supplementary numerical diagnostics' not in text[-1])
    links = [str(a.get_object()['/A']['/D']) for p in pdf.pages for a in p.get('/Annots', [])
             if a.get_object().get('/A', {}).get('/S') == '/GoTo']
    require('all_internal_links_resolve', all(d in pdf.named_destinations for d in links), len(links))
    results['compiled_pdf_sha256'] = digest(pdf_path.read_bytes())
    results['pages'] = len(pdf.pages)
    results['baseline_commit'] = BASE
    return results


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path)
    parser.add_argument('--pdf', type=Path, help='Verify a newly reproduced PDF; defaults to the committed paper/main.pdf.')
    parser.add_argument('--self-test', action='store_true', help='Reject five in-memory scientific regressions without changing any file.')
    args = parser.parse_args()
    results = check(args.pdf)
    if args.self_test:
        from unittest.mock import patch
        original = Path.read_text
        text = (PAPER / 'main.tex').read_text()
        mutations = [
            ('equation', "\\frac{E'}{2K_{\\mathrm{aux}}}", "\\frac{E'}{3K_{\\mathrm{aux}}}", 'displayed_equations_unchanged'),
            ('macro', '\\LSLengthMM}{84.15}', '\\LSLengthMM}{84.16}', 'audited_numerical_macros_unchanged'),
            ('table', '0.054025 &', '0.054026 &', 'complete_tables_unchanged_except_placement'),
            ('citation', '\\cite{Bowie1956}', '\\cite{unknown}', 'citation_keys_and_order_unchanged'),
            ('asset', 'figures/m1_mesh_sensitivity/mode_mixity.pdf', 'figures/historical/mode_mixity.pdf', 'figure_paths_and_order_unchanged'),
        ]
        rejected = {}
        for name, old, new, expected in mutations:
            assert old in text
            changed = text.replace(old, new, 1)
            def read(path, *a, **kw):
                return changed if path == PAPER / 'main.tex' else original(path, *a, **kw)
            with patch.object(Path, 'read_text', read):
                try:
                    check()
                except AssertionError as exc:
                    assert str(exc) == expected, (name, exc)
                    rejected[name] = expected
                else:
                    raise AssertionError(f'Mutation {name} was accepted')
        results['in_memory_regression_rejections'] = rejected
    if args.output:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(json.dumps(results, indent=2) + '\n')
    print(f"PASS: {len(results)} checks; {results['pages']} pages; no physical calculations.")
