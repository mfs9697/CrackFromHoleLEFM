"""Copy the IAM-UA8-F02 run evidence (exact scripts that ran, logs incl. failed attempts, verification JSON) into
paper/ukrainian_review/audit/iam_f02_evidence/ with local paths sanitized; UTF-16 logs to UTF-8; MATLAB startup path warnings removed.
Writes SANITIZATION.json (source and published SHA-256 per file)."""
import hashlib, json, os
SCR = r'<scratch>'
RD = os.path.join(SCR, 'f02')
DST = r'<repo>\paper\ukrainian_review\audit\iam_f02_evidence'
PAT = [(SCR, '<scratch>'), (r'<temp>', '<temp>'), (r'<generator_cfdf0f1>', '<generator_cfdf0f1>'),
       (r'<repo>', '<repo>'), (r'<investigator_codex_work>', '<investigator_codex_work>'),
       (r'<matlabroot>', '<matlabroot>'), (r'<drive>', '<drive>'), (r'<home>', '<home>')]


def variants(p):
    v = {p, p.replace('\\', '/'), p.replace('\\', '\\\\'), p.lower(), p.replace('\\', '/').lower(), p.replace('\\', '\\\\').lower()}
    if p.startswith('C:'):
        v |= {'/c' + p[2:].replace('\\', '/')}
    return sorted(v, key=len, reverse=True)


def sanitize(t):
    n = {}
    for p, r in PAT:
        for v in variants(p):
            k = t.count(v)
            if k:
                t = t.replace(v, r); n[r] = n.get(r, 0) + k
    return t, n


files = ['f02_regenerate.m', 'f02_regen.log', 'f02_verify_export.m', 'f02_verify_attempt1_FAILED.log', 'f02_verify.log', 'F02_VERIFY.json',
         'f02_verify_p5_t3.m', 'f02_verify_p5.log', 'F02_VERIFY_P5_T3.json', 'f02_design_diff.m', 'f02_design_diff_attempt1_FAILED.log',
         'f02_design_diff.log', 'F02_VERIFY_P5_DESIGN_DIFF.json', 'f02_probe.py', 'f02_repro_run.m', 'f02_repro.log', 'repro/IAM_F02_REGENERATION.json',
         'publish_evidence.py']
os.makedirs(DST, exist_ok=True); rec = []
for f in files:
    b = open(os.path.join(RD, *f.split('/')), 'rb').read()
    t = b.decode('utf-16') if b[:2] in (b'\xff\xfe', b'\xfe\xff') else b.decode('utf-8')
    t = t.lstrip('\ufeff').replace('\r\n', '\n'); t, n = sanitize(t)
    if f.endswith('.log'):
        t = '\n'.join(l for l in t.split('\n') if not l.startswith('[Warning: Name is nonexistent'))
    name = f.replace('/', '_'); out = os.path.join(DST, name); assert not os.path.exists(out), out
    ob = t.encode('utf-8'); open(out, 'wb').write(ob)
    rec.append({'file': name, 'source': f, 'sourceSHA256': hashlib.sha256(b).hexdigest(), 'publishedSHA256': hashlib.sha256(ob).hexdigest(), 'replacements': n})
json.dump({'schema': 'IAM-UA8-F02.Sanitization.v1', 'rule': 'local paths replaced by placeholders; UTF-16 logs converted to UTF-8; MATLAB startup path warnings removed',
           'files': rec}, open(os.path.join(DST, 'SANITIZATION.json'), 'w', encoding='utf-8', newline='\n'), indent=1)
print('published', len(rec))
