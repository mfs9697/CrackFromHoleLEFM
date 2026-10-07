# Stage-I source provenance audit

Status: **open pending local-source identification**.

This audit is intentionally separate from the manuscript prose. No numerical
solver, mesh generator, SIF extractor, or manuscript text is changed here.

## What is established from Git/GitHub

The pilot branch contains one committed Stage-I MAT snapshot:

- path: `paper/data/accepted_stage1_source.mat`;
- size: **42401 bytes**;
- Git blob SHA: **`9c6904005b5880b9ed32174bbf185865f940f799`**;
- SHA-256 recorded in `paper/data/source_manifest.csv`:
  **`93bf96d515b8719258f8fbd5fdfeffa1d009e17c286f10fe4dc4e19c02e64a17`**.

The file first entered Git in commit
`21e4cc3198514176b05b3d3280b6e2a31b5708e6`
(`Add audited pilot manuscript on incremental crack trajectory`), whose
parent is the plotting snapshot
`cfdf0f110010d688a4e3c48f6d88a00fd17dc698`.

The current branch tree contains neither
`verification/crack_path/stage1_starting_state.mat` nor a separately
committed `accepted_R0.mat`. Therefore Git history cannot by itself prove
the path or byte identity of the local file from which the manuscript copy
was made.

The pilot evidence layer reports that the copy contains an `R0` whose
configuration, initiation summary, gates, mouth point, and local frame agree
with the clean-run configuration and frozen Stage-I fingerprints. That is a
strong **content-consistency check**, but it is not yet a complete
**source-provenance chain**.

## Claim currently requiring closure

`paper/EVIDENCE_AND_PROVENANCE.md` states that an exact investigator-saved
`accepted_R0.mat`, saved for the preceding profiling task, was copied
byte-for-byte to `paper/data/accepted_stage1_source.mat`.

The destination fingerprint is preserved, but the exact original local path
was not recorded in Git. Until that original path is independently located
and hash-matched, treat the claim as **provisional provenance**, not as a
closed recovery of the formerly missing standard
`stage1_starting_state.mat`.

In particular, the two statements are compatible:

1. the standard file `stage1_starting_state.mat` was missing during the
   earlier home/Drive audit; and
2. an accepted `R0` may have survived under a different filename in a
   profiling/output location.

The second statement still needs local filesystem confirmation.

## Local closure test

Run from the repository root on the computer where Codex prepared the pilot:

```powershell
python paper/audit_stage1_source_provenance.py
```

The script is read-only. It verifies the committed destination fingerprint
and scans repository-local MAT files for a byte-identical copy, with special
reporting of files named `accepted_R0.mat`.

If Codex used a source outside the repository tree, either provide its exact
path:

```powershell
python paper/audit_stage1_source_provenance.py --original-file "C:\path\to\accepted_R0.mat"
```

or widen the scan explicitly:

```powershell
python paper/audit_stage1_source_provenance.py --search-root "C:\Users\Mikhailo"
```

The wide home-directory scan can take longer because every MAT file of the
same size is hashed.

### PASS criterion

The provenance audit closes only when:

- the committed destination has the recorded size, SHA-256, and Git blob;
- at least one separately located local source is byte-identical to it; and
- the exact source path is recorded in the audit record.

Expected closed status from the helper:

```text
STATUS: PASS_LOCAL_ORIGINAL_IDENTIFIED
```

If the helper reports

```text
STATUS: DESTINATION_VERIFIED_SOURCE_NOT_IDENTIFIED
```

the Stage-I values remain usable as an internally consistent preserved
snapshot, but the exact origin of that snapshot is still unresolved.

## What not to do

Do not regenerate Stage I merely to satisfy this provenance check, do not
rename a reconstructed MAT file as the historical source, and do not change
the manuscript's scientific claims until the source-path audit is closed.
