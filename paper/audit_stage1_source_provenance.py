#!/usr/bin/env python3
"""Audit provenance of the recovered Stage-I R0 snapshot.

This script is read-only.  It does not run MATLAB, regenerate Stage I, or
modify any numerical archive.  It checks the committed manuscript copy and
looks for a byte-identical local source file.

Examples
--------
python paper/audit_stage1_source_provenance.py
python paper/audit_stage1_source_provenance.py --original-file "C:\\path\\to\\accepted_R0.mat"
python paper/audit_stage1_source_provenance.py --search-root "C:\\Users\\Mikhailo"
"""

from __future__ import annotations

import argparse
import hashlib
import json
import os
import subprocess
from pathlib import Path

EXPECTED_SIZE = 42401
EXPECTED_SHA256 = "93bf96d515b8719258f8fbd5fdfeffa1d009e17c286f10fe4dc4e19c02e64a17"
EXPECTED_GIT_BLOB = "9c6904005b5880b9ed32174bbf185865f940f799"


def sha256_file(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        for chunk in iter(lambda: f.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def git_blob_sha(repo: Path, path: Path) -> str | None:
    try:
        return subprocess.check_output(
            ["git", "hash-object", str(path)],
            cwd=repo,
            text=True,
            stderr=subprocess.DEVNULL,
        ).strip()
    except Exception:
        return None


def iter_mat_files(root: Path):
    skip_dirs = {".git", "__pycache__", "build"}
    for dirpath, dirnames, filenames in os.walk(root):
        dirnames[:] = [d for d in dirnames if d not in skip_dirs]
        for name in filenames:
            if name.lower().endswith(".mat"):
                yield Path(dirpath) / name


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument(
        "--search-root",
        type=Path,
        default=None,
        help="Root to scan for a byte-identical MAT file; defaults to repository root.",
    )
    ap.add_argument(
        "--original-file",
        type=Path,
        default=None,
        help="Known original file to compare directly.",
    )
    ap.add_argument(
        "--json-out",
        type=Path,
        default=None,
        help="Optional path for a machine-readable audit report.",
    )
    args = ap.parse_args()

    paper = Path(__file__).resolve().parent
    repo = paper.parent
    dest = paper / "data" / "accepted_stage1_source.mat"
    search_root = (args.search_root or repo).resolve()

    if not dest.is_file():
        raise SystemExit(f"FAIL: committed destination is missing: {dest}")

    dest_size = dest.stat().st_size
    dest_sha = sha256_file(dest)
    dest_blob = git_blob_sha(repo, dest)

    destination_ok = (
        dest_size == EXPECTED_SIZE
        and dest_sha == EXPECTED_SHA256
        and (dest_blob is None or dest_blob == EXPECTED_GIT_BLOB)
    )

    exact_name = []
    identical = []

    if args.original_file is not None:
        candidates = [args.original_file.resolve()]
    else:
        candidates = list(iter_mat_files(search_root))

    dest_resolved = dest.resolve()
    for p in candidates:
        try:
            if not p.is_file() or p.resolve() == dest_resolved:
                continue
            if p.name.lower() == "accepted_r0.mat":
                exact_name.append(str(p))
            if p.stat().st_size != dest_size:
                continue
            if sha256_file(p) == dest_sha:
                identical.append(str(p))
        except (OSError, PermissionError):
            continue

    if not destination_ok:
        status = "FAIL_DESTINATION_FINGERPRINT"
    elif identical:
        status = "PASS_LOCAL_ORIGINAL_IDENTIFIED"
    else:
        status = "DESTINATION_VERIFIED_SOURCE_NOT_IDENTIFIED"

    report = {
        "status": status,
        "repository_root": str(repo),
        "search_root": str(search_root),
        "destination": str(dest),
        "destination_size": dest_size,
        "destination_sha256": dest_sha,
        "destination_git_blob": dest_blob,
        "expected_size": EXPECTED_SIZE,
        "expected_sha256": EXPECTED_SHA256,
        "expected_git_blob": EXPECTED_GIT_BLOB,
        "destination_fingerprint_ok": destination_ok,
        "exact_name_candidates": exact_name,
        "byte_identical_local_copies": identical,
        "original_file_argument": str(args.original_file.resolve()) if args.original_file else None,
    }

    print("=" * 68)
    print("STAGE-I R0 PROVENANCE AUDIT")
    print("=" * 68)
    print(f"Destination : {dest}")
    print(f"Size        : {dest_size} bytes")
    print(f"SHA-256     : {dest_sha}")
    print(f"Git blob    : {dest_blob or '(unavailable)'}")
    print(f"Fingerprint : {'PASS' if destination_ok else 'FAIL'}")
    print()
    print(f"Search root : {search_root}")
    print(f"Exact-name candidates outside paper copy : {len(exact_name)}")
    for p in exact_name:
        print(f"  NAME  {p}")
    print(f"Byte-identical local copies              : {len(identical)}")
    for p in identical:
        print(f"  MATCH {p}")
    print()
    print(f"STATUS: {status}")

    if status == "PASS_LOCAL_ORIGINAL_IDENTIFIED":
        print("The committed manuscript copy is byte-identical to at least one")
        print("separately located local source. Preserve the reported path(s).")
    elif status == "DESTINATION_VERIFIED_SOURCE_NOT_IDENTIFIED":
        print("The committed copy is internally verified, but the original local")
        print("source path is not established by this scan. Do not promote the")
        print("copy to canonical recovered provenance without locating the source")
        print("or recording an independent archival chain.")
    else:
        print("The committed Stage-I copy does not match its recorded fingerprint.")

    if args.json_out is not None:
        out = args.json_out.resolve()
        out.parent.mkdir(parents=True, exist_ok=True)
        out.write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")
        print(f"JSON report : {out}")

    return 0 if destination_ok else 2


if __name__ == "__main__":
    raise SystemExit(main())
