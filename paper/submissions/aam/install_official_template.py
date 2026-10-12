r"""Install the exact official Springer Nature template files from an author ZIP.

Usage:
    python paper/submissions/aam/install_official_template.py "C:\path\to\official-template.zip"

Keeps publisher-owned files out of the Git repository history. It never
silently overwrites a different class or BibTeX style file.
"""
from __future__ import annotations

import argparse
from pathlib import Path
from zipfile import ZipFile, BadZipFile

VERSION = "%Version 3.1 December 2024"
FILES = {
    "sn-jnl.cls": "sn-article-template/sn-jnl.cls",
    "sn-mathphys-num.bst": "sn-article-template/bst/sn-mathphys-num.bst",
}


def install(zip_path: Path, destination: Path) -> None:
    with ZipFile(zip_path) as package:
        sample = package.read("sn-article-template/sn-article.tex")
        if not sample.decode("latin-1").startswith(VERSION):
            raise ValueError("Expected the official December 2024, version 3.1 template")
        for name, member in FILES.items():
            contents = package.read(member)
            if not contents:
                raise ValueError(f"Empty official template member: {member}")
            target = destination / name
            if target.exists():
                if target.read_bytes() == contents:
                    print(f"UNCHANGED {target}")
                    continue
                raise FileExistsError(
                    f"Different file already exists: {target}. "
                    "Review it manually before replacing the publisher file."
                )
            target.write_bytes(contents)
            print(f"INSTALLED {target}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("template_zip", type=Path, help="Official Springer Nature template ZIP")
    args = parser.parse_args()
    try:
        install(args.template_zip, Path(__file__).resolve().parent)
    except (OSError, ValueError, KeyError, BadZipFile) as exc:
        parser.exit(1, f"ERROR: {exc}\n")
