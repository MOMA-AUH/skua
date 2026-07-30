#!/usr/bin/env python3
"""Bump skua version across code and packaging metadata.

Usage:
    python tools/set_version.py 0.1.2
"""

from pathlib import Path
import re
import sys


ROOT = Path(__file__).resolve().parents[1]
VERSION_PATTERN = re.compile(r"^\d+\.\d+\.\d+(?:[a-zA-Z0-9_.+-]*)?$")


VERSION_FILE = ROOT / "src" / "skua" / "_version.py"
VERSION_PATTERN_IN_FILE = re.compile(
    r'^(__version__\s*=\s*")([^"]+)(")$', flags=re.MULTILINE
)
CONDA_RECIPE = ROOT / "conda-recipe" / "meta.yaml"
CONDA_VERSION_PATTERN = re.compile(
    r'^(\{\%\s*set\s+version\s*=\s*")([^"]+)("\s*\%\})$', flags=re.MULTILINE
)


def _replace_version_once(text: str, pattern: re.Pattern[str], new_version: str, path: Path) -> tuple[str, str]:
    matches = list(pattern.finditer(text))
    if len(matches) != 1:
        raise RuntimeError(f"Expected exactly one version match in {path}, found {len(matches)}")

    old_version = matches[0].group(2)
    updated = pattern.sub(rf"\g<1>{new_version}\g<3>", text, count=1)
    return updated, old_version


def main(argv: list[str]) -> int:
    if len(argv) != 2:
        print("Usage: python tools/bump_version.py <new-version>")
        return 2

    new_version = argv[1].strip()
    if not VERSION_PATTERN.match(new_version):
        print(f"Invalid version: {new_version}")
        return 2

    targets = (
        (VERSION_FILE, VERSION_PATTERN_IN_FILE),
        (CONDA_RECIPE, CONDA_VERSION_PATTERN),
    )
    old_versions: set[str] = set()
    for path, pattern in targets:
        text = path.read_text(encoding="utf-8")
        updated, old_version = _replace_version_once(text, pattern, new_version, path)
        path.write_text(updated, encoding="utf-8")
        old_versions.add(old_version)

    old_versions_display = ", ".join(sorted(old_versions))
    print(f"Updated version(s) [{old_versions_display}] -> {new_version}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv))
