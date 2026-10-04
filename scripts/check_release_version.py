"""Fail publication if source, recipe, or optional release tag disagree."""

import argparse
from pathlib import Path
import re
import sys


def extract_version(path: Path, pattern: str) -> str:
    match = re.search(pattern, path.read_text(encoding="utf-8"), re.MULTILINE)
    if match is None:
        raise ValueError(f"Cannot read version from {path}")
    return match.group(1)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--tag", help="Release tag, including the v prefix")
    parser.add_argument("--source-root", type=Path, default=Path(__file__).resolve().parents[1])
    args = parser.parse_args()
    try:
        version = extract_version(
            args.source_root / "src/skua/_version.py", r'^__version__\s*=\s*"([^"]+)"$',
        )
        recipe_version = extract_version(
            args.source_root / "conda-recipe/meta.yaml",
            r'^\{%\s*set\s+version\s*=\s*"([^"]+)"\s*%\}$',
        )
        if recipe_version != version:
            raise ValueError(f"Conda version {recipe_version} does not match package {version}")
        if args.tag is not None and args.tag != f"v{version}":
            raise ValueError(f"Release tag {args.tag} does not match package v{version}")
    except (OSError, ValueError) as exc:
        print(str(exc), file=sys.stderr)
        return 1
    print(version)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
