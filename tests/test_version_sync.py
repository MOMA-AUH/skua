from pathlib import Path
import re
import tomllib


ROOT = Path(__file__).resolve().parents[1]


def _extract_version(path: Path, pattern: str) -> str:
    text = path.read_text(encoding="utf-8")
    match = re.search(pattern, text, flags=re.MULTILINE)
    assert match is not None, f"Could not find version in {path}"
    return match.group(1)


def test_version_is_synced_across_code_and_packaging() -> None:
    package_version = _extract_version(
        ROOT / "src" / "skua" / "_version.py",
        r'^__version__\s*=\s*"([^"]+)"$',
    )
    pyproject = tomllib.loads((ROOT / "pyproject.toml").read_text(encoding="utf-8"))
    conda_recipe = (ROOT / "conda-recipe" / "meta.yaml").read_text(encoding="utf-8")

    assert pyproject["project"]["dynamic"] == ["version"]
    assert "version" not in pyproject["project"]
    assert pyproject["tool"]["setuptools"]["dynamic"]["version"] == {
        "attr": "skua._version.__version__"
    }
    conda_version = _extract_version(
        ROOT / "conda-recipe" / "meta.yaml",
        r'^\{\%\s*set\s+version\s*=\s*"([^"]+)"\s*\%\}$',
    )
    assert conda_version == package_version
