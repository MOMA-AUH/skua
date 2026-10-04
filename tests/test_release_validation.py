from pathlib import Path
import subprocess
import sys

from skua import __version__


ROOT = Path(__file__).resolve().parents[1]


def test_release_version_check_rejects_a_mismatched_tag() -> None:
    result = subprocess.run(
        [sys.executable, str(ROOT / "scripts/check_release_version.py"), "--tag", "v99.0.0"],
        capture_output=True, text=True,
    )
    assert result.returncode == 1
    assert "Release tag v99.0.0 does not match" in result.stderr


def test_release_version_check_accepts_the_matching_tag() -> None:
    result = subprocess.run(
        [sys.executable, str(ROOT / "scripts/check_release_version.py"),
         "--tag", f"v{__version__}"],
        capture_output=True, text=True,
    )
    assert result.returncode == 0, result.stderr
    assert result.stdout.strip() == __version__


def test_release_version_check_rejects_a_mismatched_recipe(tmp_path: Path) -> None:
    for relative in ("src/skua/_version.py", "conda-recipe/meta.yaml"):
        destination = tmp_path / relative
        destination.parent.mkdir(parents=True)
        destination.write_bytes((ROOT / relative).read_bytes())
    (tmp_path / "conda-recipe/meta.yaml").write_text('{% set version = "99.0.0" %}\n')
    result = subprocess.run(
        [sys.executable, str(ROOT / "scripts/check_release_version.py"),
         "--source-root", str(tmp_path)],
        capture_output=True, text=True,
    )
    assert result.returncode == 1
    assert "Conda version 99.0.0 does not match" in result.stderr
