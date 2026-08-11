"""Transactional publication helpers for filesystem outputs."""

import os
from pathlib import Path
import tempfile
from typing import Iterable
from uuid import uuid4


def _lexists(path: Path) -> bool:
    """Return whether a path entry exists, including a dangling symlink."""
    return os.path.lexists(path)


def ensure_outputs_available(paths: Iterable[Path], *, force: bool) -> None:
    """Reject existing output entries unless replacement was requested."""
    if force:
        return
    existing = [path for path in paths if _lexists(path)]
    if existing:
        rendered = ", ".join(str(path) for path in existing)
        raise FileExistsError(
            f"Output already exists: {rendered}; pass force=True or --force to replace it"
        )


def sibling_temporary_path(target: Path, *, suffix: str) -> Path:
    """Create and return a temporary path in the target's directory."""
    descriptor, name = tempfile.mkstemp(
        prefix=f".{target.name}.",
        suffix=suffix,
        dir=target.parent,
    )
    os.close(descriptor)
    return Path(name)


def _backup_path(target: Path) -> Path:
    """Return an unused sibling path suitable for a hard-link backup."""
    while True:
        candidate = target.with_name(f".{target.name}.{uuid4().hex}.backup")
        if not _lexists(candidate):
            return candidate


def publish_outputs(
    replacements: Iterable[tuple[Path, Path]],
    *,
    force: bool,
) -> None:
    """Atomically publish completed sibling files, rolling back a failed pair.

    Each individual publication is atomic. For multiple related outputs, all
    temporary files must already be complete before publication begins, and a
    failure during publication restores the previous set.
    """
    pairs = tuple(replacements)
    ensure_outputs_available((target for _temporary, target in pairs), force=force)

    if not force:
        published: list[Path] = []
        try:
            for temporary, target in pairs:
                os.link(temporary, target)
                published.append(target)
        except BaseException:
            for target in reversed(published):
                target.unlink(missing_ok=True)
            raise
        else:
            for temporary, _target in pairs:
                temporary.unlink(missing_ok=True)
        return

    backups: dict[Path, Path | None] = {}
    published: list[Path] = []
    try:
        for _temporary, target in pairs:
            if _lexists(target):
                backup = _backup_path(target)
                os.link(target, backup, follow_symlinks=False)
                backups[target] = backup
            else:
                backups[target] = None

        for temporary, target in pairs:
            os.replace(temporary, target)
            published.append(target)
    except BaseException:
        for target in reversed(published):
            backup = backups.get(target)
            if backup is None:
                target.unlink(missing_ok=True)
            else:
                os.replace(backup, target)
        raise
    finally:
        for backup in backups.values():
            if backup is not None:
                backup.unlink(missing_ok=True)


def cleanup_paths(paths: Iterable[Path]) -> None:
    """Best-effort removal of unpublished temporary files."""
    for path in paths:
        path.unlink(missing_ok=True)
