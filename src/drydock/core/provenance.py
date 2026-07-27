"""Recording what produced a result.

A screen produces a CSV of numbers. Six months later, or in review, the question
is what produced them: which receptor, which box, which engine, which version
of Vina, which seed. Reconstructing that from memory is unreliable, and from a
shell history worse.

So every run writes ``provenance.json``, and it records enough to answer that
question without reference to anything outside the run directory.

Checksums rather than paths
---------------------------

Inputs are identified by SHA-256 as well as by path. A path says where a file was
on one machine on one day; a checksum says whether the file someone has now is
the file that was actually screened. Receptors in particular get re-prepared,
re-minimised and re-saved under the same name, and a run that cannot prove which
version it used cannot really be reproduced.

Large files are hashed by streaming, so a 174 MB library costs a second and no
memory worth mentioning.
"""

from __future__ import annotations

import hashlib
import platform
import subprocess
import sys
import time
from pathlib import Path
from typing import Any

from drydock import __version__

# Packages whose version can change a result. Recorded so a later discrepancy can
# be traced to a dependency bump rather than argued about.
TRACKED_PACKAGES: tuple[str, ...] = (
    "vina",
    "rdkit",
    "meeko",
    "molscrub",
    "numpy",
)

# Read in chunks rather than whole: a compound library can be hundreds of MB, and
# nothing here needs it resident.
_HASH_CHUNK = 1 << 20


def file_sha256(path: str | Path) -> str | None:
    """Hash a file by streaming it. Returns None if it does not exist."""
    path = Path(path)
    if not path.is_file():
        return None

    digest = hashlib.sha256()
    with open(path, "rb") as fh:
        while chunk := fh.read(_HASH_CHUNK):
            digest.update(chunk)
    return digest.hexdigest()


def directory_digest(path: str | Path, pattern: str = "*.pdbqt") -> dict[str, Any]:
    """Summarise a directory of inputs.

    Hashing every ligand in a 47,000-file directory would take longer than some
    screens, so this records the file count, total size and a hash over the
    sorted (name, size) pairs. That catches files being added, removed, renamed
    or resized, which covers how a prepared library actually changes, and stays
    instant.
    """
    path = Path(path)
    if not path.is_dir():
        return {"path": str(path), "exists": False}

    entries = sorted((p.name, p.stat().st_size) for p in path.glob(pattern))
    digest = hashlib.sha256()
    for name, size in entries:
        digest.update(f"{name}:{size}\n".encode())

    return {
        "path": str(path.resolve()),
        "n_files": len(entries),
        "total_bytes": sum(size for _, size in entries),
        "listing_sha256": digest.hexdigest(),
    }


def package_versions() -> dict[str, str]:
    """Resolved versions of the packages that can change a result."""
    from importlib.metadata import PackageNotFoundError, version

    found: dict[str, str] = {}
    for name in TRACKED_PACKAGES:
        try:
            found[name] = version(name)
        except PackageNotFoundError:
            found[name] = "not installed"
    return found


def engine_versions() -> dict[str, str]:
    """Versions of the external binaries, which are not Python packages.

    autogrid4 and autodock4 are separate executables whose versions are invisible
    to Python's package metadata, yet they determine the numbers on the AD4 path.
    """
    versions: dict[str, str] = {}
    for binary, args in (("autogrid4", ["--version"]), ("autodock4", ["--version"])):
        try:
            result = subprocess.run(
                [binary, *args], capture_output=True, text=True, timeout=10
            )
            text = (result.stdout + result.stderr).strip().splitlines()
            versions[binary] = text[0][:80] if text else "unknown"
        except (OSError, subprocess.SubprocessError):
            versions[binary] = "not available"
    return versions


def environment() -> dict[str, Any]:
    """The machine and interpreter a run happened on."""
    return {
        "python": sys.version.split()[0],
        "executable": sys.executable,
        "platform": platform.platform(),
        "machine": platform.machine(),
        "hostname": platform.node(),
    }


def build(
    *,
    receptor: str | Path,
    ligand_dir: str | Path,
    config: Any,
    manifest: str | Path | None = None,
    extra: dict[str, Any] | None = None,
) -> dict[str, Any]:
    """Assemble the provenance record for a run.

    Args:
        receptor: Receptor PDBQT. Hashed in full.
        ligand_dir: Prepared ligand directory. Summarised, not hashed per file.
        config: The run's :class:`~drydock.engines.base.DockConfig`.
        manifest: Ligand manifest, if one exists.
        extra: Anything else worth recording.

    Returns:
        A JSON-serialisable dictionary.
    """
    receptor = Path(receptor)

    record: dict[str, Any] = {
        "drydock_version": __version__,
        "created_at": time.time(),
        "created_at_iso": time.strftime("%Y-%m-%dT%H:%M:%S%z"),
        "environment": environment(),
        "packages": package_versions(),
        "engines": engine_versions(),
        "receptor": {
            "path": str(receptor.resolve()),
            "sha256": file_sha256(receptor),
            "bytes": receptor.stat().st_size if receptor.is_file() else None,
        },
        "ligands": directory_digest(ligand_dir),
        "config": config.to_dict() if hasattr(config, "to_dict") else dict(config),
    }

    if manifest:
        manifest = Path(manifest)
        record["manifest"] = {
            "path": str(manifest.resolve()),
            "sha256": file_sha256(manifest),
        }

    # Stated explicitly rather than left implicit in cpu=1, because it is the
    # claim a reader most needs to evaluate.
    record["reproducibility"] = {
        "global_seed": getattr(config, "seed", None),
        "threads_per_job": getattr(config, "cpu", None),
        "note": (
            "AutoDock Vina is not reproducible across threads even with a fixed "
            "seed, because the parallel search makes thread scheduling affect "
            "which minima are found. Runs are therefore single-threaded per "
            "ligand, parallel across ligands."
        ),
    }

    if extra:
        record.update(extra)
    return record
