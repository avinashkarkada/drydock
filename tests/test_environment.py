"""Guards on environment isolation.

Drydock's reproducibility claim rests on the lockfile describing the whole
environment. Anything importable that the lockfile does not name breaks that
claim, and does so silently -- the failure surfaces later, on someone else's
machine, as a version skew nobody can reproduce.

These tests are cheap and they fail loudly, which is the point.
"""

from __future__ import annotations

import shutil
import site
import subprocess
import sys
from pathlib import Path


def test_user_site_packages_is_disabled():
    """~/.local must not be on sys.path.

    Python enables the per-user site directory by default, which merges whatever
    the user pip-installed globally into this environment. Real example from
    development: a stray `zarr` in ~/.local pulled in `donfig`, which imported
    `yaml`, which was not installed -- and collapsed the entire test run before
    a single test executed.
    """
    assert not site.ENABLE_USER_SITE, (
        "user site-packages is enabled; set PYTHONNOUSERSITE=1 in [activation.env]"
    )


def test_no_foreign_paths_on_sys_path():
    """Every sys.path entry must be the pinned env or this checkout.

    Drydock is installed editable, so the repo root and ``src/`` are legitimately
    on the path. Anything outside both is an environment leak.
    """
    prefix = Path(sys.prefix).resolve()
    repo_root = Path(__file__).resolve().parent.parent

    foreign = []
    for entry in sys.path:
        if not entry or entry.endswith(".zip"):
            continue
        resolved = Path(entry).resolve()
        if resolved.is_relative_to(prefix) or resolved.is_relative_to(repo_root):
            continue
        foreign.append(entry)

    assert not foreign, f"sys.path escapes the pinned environment: {foreign}"


def test_interpreter_comes_from_the_pinned_environment():
    assert ".pixi" in sys.executable, f"unexpected interpreter: {sys.executable}"


def test_engine_binaries_resolve_inside_the_environment():
    """A binary found on the host PATH is not the one the lockfile pinned."""
    for binary in ("vina", "autogrid4", "autodock4"):
        resolved = shutil.which(binary)
        assert resolved, f"{binary} not found on PATH"
        assert ".pixi" in resolved, f"{binary} resolves outside the environment: {resolved}"


def test_autogrid_is_the_pinned_version():
    """autogrid ships separately from autodock and is easy to omit entirely."""
    out = subprocess.run(["autogrid4", "--version"], capture_output=True, text=True)
    assert "4.2.6" in (out.stdout + out.stderr)


def test_core_scientific_stack_imports():
    """One import each, so a broken install fails here rather than mid-screen."""
    import meeko  # noqa: F401
    import molscrub  # noqa: F401
    import rdkit  # noqa: F401
    from vina import Vina

    # Both scoring paths must construct: `ad4` in particular depends on parts of
    # the build that a Python-only smoke test would not touch.
    Vina(sf_name="vina", verbosity=0)
    Vina(sf_name="ad4", verbosity=0)
