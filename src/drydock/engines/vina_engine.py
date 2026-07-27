"""AutoDock Vina engine, covering the vina, vinardo and ad4 scoring functions.

Uses Vina's Python API rather than shelling out to the binary. That avoids
parsing terminal output, and more importantly it allows the receptor and its
grid maps to be set up **once per worker process** and reused for every ligand.

Why reuse is safe
-----------------

Map computation costs ~1.5 s against a typical receptor, against ~1-4 s to dock
one ligand. Recomputing per ligand would therefore add something like 50% to the
total runtime of a screen.

The reason it can be shared is that ``Vina.dock()`` re-seeds from the object's
seed on every call, so a ligand's result does not depend on what the object
docked before it. Verified directly: a ligand docked first on a fresh object and
the same ligand docked after three others on a reused object give identical
affinities to four decimal places. That is what makes one Vina object per worker
compatible with reproducible results -- without it, output would depend on how
work happened to be distributed across workers.
"""

from __future__ import annotations

import contextlib
import os
import sys
from collections.abc import Iterator
from pathlib import Path

from drydock.engines.base import DockConfig, DockOutcome, EngineError, parse_vina_poses


@contextlib.contextmanager
def _silenced_stdout() -> Iterator[None]:
    """Silence writes to file descriptor 1, including from C++.

    ``contextlib.redirect_stdout`` is not sufficient here. It rebinds Python's
    ``sys.stdout`` object, but Vina is a C++ extension that writes to the
    descriptor directly, so its output sails past any Python-level redirection.

    Left unhandled this is not merely untidy. Twelve worker processes writing to
    a descriptor they inherited from the parent share one file offset, and their
    interleaved writes punch holes that the kernel fills with NUL bytes: a single
    benchmark produced a 32 MB log of which 3.2 million bytes were NUL. Replacing
    the descriptor itself is the only redirection the extension cannot bypass.
    """
    sys.stdout.flush()
    saved = os.dup(1)
    devnull = os.open(os.devnull, os.O_WRONLY)
    try:
        os.dup2(devnull, 1)
        yield
    finally:
        os.dup2(saved, 1)
        os.close(devnull)
        os.close(saved)


class VinaEngine:
    """Docks ligands with AutoDock Vina, holding one prepared Vina instance."""

    def __init__(self, config: DockConfig) -> None:
        self.config = config
        self.name = config.engine
        self._vina = None

    def _ensure_ready(self):
        """Build the Vina object and its maps, once.

        Deferred rather than done in ``__init__`` so the engine can be pickled
        and sent to a worker: the underlying Vina object is a C++ handle and
        cannot cross a process boundary.
        """
        if self._vina is not None:
            return self._vina

        from vina import Vina

        config = self.config
        vina = Vina(
            sf_name=config.scoring_function,
            cpu=config.cpu,
            seed=config.seed,
            verbosity=0,
        )

        try:
            # Receptor loading and map computation are also chatty at the C++
            # level, and happen once per worker.
            with _silenced_stdout():
                if config.scoring_function == "ad4":
                    # AD4 scoring runs against pre-computed AutoGrid maps rather
                    # than maps Vina derives itself, because the AD4 force field
                    # -- and AutoDock4Zn's zinc terms in particular -- are defined
                    # by the grid parameter file, not by anything Vina can infer.
                    #
                    # Checked here rather than left to Vina: its own message is
                    # "Cannot find affinity maps with <path>", repeated once per
                    # ligand, which describes the symptom and not the fix.
                    from drydock.core.zinc import MAP_PREFIX, maps_status

                    ok, detail = maps_status(config.maps_dir)
                    if not ok:
                        raise EngineError(detail, setup=True)
                    vina.load_maps(str(Path(config.maps_dir) / MAP_PREFIX))
                else:
                    vina.set_receptor(config.receptor)
                    vina.compute_vina_maps(
                        center=list(config.box.center),
                        box_size=list(config.box.size),
                    )
        except EngineError:
            raise
        except Exception as exc:  # noqa: BLE001
            raise EngineError(_explain_setup_failure(exc, config), setup=True) from exc

        self._vina = vina
        return vina

    def dock(self, ligand_path: str, ligand_id: str, want_poses: bool = False) -> DockOutcome:
        """Dock one ligand and return its scored modes."""
        vina = self._ensure_ready()
        config = self.config

        try:
            vina.set_ligand_from_file(str(ligand_path))
        except Exception as exc:  # noqa: BLE001
            raise EngineError(f"could not load ligand: {_brief(exc)}", ligand_id) from exc

        try:
            with _silenced_stdout():
                vina.dock(
                    exhaustiveness=config.exhaustiveness,
                    n_poses=config.n_modes,
                )
                pdbqt = vina.poses(n_poses=config.n_modes)
        except Exception as exc:  # noqa: BLE001
            raise EngineError(f"docking failed: {_brief(exc)}", ligand_id) from exc

        modes = parse_vina_poses(pdbqt)
        if not modes:
            raise EngineError("engine returned no poses", ligand_id)

        return DockOutcome(modes=modes, poses_pdbqt=pdbqt if want_poses else None)


def _brief(exc: Exception) -> str:
    """First line of an exception, trimmed. Vina's errors run to many lines."""
    text = str(exc).strip().splitlines()
    return (text[0] if text else exc.__class__.__name__)[:200]


def _explain_setup_failure(exc: Exception, config: DockConfig) -> str:
    """Turn a receptor-loading failure into something actionable.

    Vina's most common rejection is an atom type in the wrong case -- ``ZN``
    rather than ``Zn``. The message it produces does say "atom types are
    case-sensitive", but buried in a C++ overload-resolution error that reads as
    a bug in Drydock rather than a fixable problem with the receptor.
    """
    text = str(exc)
    if "not a valid AutoDock type" in text:
        offending = ""
        for token in text.split():
            if token.strip(".,") and "Atom type" in text:
                break
        with contextlib.suppress(IndexError):
            offending = text.split("Atom type", 1)[1].split()[0]
        return (
            f"receptor {config.receptor} has an invalid AutoDock atom type"
            + (f" ({offending})" if offending else "")
            + ". AutoDock types are case-sensitive: metals use element casing "
            "(Zn, Mg, Mn, Fe, Ca, Cu), not upper case. "
            "Run 'drydock check-receptor' for details."
        )
    return f"could not set up receptor: {_brief(exc)}"


def build_engine(config: DockConfig) -> VinaEngine:
    """Create the engine for a configuration.

    Classic AutoDock4 is a separate executable with a different workflow and is
    handled elsewhere; everything Vina-driven comes through here.
    """
    if config.engine == "autodock4":
        raise EngineError(
            "the autodock4 engine is not wired up yet; use 'ad4' for AutoDock4 "
            "scoring at Vina speed, which is the right choice for library screening"
        )
    return VinaEngine(config)
