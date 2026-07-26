"""Drydock command line interface.

Commands are thin: they parse arguments and hand off to :mod:`drydock.core`.
Everything reachable here is reachable from the API too, which is what lets the
GUI drive the same code rather than reimplementing it.

Heavy imports (RDKit, Vina) are deferred into the command bodies so that
``drydock --help`` and ``drydock status`` stay instant.
"""

from __future__ import annotations

import json
from pathlib import Path

import click

from drydock import __version__


@click.group(context_settings={"help_option_names": ["-h", "--help"]})
@click.version_option(__version__, prog_name="drydock")
def main() -> None:
    """Reproducible high-throughput virtual screening.

    Prepare a compound library, then dock it against a receptor you prepared
    yourself, and get a ranked CSV out.
    """


@main.command("prep-ligands")
@click.option(
    "-i", "--input", "library", required=True, type=click.Path(exists=True, path_type=Path),
    help="Compound library (.sdf/.smi/.mol2, optionally .gz).",
)
@click.option(
    "-o", "--out", required=True, type=click.Path(path_type=Path),
    help="Directory to write PDBQTs and the manifest into.",
)
@click.option("--ph", default=7.4, show_default=True, help="pH for protonation states.")
@click.option(
    "--optimize/--no-optimize", default=True, show_default=True,
    help="Embed in 3D and minimise with MMFF94s, rather than keeping input coordinates.",
)
@click.option(
    "--conformers", default=1, show_default=True,
    help="Starting geometries per compound. >1 is the only setting that samples ring space.",
)
@click.option(
    "--macrocycles/--rigid-macrocycles", default=True, show_default=True,
    help="Let Vina open and re-close macrocyclic rings during the search.",
)
@click.option("--tautomers", is_flag=True, help="Enumerate tautomers (multiplies the library).")
@click.option("--max-torsions", type=int, default=None, help="Reject floppier ligands than this.")
@click.option("--id-field", default=None, help="SDF property holding the compound identifier.")
@click.option("-j", "--workers", default=0, help="Processes to use. 0 means one per CPU.")
@click.option("--limit", type=int, default=None, help="Only read this far into the library.")
@click.option("--seed", default=0, show_default=True, help="Global seed for reproducibility.")
@click.option("--no-resume", is_flag=True, help="Re-prepare ligands already in the manifest.")
def prep_ligands(
    library: Path,
    out: Path,
    ph: float,
    optimize: bool,
    conformers: int,
    macrocycles: bool,
    tautomers: bool,
    max_torsions: int | None,
    id_field: str | None,
    workers: int,
    limit: int | None,
    seed: int,
    no_resume: bool,
) -> None:
    """Convert a compound library into docking-ready PDBQT files.

    Streams the input, so library size is bounded by disk rather than memory, and
    resumes by default: re-running the same command after an interruption picks
    up where it stopped.
    """
    from drydock.core.ligprep import PrepConfig
    from drydock.core.prep_runner import run_prep

    config = PrepConfig(
        ph=ph,
        optimize=optimize,
        n_conformers=conformers,
        skip_tautomers=not tautomers,
        macrocycles=macrocycles,
        seed=seed,
        max_torsions=max_torsions,
    )

    with click.progressbar(length=100, label="preparing", show_pos=False) as bar:
        state = {"last": 0}

        def on_progress(summary) -> None:
            done = summary.prepared + summary.failed + summary.skipped
            pct = int(100 * done / summary.total) if summary.total else 0
            bar.update(max(0, pct - state["last"]))
            state["last"] = pct

        summary = run_prep(
            library,
            out,
            config,
            n_workers=workers,
            id_field=id_field,
            limit=limit,
            resume=not no_resume,
            progress=on_progress,
        )

    click.echo(
        f"prepared {summary.prepared}, failed {summary.failed}, "
        f"skipped {summary.skipped} in {_format_duration(summary.elapsed_s)} "
        f"({summary.rate_per_s:.1f}/s)"
    )
    click.echo(f"manifest: {out}/manifest.csv")
    if summary.failed:
        click.echo(f"failures: {out}/failures.csv")


@main.command("survey")
@click.argument("library", type=click.Path(exists=True, path_type=Path))
@click.option("--id-field", default=None, help="SDF property holding the compound identifier.")
def survey_cmd(library: Path, id_field: str | None) -> None:
    """Summarise a library without preparing it.

    Reports records against distinct compounds. These differ whenever a source
    enumerates stereoisomers under one identifier, and the gap matters: records
    determine how long a screen takes, distinct compounds determine how many
    things you are actually testing.
    """
    from drydock.core.library import survey

    info = survey(library, id_field=id_field)
    click.echo(f"records:                 {info['records']}")
    click.echo(f"distinct compounds:      {info['distinct_compounds']}")
    click.echo(f"compounds with variants: {info['compounds_with_variants']}")
    click.echo(f"most variants:           {info['max_variants']}")
    if info["most_repeated"]:
        click.echo("most repeated identifiers:")
        for cid, n in info["most_repeated"]:
            click.echo(f"  {cid}: {n}")


@main.command()
@click.argument("run_dir", type=click.Path(path_type=Path))
@click.option("--json", "as_json", is_flag=True, help="Emit machine-readable JSON.")
def status(run_dir: Path, as_json: bool) -> None:
    """Report progress of a run.

    Reads the run directory only, so it is safe to call against a run in flight.
    """
    from drydock.core.rundir import RunDir

    run = RunDir(run_dir)
    if not run.exists():
        raise click.ClickException(f"no run directory at {run_dir}")

    # Prefer the cached summary, but fall back to the journal when it is absent
    # or unreadable, since the journal is the authority.
    current = run.read_status() or run.rebuild_status()

    if as_json:
        click.echo(json.dumps(current.to_dict(), indent=2))
        return

    click.echo(f"run:       {run.path}")
    click.echo(f"state:     {current.state}")
    if current.engine:
        click.echo(f"engine:    {current.engine}")
    click.echo(
        f"progress:  {current.done}/{current.total} "
        f"({current.fraction:.1%})  ok={current.completed} "
        f"failed={current.failed} skipped={current.skipped}"
    )
    if (rate := current.rate_per_s) is not None:
        click.echo(f"rate:      {rate:.2f} ligands/s")
    if (eta := current.eta_s) is not None:
        click.echo(f"eta:       {_format_duration(eta)}")


@main.group()
def dev() -> None:
    """Development helpers. Not part of the screening workflow."""


@dev.command("synthesize-run")
@click.argument("run_dir", type=click.Path(path_type=Path))
@click.option("-n", "--n-ligands", default=500, show_default=True, help="Ligands to fabricate.")
@click.option("--seed", default=0, show_default=True, help="Makes the output reproducible.")
@click.option("--failure-rate", default=0.02, show_default=True, help="Fraction marked failed.")
@click.option("--live", is_flag=True, help="Write slowly, so a watcher sees it progress.")
@click.option("--delay", default=0.05, show_default=True, help="Seconds between records if --live.")
@click.option("--incomplete", is_flag=True, help="Stop partway, leaving the run marked running.")
def synthesize_run(
    run_dir: Path,
    n_ligands: int,
    seed: int,
    failure_rate: float,
    live: bool,
    delay: float,
    incomplete: bool,
) -> None:
    """Fabricate a run directory with plausible results.

    Lets the GUI be built and exercised before any engine exists, and makes
    awkward states -- a killed run, a run full of failures, a run large enough to
    stress the results table -- reproducible on demand.
    """
    from drydock.core.synthetic import synthesize_run as _synthesize

    run = _synthesize(
        run_dir,
        n_ligands=n_ligands,
        seed=seed,
        failure_rate=failure_rate,
        live=live,
        delay_s=delay,
        complete=not incomplete,
    )
    current = run.read_status()
    click.echo(f"wrote {current.done} records to {run.path}")


def _format_duration(seconds: float) -> str:
    seconds = int(seconds)
    hours, rem = divmod(seconds, 3600)
    minutes, secs = divmod(rem, 60)
    if hours:
        return f"{hours}h {minutes}m"
    if minutes:
        return f"{minutes}m {secs}s"
    return f"{secs}s"


if __name__ == "__main__":
    main()
