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
