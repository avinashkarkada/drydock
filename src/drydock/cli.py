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


@main.command("screen")
@click.option(
    "-r", "--receptor", required=True, type=click.Path(exists=True, path_type=Path),
    help="Prepared receptor PDBQT.",
)
@click.option(
    "-l", "--ligands", required=True, type=click.Path(exists=True, path_type=Path),
    help="Prepared ligand directory (from prep-ligands).",
)
@click.option(
    "-o", "--run", "run_dir", required=True, type=click.Path(path_type=Path),
    help="Run directory to create or resume.",
)
@click.option("--residues", default=None, help="Residues defining the box, e.g. 187,188,401.")
@click.option("--chain", default=None, help="Restrict residue selection to one chain.")
@click.option("--pad", default=5.0, show_default=True, help="Box padding in Angstroms.")
@click.option("--center", default=None, help="Explicit box centre as x,y,z.")
@click.option("--size", default=None, help="Explicit box size as x,y,z.")
@click.option(
    "-e", "--engine", default="vina", show_default=True,
    type=click.Choice(["vina", "vinardo", "ad4", "autodock4"]),
    help="Scoring function. Use ad4 for metalloproteins (AutoDock4Zn).",
)
@click.option("--exhaustiveness", default=8, show_default=True, help="Search effort per ligand.")
@click.option("--modes", default=9, show_default=True, help="Binding modes to report.")
@click.option("--seed", default=0, show_default=True, help="Global seed for reproducibility.")
@click.option("-j", "--workers", default=0, help="Parallel docking jobs. 0 means one per CPU.")
@click.option("--limit", type=int, default=None, help="Only screen this many ligands.")
@click.option("--maps", type=click.Path(path_type=Path), default=None, help="AutoGrid maps (ad4).")
@click.option("--no-resume", is_flag=True, help="Re-dock ligands already journalled.")
@click.option("--flat", is_flag=True, help="Do not group stereoisomers in results.csv.")
def screen_cmd(
    receptor: Path,
    ligands: Path,
    run_dir: Path,
    residues: str | None,
    chain: str | None,
    pad: float,
    center: str | None,
    size: str | None,
    engine: str,
    exhaustiveness: int,
    modes: int,
    seed: int,
    workers: int,
    limit: int | None,
    maps: Path | None,
    no_resume: bool,
    flat: bool,
) -> None:
    """Dock a prepared library against a prepared receptor.

    Writes every finished ligand to the run directory immediately, so an
    interrupted screen resumes without losing work, and the GUI can watch a run
    it does not own.
    """
    from drydock.core.box import Box
    from drydock.core.prep_runner import PrepDir
    from drydock.core.receptor import box_from_residues, inspect
    from drydock.core.results import write_results
    from drydock.core.rundir import RunDir
    from drydock.core.screen import run_screen, write_config
    from drydock.engines.base import DockConfig

    report = inspect(receptor)
    for problem in report.problems:
        click.secho(f"receptor problem: {problem}", fg="red")
    if report.problems:
        click.secho(
            "continuing anyway -- results will be affected. "
            "Run 'drydock check-receptor' for detail.",
            fg="yellow",
        )

    if center and size:
        box = Box(
            tuple(float(v) for v in center.split(",")),
            tuple(float(v) for v in size.split(",")),
        )
    elif residues:
        numbers = [int(v) for v in residues.replace(" ", "").split(",") if v]
        try:
            box, _ = box_from_residues(receptor, numbers, padding=pad, chain=chain)
        except ValueError as exc:
            raise click.ClickException(str(exc)) from exc
    else:
        raise click.ClickException("give either --residues or both --center and --size")

    click.echo(f"box: {box}")
    for warning in box.warnings():
        click.secho(f"warning: {warning}", fg="yellow")

    prep = PrepDir(ligands)
    pdbqt_dir = prep.pdbqt_dir if prep.pdbqt_dir.is_dir() else ligands

    config = DockConfig(
        receptor=str(receptor),
        box=box,
        engine=engine,
        exhaustiveness=exhaustiveness,
        n_modes=modes,
        seed=seed,
        cpu=1,
        maps_dir=str(maps) if maps else None,
    )

    run = RunDir(run_dir).create()
    write_config(run, config, str(pdbqt_dir))

    with click.progressbar(length=100, label=f"docking ({engine})") as bar:
        state = {"last": 0}

        def on_progress(status) -> None:
            pct = int(100 * status.fraction)
            bar.update(max(0, pct - state["last"]))
            state["last"] = pct

        status = run_screen(
            run_dir, pdbqt_dir, config,
            n_workers=workers, limit=limit, resume=not no_resume, progress=on_progress,
        )

    click.echo(
        f"docked {status.completed}, failed {status.failed} "
        f"in {_format_duration((status.finished_at or 0) - (status.started_at or 0))}"
    )

    results, all_modes, n = write_results(
        run_dir, prep.manifest_file, group_stereoisomers=not flat
    )
    click.echo(f"results:   {results}  ({n} rows)")
    click.echo(f"all modes: {all_modes}")


@main.command("benchmark")
@click.option(
    "-r", "--receptor", required=True, type=click.Path(exists=True, path_type=Path),
    help="Prepared receptor PDBQT.",
)
@click.option(
    "-l", "--ligands", required=True, type=click.Path(exists=True, path_type=Path),
    help="Prepared ligand directory.",
)
@click.option("--residues", default=None, help="Residues defining the box.")
@click.option("--chain", default=None, help="Restrict residue selection to one chain.")
@click.option("--pad", default=5.0, show_default=True, help="Box padding in Angstroms.")
@click.option("--center", default=None, help="Explicit box centre as x,y,z.")
@click.option("--size", default=None, help="Explicit box size as x,y,z.")
@click.option(
    "-e", "--engine", default="vina", show_default=True,
    type=click.Choice(["vina", "vinardo", "ad4"]),
)
@click.option(
    "--exhaustiveness", "exhaustiveness_values", multiple=True, type=int,
    help="Search effort to test. Repeat to compare several.",
)
@click.option("-n", "--sample", default=25, show_default=True, help="Ligands to time.")
@click.option("-j", "--workers", default=0, help="Parallel jobs. 0 means one per CPU.")
@click.option(
    "--library-size", default=None, type=int,
    help="Library size to project onto. Defaults to the prepared directory's size.",
)
@click.option("--seed", default=0, show_default=True, help="Seed for sampling and docking.")
def benchmark_cmd(
    receptor: Path,
    ligands: Path,
    residues: str | None,
    chain: str | None,
    pad: float,
    center: str | None,
    size: str | None,
    engine: str,
    exhaustiveness_values: tuple[int, ...],
    sample: int,
    workers: int,
    library_size: int | None,
    seed: int,
) -> None:
    """Measure throughput before committing to a long run.

    Docks a representative sample and projects the full run. Worth doing every
    time: cost per ligand varies by more than an order of magnitude with box
    volume and ligand flexibility, so an estimate from another target is worth
    very little.
    """
    import os

    from drydock.core.benchmark import benchmark, format_comparison
    from drydock.core.box import Box
    from drydock.core.prep_runner import PrepDir
    from drydock.core.receptor import box_from_residues
    from drydock.core.screen import iter_ligands
    from drydock.engines.base import DockConfig

    if center and size:
        box = Box(
            tuple(float(v) for v in center.split(",")),
            tuple(float(v) for v in size.split(",")),
        )
    elif residues:
        numbers = [int(v) for v in residues.replace(" ", "").split(",") if v]
        try:
            box, _ = box_from_residues(receptor, numbers, padding=pad, chain=chain)
        except ValueError as exc:
            raise click.ClickException(str(exc)) from exc
    else:
        raise click.ClickException("give either --residues or both --center and --size")

    prep = PrepDir(ligands)
    pdbqt_dir = prep.pdbqt_dir if prep.pdbqt_dir.is_dir() else ligands

    n_workers = workers or (os.cpu_count() or 1)
    n_library = library_size or len(list(iter_ligands(pdbqt_dir)))

    click.echo(f"box: {box}")
    for warning in box.warnings():
        click.secho(f"warning: {warning}", fg="yellow")
    click.echo(f"timing {sample} ligands on {n_workers} workers…\n")

    results = []
    for exhaustiveness in exhaustiveness_values or (8,):
        config = DockConfig(
            receptor=str(receptor), box=box, engine=engine,
            exhaustiveness=exhaustiveness, seed=seed, cpu=1,
        )
        results.append(
            benchmark(
                pdbqt_dir, config,
                label=f"{engine} exh={exhaustiveness}",
                n=sample, n_workers=n_workers, seed=seed,
            )
        )

    click.echo(format_comparison(results, n_library, n_workers))


@main.command("report")
@click.argument("run_dir", type=click.Path(exists=True, path_type=Path))
@click.option(
    "-m", "--manifest", type=click.Path(exists=True, path_type=Path), default=None,
    help="Ligand manifest.csv, to join descriptors onto the results.",
)
@click.option("--flat", is_flag=True, help="Do not group stereoisomers.")
@click.option("--top", default=20, show_default=True, help="Rows to print.")
def report_cmd(run_dir: Path, manifest: Path | None, flat: bool, top: int) -> None:
    """Rebuild results.csv from a run's journal.

    Safe to run against a screen still in progress: it reports on whatever has
    finished so far.
    """
    from drydock.core.results import write_results

    results, all_modes, n = write_results(run_dir, manifest, group_stereoisomers=not flat)
    click.echo(f"wrote {results} ({n} rows) and {all_modes}")

    import csv as _csv

    with open(results, newline="", encoding="utf-8") as fh:
        rows = list(_csv.DictReader(fh))[:top]
    if not rows:
        return

    click.echo()
    click.echo(f"{'rank':>4}  {'compound':<16} {'affinity':>9} {'LE':>7} {'MW':>8}  formula")
    for row in rows:
        le = row.get("ligand_efficiency") or ""
        click.echo(
            f"{row['rank']:>4}  {row['compound_id']:<16} "
            f"{row['best_affinity'] or '-':>9} {le[:6]:>7} "
            f"{(row.get('mw') or '')[:7]:>8}  {row.get('formula', '')}"
        )


@main.command("add-zinc-pseudo")
@click.argument("receptor", type=click.Path(exists=True, path_type=Path))
@click.option("-o", "--out", type=click.Path(path_type=Path), default=None, help="Output PDBQT.")
def add_zinc_pseudo(receptor: Path, out: Path | None) -> None:
    """Add tetrahedral zinc pseudo-atoms for AutoDock4Zn.

    This is the one receptor modification Drydock performs, and it is opt-in
    because the AutoDock4Zn scoring path cannot run without it.

    AutoDock represents zinc coordination through pseudo-atoms (type TZ) placed
    at the *vacant* tetrahedral positions around each zinc. A zinc already
    saturated by the protein gets none, correctly -- there is no site there for a
    ligand to reach.
    """
    from drydock.core.zinc import ZincError, add_zinc_pseudo_atoms

    try:
        result = add_zinc_pseudo_atoms(receptor, out)
    except ZincError as exc:
        raise click.ClickException(str(exc)) from exc

    click.echo(f"zinc atoms found:     {result.n_zinc}")
    click.echo(f"pseudo-atoms placed:  {result.n_pseudo_atoms}")
    click.echo(f"wrote:                {result.receptor_tz}")
    if result.n_pseudo_atoms < result.n_zinc:
        click.echo(
            "note: fewer pseudo-atoms than zincs, which is expected when a zinc is "
            "fully coordinated by the protein (a structural rather than catalytic "
            "site) and so offers a ligand nowhere to bind."
        )


@main.command("maps")
@click.option(
    "-r", "--receptor", required=True, type=click.Path(exists=True, path_type=Path),
    help="Receptor PDBQT, including TZ atoms if using AutoDock4Zn.",
)
@click.option(
    "-o", "--out", required=True, type=click.Path(path_type=Path),
    help="Directory to write maps into.",
)
@click.option("--residues", default=None, help="Residues defining the box.")
@click.option("--chain", default=None, help="Restrict residue selection to one chain.")
@click.option("--pad", default=5.0, show_default=True, help="Box padding in Angstroms.")
@click.option("--center", default=None, help="Explicit box centre as x,y,z.")
@click.option("--size", default=None, help="Explicit box size as x,y,z.")
@click.option(
    "--parameters", type=click.Path(exists=True, path_type=Path), default=None,
    help="AD4 parameter file. Defaults to the bundled AD4Zn.dat.",
)
def maps_cmd(
    receptor: Path,
    out: Path,
    residues: str | None,
    chain: str | None,
    pad: float,
    center: str | None,
    size: str | None,
    parameters: Path | None,
) -> None:
    """Compute AutoGrid maps for the ad4 engine.

    Maps depend on the receptor and the box, not on any individual ligand, so
    they are computed once and reused for the whole screen.
    """
    import shutil

    from drydock.core.box import Box
    from drydock.core.receptor import box_from_residues, inspect
    from drydock.core.zinc import ZincError, run_autogrid, write_gpf

    if center and size:
        box = Box(
            tuple(float(v) for v in center.split(",")),
            tuple(float(v) for v in size.split(",")),
        )
    elif residues:
        numbers = [int(v) for v in residues.replace(" ", "").split(",") if v]
        try:
            box, _ = box_from_residues(receptor, numbers, padding=pad, chain=chain)
        except ValueError as exc:
            raise click.ClickException(str(exc)) from exc
    else:
        raise click.ClickException("give either --residues or both --center and --size")

    report = inspect(receptor)
    if report.metals and not report.has_zinc_pseudo_atoms:
        click.secho(
            "note: this receptor has metals but no TZ pseudo-atoms. For zinc "
            "metalloproteins run 'drydock add-zinc-pseudo' first, or the maps will "
            "not carry AutoDock4Zn's coordination geometry.",
            fg="yellow",
        )

    out.mkdir(parents=True, exist_ok=True)
    local_receptor = out / receptor.name
    if receptor.resolve() != local_receptor.resolve():
        shutil.copy(receptor, local_receptor)

    click.echo(f"box: {box}")
    gpf = write_gpf(local_receptor, box, out / "receptor.gpf", parameter_file=parameters)
    click.echo(f"grid parameter file: {gpf}")

    click.echo("running autogrid4…")
    try:
        maps_dir = run_autogrid(gpf, out)
    except ZincError as exc:
        raise click.ClickException(str(exc)) from exc

    written = sorted(maps_dir.glob("*.map"))
    total_mb = sum(m.stat().st_size for m in written) / 1e6
    click.echo(f"wrote {len(written)} maps ({total_mb:.0f} MB) in {maps_dir}")
    click.echo(f"screen with: drydock screen --engine ad4 --maps {maps_dir} …")


@main.command("check-receptor")
@click.argument("receptor", type=click.Path(exists=True, path_type=Path))
def check_receptor(receptor: Path) -> None:
    """Check a prepared receptor for problems that would distort results.

    Receptor preparation fails in ways that produce a perfectly valid file. The
    two worth checking every time: missing polar hydrogens, which removes every
    hydrogen-bond donor from the protein, and atom types in the wrong case, which
    AutoDock Vina rejects outright.
    """
    from drydock.core.receptor import inspect

    report = inspect(receptor)

    click.echo(f"atoms:          {report.n_atoms}")
    click.echo(f"chains:         {', '.join(report.chains) or '-'}")
    if report.residue_range:
        click.echo(f"residues:       {report.residue_range[0]}-{report.residue_range[1]}")
    click.echo(f"polar H (HD):   {report.n_polar_hydrogens}")
    click.echo(f"metals:         {', '.join(report.metals) or 'none'}")
    click.echo(
        "atom types:     "
        + ", ".join(f"{t or '(none)'}={n}" for t, n in report.atom_types.items())
    )

    for note in report.notes:
        click.echo(f"\nnote: {note}")
    for problem in report.problems:
        click.secho(f"\nPROBLEM: {problem}", fg="red")

    if report.ok:
        click.secho("\nreceptor looks usable", fg="green")
    else:
        raise SystemExit(1)


@main.command("box")
@click.option(
    "-r", "--receptor", required=True, type=click.Path(exists=True, path_type=Path),
    help="Prepared receptor PDBQT.",
)
@click.option(
    "--residues", default=None,
    help="Comma-separated residue numbers lining the site, e.g. 187,188,401,405.",
)
@click.option("--chain", default=None, help="Restrict residue selection to one chain.")
@click.option("--pad", default=5.0, show_default=True, help="Angstroms added on every side.")
@click.option("--sidechains-only", is_flag=True, help="Ignore backbone atoms, tightening the box.")
@click.option("--cubic", is_flag=True, help="Force equal box dimensions.")
@click.option("--center", default=None, help="Explicit centre as x,y,z (instead of --residues).")
@click.option("--size", default=None, help="Explicit size as x,y,z (instead of --residues).")
@click.option("-o", "--out", type=click.Path(path_type=Path), default=None, help="Write a config.")
def box_cmd(
    receptor: Path,
    residues: str | None,
    chain: str | None,
    pad: float,
    sidechains_only: bool,
    cubic: bool,
    center: str | None,
    size: str | None,
    out: Path | None,
) -> None:
    """Define the search box, explicitly or from active-site residues.

    There is no pocket finder: when the site is known, guessing at it adds a
    failure mode without adding information.
    """
    from drydock.core.box import Box
    from drydock.core.receptor import box_from_residues, read_pdbqt

    if center and size:
        box = Box(
            tuple(float(v) for v in center.split(",")),
            tuple(float(v) for v in size.split(",")),
        )
        selected = []
    elif residues:
        numbers = [int(v) for v in residues.replace(" ", "").split(",") if v]
        try:
            box, selected = box_from_residues(
                receptor, numbers, padding=pad, chain=chain,
                sidechains_only=sidechains_only, cubic=cubic,
            )
        except ValueError as exc:
            raise click.ClickException(str(exc)) from exc
    else:
        raise click.ClickException("give either --residues or both --center and --size")

    click.echo(str(box))
    if selected:
        found = sorted({a.residue_seq for a in selected})
        click.echo(f"derived from {len(selected)} atoms of residues {found}")

    # Anything else inside the box competes with the intended site for poses.
    metals_inside = [
        a for a in read_pdbqt(receptor) if a.is_metal and box.contains(a.coordinates)
    ]
    if metals_inside:
        click.echo(
            "metals inside the box: "
            + ", ".join(f"{a.atom_type}{a.residue_seq}" for a in metals_inside)
        )
        if len(metals_inside) > 1:
            click.secho(
                "note: more than one metal site falls inside this box, so ligands "
                "may dock at a site other than the intended one. Consider reducing "
                "--pad or using --sidechains-only.",
                fg="yellow",
            )

    for warning in box.warnings():
        click.secho(f"warning: {warning}", fg="yellow")

    if out:
        out.write_text(box.to_vina_config(receptor=str(receptor)), encoding="utf-8")
        click.echo(f"wrote {out}")


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
