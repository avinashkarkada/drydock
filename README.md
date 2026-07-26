# Drydock

Reproducible high-throughput virtual screening for Linux: ligand preparation and
batch molecular docking, driven from a CLI or a GUI.

Drydock does two things.

**Prepare ligands.** Stream a compound library (SDF, SMILES, MOL2), split it,
sanitise it, set protonation states at a target pH, optimise geometry, and emit
docking-ready PDBQT files alongside a descriptor manifest.

**Screen.** Dock that library against a receptor you prepared, and produce a
ranked CSV with affinities joined to compound descriptors.

Receptor preparation is deliberately **not** part of Drydock. Point it at a
receptor PDBQT you prepared with your tool of choice. The single exception is
zinc pseudo-atom placement, which is opt-in, because the AutoDock4Zn scoring
path cannot work without it.

## Why another docking wrapper

Existing graphical tools — [AMDock](https://github.com/Valdes-Tresanco-MS/AMDock)
in particular, which Drydock borrows from and credits below — are built around
docking one complex interactively. That design does not survive a library of
tens of thousands of compounds: results accumulate in memory, closing the window
kills the run, and a crash at hour fourteen loses everything.

Drydock inverts the relationship. Screening runs as a **detached process** that
writes every finished ligand straight to disk. The GUI is a read-only monitor
that can be closed and reopened at will while the run continues. A killed run
resumes and loses at most one ligand.

## Reproducibility

The goal is that a run on your machine and a run on someone else's produce the
same numbers.

- **Pinned environment.** Every dependency is locked to an exact version, build
  and hash for `linux-64` in `pixi.lock`. Pixi bootstraps its own Python, so
  nothing depends on the host system's interpreter.
- **Deterministic search.** AutoDock Vina is stochastic, and is *not*
  bit-reproducible across threads even with a fixed seed, because thread
  scheduling perturbs the search. Drydock therefore runs every docking job
  single-threaded and parallelises across ligands instead. Per-ligand seeds are
  derived deterministically from one global seed and recorded in the output.
- **Provenance.** Each run writes a manifest recording tool version, resolved
  package versions, engine, scoring function, seed, box, and SHA-256 checksums
  of the receptor and ligand library.
- **Verification.** `drydock selftest` runs the full pipeline against a bundled
  mini test case and checks the affinities against known-good values. CI runs it
  on a clean machine on every commit.

## Install

```bash
git clone https://github.com/<user>/drydock.git
cd drydock
pixi install
pixi run selftest
```

Pixi is the supported install path ([install
pixi](https://pixi.sh/latest/#installation) if you do not have it). It needs no
root and no daemon.

For HPC clusters, a container is available. Build it once and convert to an
Apptainer image, which — unlike Docker — runs rootless on shared systems:

```bash
docker build -t drydock:0.1.0 .
apptainer build drydock.sif docker-daemon://drydock:0.1.0
```

## Quickstart

```bash
# 1. Prepare the ligand library (streams; never loads the whole file)
drydock prep-ligands --input library.sdf --out ligands/ --ph 7.4 --optimize mmff94

# 2. Define the search box, either explicitly or from active-site residues
drydock box --receptor receptor.pdbqt --residues 187,188,401,405,411,420,421,423 \
            --pad 8 --out config.toml

# 3. Measure throughput before committing to a long run
drydock benchmark --config config.toml --n 100

# 4. Screen
drydock screen --config config.toml

# or drive the whole thing graphically
drydock-gui
```

## Engines

| Engine | Scoring | Zinc handling | Speed | Suits |
|---|---|---|---|---|
| `vina` | Vina | Generic; metals typed as H-bond donors | Fastest | General screening |
| `vinardo` | Vinardo | As above | Fastest | General screening |
| `ad4` | AutoDock4 / AD4Zn | **Full tetrahedral zinc pseudo-atoms** | ~Vina | **Metalloproteins** |
| `autodock4` | AutoDock4 / AD4Zn | Full | 10–50× slower | Small focused sets |

`ad4` runs the AutoDock4 scoring function over pre-computed autogrid4 maps while
using Vina's Monte Carlo search. For a metalloprotein target that is usually the
right trade: AutoDock4Zn's coordination geometry at roughly Vina's throughput.

Custom atom parameters (ions and other non-standard atoms) only exist on the
AutoDock4 scoring path — Vina's native function has no user-editable parameter
file. Supplying a custom parameter file therefore implies `ad4`.

`autodock4` is included for small focused sets. Drydock warns if you point it at
a large library.

## Output

`results.csv` — one ranked row per compound:

```
compound_id, rank, best_affinity, ligand_efficiency, smiles, mw, formula,
heavy_atoms, rot_bonds, clogp, tpsa, hbd, hba, formal_charge, n_modes,
seed, engine, status
```

`results_all_modes.csv` — every pose, in the schema PaDEL-ADV emitted, so
existing downstream analyses keep working:

```
Ligand, Mode, Affinity, RMSD_LB, RMSD_UB, RandomSeed
```

## Ring conformations: a caveat worth reading

Vina samples rotatable bonds but **never samples ring geometry**. Whatever ring
conformation is in the input PDBQT is held rigid for the entire run. A sugar
that arrives in a strained boat stays a boat; an extended macrocycle can never
curl to fit a pocket.

`--optimize mmff94` (the default) relaxes each molecule to its *nearest* energy
minimum, which fixes obviously bad geometry but will not flip a chair to a
twist-boat. Meeko's macrocycle handling, also on by default, lets Vina open and
re-close macrocyclic rings. For libraries where ring conformation is critical,
`--conformers N` generates N distinct starting geometries per compound and docks
each, at N times the cost.

## Credits

Drydock stands on:

- **AutoDock Vina** — Trott & Olson, *J Comput Chem* 31:455 (2010); Eberhardt,
  Santos-Martins, Tillack & Forli, *J Chem Inf Model* 61:3891 (2021).
- **AutoDock4 / AutoGrid4** — Morris *et al.*, *J Comput Chem* 30:2785 (2009).
- **AutoDock4Zn** — Santos-Martins, Forli, Ramos & Olson, *J Chem Inf Model*
  54:2371 (2014). The `zinc_pseudo.py` implementation and `AD4Zn.dat` parameter
  file are vendored from AMDock under GPL-3.
- **AMDock** — Valdés-Tresanco, Valdés-Tresanco, Valiente & Moreno, *Biol Direct*
  15:12 (2020). Drydock's zinc handling comes directly from AMDock's data files,
  and its engine layout was informed by AMDock's design.
- **Meeko** and **Scrubber** — Forli lab, Scripps Research.
- **RDKit** — open-source cheminformatics.

## License

GPL-3.0-or-later. See [LICENSE](LICENSE).

Drydock vendors GPL-3 code from AMDock, so the combined work is GPL-3. In
practice this constrains nothing about research use: use it, publish with it,
and share it — the only obligation is that anyone you distribute it to also gets
the source.
