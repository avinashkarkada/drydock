# Drydock

Reproducible high-throughput virtual screening for Linux: ligand preparation and
batch molecular docking, driven from a CLI or a GUI.

Drydock does three things.

**Prepare receptors.** Read a PDB or mmCIF, add hydrogens, assign AutoDock atom
types and charges, then check the result. A receptor can prepare cleanly and
still be unusable, so it is worth checking.

**Prepare ligands.** Stream a compound library (SDF, SMILES, MOL2), split it,
sanitise it, set protonation states at a target pH, optimise geometry and write
docking-ready PDBQT files alongside a descriptor manifest.

**Screen.** Dock that library against the receptor and produce a ranked CSV with
affinities joined to compound descriptors.

## Why another docking wrapper

Existing graphical tools are built around docking one complex interactively.
[AMDock](https://github.com/Valdes-Tresanco-MS/AMDock), which Drydock borrows
from and credits below, is the obvious example. That design does not survive a
library of tens of thousands of compounds: results accumulate in memory, closing
the window kills the run, and a crash at hour fourteen loses everything.

Drydock turns that around. Screening runs as a **detached process** that writes
every finished ligand straight to disk. The GUI is a read-only monitor you can
close and reopen while the run continues. A killed run resumes and loses at most
one ligand.

## Reproducibility

The goal is that a run on your machine and a run on someone else's produce the
same numbers.

- **Pinned environment.** Every dependency is locked to an exact version, build
  and hash for `linux-64` in `pixi.lock`. Pixi bootstraps its own Python, so
  nothing depends on the host system's interpreter.
- **Deterministic search.** AutoDock Vina is stochastic, and is *not*
  bit-reproducible across threads even with a fixed seed, because thread
  scheduling perturbs the search. Drydock runs every docking job single-threaded
  and parallelises across ligands instead. Seeds are recorded in the output.
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
Apptainer image, which, unlike Docker, runs rootless on shared systems:

```bash
docker build -t drydock:0.1.0 .
apptainer build drydock.sif docker-daemon://drydock:0.1.0
```

## Quickstart

```bash
# 0. Prepare the receptor (PDB or mmCIF), and check what came out
drydock prep-receptor -i 2OVX.cif -o receptor
drydock check-receptor receptor.pdbqt

# 1. Look at the library before committing to it
drydock survey library.sdf

# 2. Prepare it (streams; never loads the whole file into memory)
drydock prep-ligands --input library.sdf --out ligands/ --ph 7.4

# 3. Define the search box from the residues lining your site
drydock box --receptor receptor.pdbqt --residues 187,188,401,405,411,420,421,423

# 4. Measure throughput on your hardware before starting a long run
drydock benchmark --receptor receptor.pdbqt --ligands ligands/ \
                  --residues 187,188,401,405,411,420,421,423

# 5. Screen
drydock screen --receptor receptor.pdbqt --ligands ligands/ --run runs/mmp9 \
               --residues 187,188,401,405,411,420,421,423

# progress, from anywhere, while it runs
drydock status runs/mmp9

# or drive the whole thing graphically
drydock-gui
```

The GUI has four tabs in the order the work is done: **1. Receptor**,
**2. Ligands**, **3. Screen** and **4. Results**. Each stage hands its output to
the next so paths do not have to be retyped. Ligand preparation and screening
run as detached processes watched through the directories they write, so the
window can be closed and reopened at any point without disturbing a run.

For a zinc metalloprotein, add the AutoDock4Zn steps before screening:

```bash
drydock add-zinc-pseudo receptor.pdbqt          # -> receptor_TZ.pdbqt
drydock maps --receptor receptor_TZ.pdbqt --out maps/ --residues 401,405,411
drydock screen --engine ad4 --maps maps/ ...
```

### Always check the receptor

`prep-receptor` checks its own output, and `check-receptor` will inspect a file
prepared anywhere else. Both exist because receptor preparation fails in ways
that produce a perfectly valid file. Two are common enough to be worth naming:

- **No polar hydrogens.** AutoDock encodes a hydrogen-bond donor as a heavy atom
  with an `HD` hydrogen attached. Strip the hydrogens and the donors do not just
  get weaker, they stop existing: every backbone amide, lysine and tyrosine
  hydroxyl then scores acceptor-only. Nothing errors.
- **Atom types in the wrong case.** AutoDock types are case-sensitive: `Zn`, not
  `ZN`. Vina rejects the wrong case, but reports it as a C++ overload-resolution
  error that reads like a bug in the caller.

Both were present in the first real receptor Drydock was pointed at.

## Engines

| Engine | Scoring | Zinc handling | Speed | Suits |
|---|---|---|---|---|
| `vina` | Vina | Generic; metals typed as H-bond donors | Fastest | General screening |
| `vinardo` | Vinardo | As above | Fastest | General screening |
| `ad4` | AutoDock4 / AD4Zn | **Full tetrahedral zinc pseudo-atoms** | ~Vina | **Metalloproteins** |
| `autodock4` | AutoDock4 / AD4Zn | Full | 10-50x slower | Small focused sets |

`ad4` runs the AutoDock4 scoring function over pre-computed autogrid4 maps while
using Vina's Monte Carlo search. For a metalloprotein target that is usually the
right trade: AutoDock4Zn's coordination geometry at roughly Vina's throughput.

Custom atom parameters (ions and other non-standard atoms) only exist on the
AutoDock4 scoring path. Vina's native function has no user-editable parameter
file, so supplying one implies `ad4`.

`autodock4` is included for small focused sets. Drydock warns if you point it at
a large library.

## Identifiers are not unique

Compound identifiers repeat in real libraries, and a tool that treats them as
filenames will destroy data without reporting anything.

CMNPD 1.0 is a good example. It has **47,451 records under 25,224 distinct
`COMPOUND_ID` values**, because compounds with undefined stereocentres have had
all their stereoisomers enumerated under one identifier. `CMNPD22318` appears 64
times. Each is a different molecule that docks differently, so writing them all
to `CMNPD22318.pdbqt` would keep one and lose 63.

Every record therefore gets two identifiers: a unique `ligand_id` (suffixed on
collision) and the source `compound_id` to group by. `results.csv` ranks by
compound, keeping the best-scoring variant and recording how many were tried.
Use `--flat` to report every record separately.

Run `drydock survey` on a library to see the split before you commit to it.

## Output

`results.csv` has one ranked row per compound:

```
compound_id, rank, best_affinity, ligand_efficiency, smiles, mw, formula,
heavy_atoms, rot_bonds, clogp, tpsa, hbd, hba, formal_charge, n_modes,
seed, engine, status
```

`results_all_modes.csv` has every pose, in the schema PaDEL-ADV used, so
existing downstream analyses keep working:

```
Ligand, Mode, Affinity, RMSD_LB, RMSD_UB, RandomSeed
```

## Ring conformations

Vina samples rotatable bonds but **never samples ring geometry**. Whatever ring
conformation is in the input PDBQT stays rigid for the whole run. A sugar that
arrives in a strained boat stays a boat, and an extended macrocycle can never
curl up to fit a pocket.

The default `--optimize` relaxes each molecule to its *nearest* energy minimum,
which fixes obviously bad geometry but will not flip a chair to a twist-boat.
Meeko's macrocycle handling, also on by default, lets Vina open and re-close
macrocyclic rings. If ring conformation matters for your library, `--conformers
N` generates N starting geometries per compound and docks each, at N times the
cost.

## Credits

Drydock stands on:

- **AutoDock Vina**: Trott & Olson, *J Comput Chem* 31:455 (2010); Eberhardt,
  Santos-Martins, Tillack & Forli, *J Chem Inf Model* 61:3891 (2021).
- **AutoDock4 / AutoGrid4**: Morris *et al.*, *J Comput Chem* 30:2785 (2009).
- **AutoDock4Zn**: Santos-Martins, Forli, Ramos & Olson, *J Chem Inf Model*
  54:2371 (2014). The `zinc_pseudo.py` implementation and `AD4Zn.dat` parameter
  file are vendored from AMDock under GPL-3.
- **AMDock**: Valdés-Tresanco, Valdés-Tresanco, Valiente & Moreno, *Biol Direct*
  15:12 (2020). Drydock's zinc handling comes from AMDock's data files, and its
  engine layout was informed by AMDock's design.
- **Meeko** and **Scrubber**: Forli lab, Scripps Research.
- **RDKit**: open-source cheminformatics.

## License

GPL-3.0-or-later. See [LICENSE](LICENSE).

Drydock vendors GPL-3 code from AMDock, so the combined work is GPL-3. This
constrains nothing about research use. Use it, publish with it and share it. The
only obligation is that anyone you distribute it to also gets the source.
