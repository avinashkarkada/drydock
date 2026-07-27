# Installing Drydock

Linux only. Takes about five minutes, most of it downloading.

## 1. Install pixi

Pixi manages the environment. It needs no root and no daemon, and installs
everything — Python, AutoDock Vina, AutoGrid, RDKit — into a folder inside this
directory rather than touching your system.

```bash
curl -fsSL https://pixi.sh/install.sh | bash
```

Open a new terminal afterwards, or `source ~/.bashrc`, so `pixi` is on your PATH.

## 2. Install Drydock

From inside this directory:

```bash
pixi install
```

This reads `pixi.lock` and fetches the exact versions it names — same versions,
same builds, same checksums as the machine it was developed on. Roughly 2 GB.

## 3. Check it works

```bash
pixi run drydock selftest
```

This is not a smoke test. It runs the real pipeline — prepares five ligands,
docks them against a bundled receptor, ranks them — and compares the affinities
against values recorded from a known-good run. It takes about 20 seconds and
should end with:

```
self-test passed
```

If it does, your installation produces the same numbers as everyone else's. If
it fails, the output names which check failed and by how much; send that rather
than a description.

## 4. Run it

```bash
pixi run drydock-gui        # graphical
pixi run drydock --help     # command line
```

The GUI has four tabs, in the order the work is done: **1. Receptor**,
**2. Ligands**, **3. Screen**, **4. Results**. Each hands its output to the next.

## What you need to supply

- **A receptor structure** — PDB or mmCIF. Drydock prepares it for you.
- **A compound library** — SDF, SMILES or MOL2, optionally gzipped.
- **The active-site residues**, or an explicit box centre and size.

## Things that will save you time

**Check the receptor before screening.** Preparation can succeed and still
produce a receptor that quietly ruins a screen — most commonly one with no polar
hydrogens, which removes every hydrogen-bond donor from the protein without
raising anything. `drydock check-receptor` catches that class of problem, and the
Receptor tab runs it automatically.

**Benchmark before a long run.** Cost per ligand varies by more than tenfold with
box size, ligand flexibility and exhaustiveness. `drydock benchmark` docks a
sample and projects the full run, which is worth five minutes before committing
to something that may take days. Run it on an otherwise idle machine — a
benchmark saturates every core it is given, so anything else running inflates it.

**Survey the library first.** `drydock survey` reports records against distinct
compounds. These differ whenever a source enumerates stereoisomers under one
identifier, and the gap can be large: CMNPD 1.0 has 47,451 records under 25,224
identifiers.

**Runs are detached and resumable.** Closing the GUI does not stop a screen, and
re-running the same command after an interruption picks up where it stopped. A
killed run loses at most the ligand in flight.

## For a cluster

Build a container once and convert it to Apptainer, which runs rootless where
Docker usually cannot:

```bash
docker build -t drydock:0.1.0 .
apptainer build drydock.sif docker-daemon://drydock:0.1.0
```

## Licence

GPL-3.0-or-later. Use it, publish with it, share it — the only obligation is that
anyone you pass it to also gets the source, which is what this archive is.

Please cite the underlying engines as well as Drydock; see `CITATION.cff`.
