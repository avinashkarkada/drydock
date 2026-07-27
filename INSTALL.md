# Installing Drydock

Linux only. Takes about five minutes, most of it downloading.

## 1. Install pixi

Pixi manages the environment. It needs no root and no daemon. Everything it
installs (Python, AutoDock Vina, AutoGrid, RDKit) goes into a folder inside this
directory rather than into your system.

```bash
curl -fsSL https://pixi.sh/install.sh | bash
```

Open a new terminal afterwards, or `source ~/.bashrc`, so `pixi` is on your PATH.

## 2. Install Drydock

From inside this directory:

```bash
pixi install
```

This reads `pixi.lock` and fetches the exact versions listed there: same
versions, same builds, same checksums as the machine it was developed on. About
2 GB.

## 3. Check it works

```bash
pixi run drydock selftest
```

This is not a smoke test. It runs the real pipeline, preparing five ligands,
docking them against a bundled receptor and ranking them, then compares the
affinities against values recorded from a known-good run. Takes about 20 seconds
and should end with:

```
self-test passed
```

If it does, your installation produces the same numbers as everyone else's. If
it fails, the output names which check failed and by how much. Send that rather
than a description.

## 4. Run it

```bash
pixi run drydock-gui        # graphical
pixi run drydock --help     # command line
```

The GUI has four tabs, in the order the work is done: **1. Receptor**,
**2. Ligands**, **3. Screen**, **4. Results**. Each hands its output to the next.

## What you need to supply

- **A receptor structure** in PDB or mmCIF format. Drydock prepares it for you.
- **A compound library** as SDF, SMILES or MOL2, optionally gzipped.
- **The active-site residues**, or an explicit box centre and size.

## Things that will save you time

**Check the receptor before screening.** Preparation can succeed and still give
you a receptor that ruins a screen. The usual culprit is missing polar hydrogens,
which removes every hydrogen-bond donor from the protein without raising an
error. `drydock check-receptor` catches that, and the Receptor tab runs it for
you.

**Benchmark before a long run.** Cost per ligand varies by more than tenfold with
box size, ligand flexibility and exhaustiveness. `drydock benchmark` docks a
sample and projects the full run. Five minutes well spent before committing to
something that might take days. Run it on an idle machine: a benchmark uses every
core it is given, so anything else running will inflate the numbers.

**Survey the library first.** `drydock survey` reports records against distinct
compounds. The two differ whenever a source enumerates stereoisomers under one
identifier, and the gap can be large. CMNPD 1.0 has 47,451 records under 25,224
identifiers.

**Runs are detached and resumable.** Closing the GUI does not stop a screen, and
re-running the same command after an interruption picks up where it left off. A
killed run loses at most the ligand that was in flight.

## For a cluster

Build a container once and convert it to Apptainer, which runs rootless where
Docker usually cannot:

```bash
docker build -t drydock:0.1.0 .
apptainer build drydock.sif docker-daemon://drydock:0.1.0
```

## Licence

GPL-3.0-or-later. Use it, publish with it and share it. The only obligation is
that anyone you pass it to also gets the source, which is what this archive is.

Please cite the underlying engines as well as Drydock. See `CITATION.cff`.
