# The run directory

Everything a screening run produces lands in one directory, and everything that
reads a run reads it from there. This is the interface between the screening
process and anything watching it, so it is specified rather than incidental.

```
run/
├── config.toml            the run configuration
├── provenance.json        versions, seed, box, input checksums
├── journal.jsonl          append-only, one record per finished ligand
├── status.json            aggregate summary, rewritten periodically
├── poses/<id>.pdbqt       exported poses, top-N per ligand
├── logs/                  engine stderr worth keeping
├── results.csv            ranked, one row per compound (derived)
└── results_all_modes.csv  every pose, PaDEL-ADV schema (derived)
```

## Why it is shaped this way

Screening a large library takes hours to days. Over that window the process will
be interrupted, a closed laptop lid, an OOM kill, a cluster preemption, an
impatient `Ctrl-C`. Two decisions follow.

**The journal is the authority, and it is append-only.** Every finished ligand is
written to `journal.jsonl` and `fsync`'d before the run moves on. Exactly one
process writes it: the parent, which collects finished work from its worker pool.
Workers never touch it, so there is no locking and no interleaved lines.

The cost is one `fsync` per ligand. For work measured in seconds per ligand that
is noise, and it buys an honest guarantee: **a killed run loses at most the
ligand that was in flight.**

**Nothing important lives in memory.** Because progress is on disk, the process
doing the docking and the process watching it are fully independent. The GUI is a
read-only watcher; closing it does not touch the run, and reopening it recovers
the full picture. This is the main structural difference from interactive docking
tools, where results accumulate in the UI process and closing the window ends
the run.

## `journal.jsonl`

One JSON object per line. Successful ligand:

```json
{"ligand_id":"CMNPD1","status":"ok","seed":1943860988,"elapsed_s":4.21,
 "timestamp":1753500000.0,
 "modes":[{"mode":1,"affinity":-7.5,"rmsd_lb":0.0,"rmsd_ub":0.0},
          {"mode":2,"affinity":-7.4,"rmsd_lb":2.219,"rmsd_ub":3.759}]}
```

Failed ligand:

```json
{"ligand_id":"CMNPD2","status":"failed","seed":42,"elapsed_s":0.3,
 "timestamp":1753500001.0,"error":"engine returned no poses"}
```

`status` is `ok`, `failed`, or `skipped`. Field names under `modes` match what
Vina prints, and pass through to `results_all_modes.csv` unchanged.

Two reader rules, both load-bearing:

- **A truncated final line is discarded, not an error.** That is the signature of
  a killed writer, and that run is exactly the one being recovered.
- **Failures count as done.** A compound that reliably crashes the engine is
  recorded and never retried, so it cannot trap a resumed run in a loop.

## `status.json`

A cache, never authoritative, it holds only what can be recomputed from the
journal. It exists so a watcher polling once a second does not re-parse a
47,000-line file each time.

Written atomically (temp file, `fsync`, `os.replace`), so a reader observes
either the old file or the new one, never a half-written one. A corrupt or
missing `status.json` costs a rebuild, not a run: readers fall back to
`RunDir.rebuild_status()`.

## Watching a run

Two mechanisms, matched to what a watcher needs:

```python
from drydock.core.rundir import RunDir

run = RunDir("runs/mmp9")

# Aggregate progress: one small file read.
status = run.read_status() or run.rebuild_status()
print(status.done, status.total, status.eta_s)

# Individual records, incrementally: constant time per poll regardless of
# how long the run has been going.
offset = 0
while True:
    new_records, offset = run.tail_journal(offset)
    for record in new_records:
        print(record.ligand_id, record.best_affinity)
```

`tail_journal` consumes only whole lines. A record still being written is left
for the next call rather than parsed half-formed, and the returned offset does
not advance past it. If the journal shrinks, replaced or deleted underneath the
watcher, it restarts from the top rather than reading from a now-meaningless
position.

## Resuming

`RunDir.completed_ids()` returns every ligand ID in the journal, successes and
failures alike. A resumed run docks the set difference against the library.
Re-running a finished run is therefore a no-op, and re-running an interrupted one
picks up exactly where it stopped.

## Schema version

`SCHEMA_VERSION` is stamped into `provenance.json` and `status.json`, and is
bumped when the layout changes incompatibly, so a future reader can refuse an
unfamiliar run rather than silently misinterpreting it.

## Fabricating one

Run directories can be generated without an engine, which is how the GUI is
developed and tested:

```bash
drydock dev synthesize-run /tmp/demo -n 5000            # a finished run
drydock dev synthesize-run /tmp/demo -n 500 --live      # updates as you watch
drydock dev synthesize-run /tmp/demo -n 500 --incomplete # looks like a killed run
```

Output is deterministic for a given `--seed`.
