"""Streaming readers for compound libraries.

Every reader here is a generator that holds one record in memory at a time. The
library that motivated Drydock is a 174 MB SDF of 47,451 marine natural products;
loading it with ``Chem.SDMolSupplier`` would take gigabytes and minutes before any
work started. Streaming means preparation begins on the first compound and peak
memory does not depend on library size.

Records are yielded as :class:`Record` -- identifiers plus the raw text block --
rather than as RDKit molecules, because these get handed to worker processes and
strings pickle cheaply where molecules do not.

Identifiers: two of them
------------------------

Compound identifiers in real libraries are **not unique**, and treating them as
filenames silently destroys data. CMNPD 1.0 is a working example: 47,451 records
carry only 25,224 distinct ``COMPOUND_ID`` values, because every compound with
undefined stereocentres has had all 2^n stereoisomers enumerated under one ID
(``CMNPD22318`` appears 64 times). Each is a genuinely different molecule that
docks differently, so they all deserve to be screened -- but writing them all to
``CMNPD22318.pdbqt`` would keep one and lose 63.

So every record carries both:

``compound_id``
    Whatever the file said. Shared across stereoisomers, and the right thing to
    group a hit list by.
``ligand_id``
    Unique within the library, derived by suffixing repeats (``CMNPD30``,
    ``CMNPD30_2``, ...). Safe as a filename and as a dictionary key.
"""

from __future__ import annotations

import gzip
import re
from collections import defaultdict
from collections.abc import Iterator
from dataclasses import dataclass
from pathlib import Path
from typing import IO

# Property names that commonly hold a compound identifier, most specific first.
# CMNPD uses COMPOUND_ID; other sources vary, and guessing well saves the user
# from discovering --id-field on their first run.
ID_FIELDS: tuple[str, ...] = (
    "COMPOUND_ID",
    "compound_id",
    "IDNUMBER",
    "ID",
    "id",
    "Name",
    "NAME",
    "name",
    "ChEMBL_ID",
    "chembl_id",
    "PUBCHEM_COMPOUND_CID",
    "CATALOG_ID",
    "Catalog_ID",
    "ZINC_ID",
    "zinc_id",
)

# Characters that make an identifier unusable as a filename. Replaced rather than
# rejected, since a compound should not be dropped over punctuation.
_UNSAFE_CHARS = re.compile(r"[^\w.\-+]")

# A byte-order mark leading the file is not part of the first record's title.
_BOM = "﻿"


@dataclass(frozen=True, slots=True)
class Record:
    """One compound as it appeared in the library file."""

    ligand_id: str
    """Unique within the library. Safe to use as a filename."""

    compound_id: str
    """As found in the file. May be shared by stereoisomers of one compound."""

    block: str
    """The raw record text, ready to hand to RDKit."""

    fmt: str
    """``sdf``, ``smi`` or ``mol2``."""

    index: int
    """Zero-based position in the file, for reporting parse failures."""

    @property
    def is_duplicate_id(self) -> bool:
        """True when this record had to be renamed to avoid a collision."""
        return self.ligand_id != self.compound_id


class LibraryFormatError(ValueError):
    """Raised when a library file's format cannot be determined or parsed."""


def _open_text(path: Path) -> IO[str]:
    """Open plain or gzipped text, since compound libraries often ship as .gz."""
    if path.suffix == ".gz":
        return gzip.open(path, "rt", encoding="utf-8", errors="replace")
    return open(path, encoding="utf-8", errors="replace")


def detect_format(path: str | Path) -> str:
    """Infer library format from the filename, ignoring any .gz wrapper."""
    name = Path(path).name.lower()
    if name.endswith(".gz"):
        name = name[:-3]
    if name.endswith((".sdf", ".sd", ".mol")):
        return "sdf"
    if name.endswith((".smi", ".smiles", ".ism", ".txt")):
        return "smi"
    if name.endswith(".mol2"):
        return "mol2"
    raise LibraryFormatError(
        f"cannot determine format of {path!r}; expected .sdf, .smi, .mol2 (optionally .gz)"
    )


def sanitize_id(value: str) -> str:
    """Make an identifier safe to use as a filename.

    Path separators and other awkward characters become underscores. An
    identifier that reduces to nothing is rejected by the caller rather than
    silently becoming an empty filename.
    """
    cleaned = _UNSAFE_CHARS.sub("_", value.strip().lstrip(_BOM))
    return cleaned.strip("_")


def iter_library(
    path: str | Path,
    fmt: str | None = None,
    id_field: str | None = None,
    unique_ids: bool = True,
) -> Iterator[Record]:
    """Stream records from a compound library.

    Args:
        path: Library file, optionally gzipped.
        fmt: ``sdf``, ``smi`` or ``mol2``. Inferred from the filename if omitted.
        id_field: Property holding the compound ID. If omitted, the reader tries
            the names in :data:`ID_FIELDS`, then the record title, then falls
            back to a positional ID.
        unique_ids: Suffix repeated identifiers so ``ligand_id`` is unique. Turn
            off only if you have independently guaranteed uniqueness; duplicates
            will otherwise overwrite one another downstream.

    Yields:
        One :class:`Record` per compound, in file order.
    """
    path = Path(path)
    fmt = fmt or detect_format(path)

    if fmt == "sdf":
        raw = _iter_sdf(path, id_field)
    elif fmt == "smi":
        raw = _iter_smi(path)
    elif fmt == "mol2":
        raw = _iter_mol2(path)
    else:
        raise LibraryFormatError(f"unsupported format {fmt!r}")

    if not unique_ids:
        yield from raw
        return

    # Counts of identifiers seen so far. Holding these is the one part of reading
    # that is not constant-memory, but it is a few hundred thousand short strings
    # even for a large library -- trivial next to the cost of getting it wrong.
    seen: defaultdict[str, int] = defaultdict(int)
    for record in raw:
        seen[record.compound_id] += 1
        occurrence = seen[record.compound_id]
        if occurrence == 1:
            yield record
        else:
            yield Record(
                ligand_id=f"{record.compound_id}_{occurrence}",
                compound_id=record.compound_id,
                block=record.block,
                fmt=record.fmt,
                index=record.index,
            )


def _iter_sdf(path: Path, id_field: str | None) -> Iterator[Record]:
    """Stream an SDF by splitting on the ``$$$$`` record terminator.

    Parsed as text rather than through RDKit's supplier so that a molecule RDKit
    rejects still yields a record with an identifier -- which is what lets the
    failure be reported against a name the user recognises instead of an index.
    """
    with _open_text(path) as fh:
        buffer: list[str] = []
        index = 0
        for line in fh:
            if line.startswith("$$$$"):
                block = "".join(buffer)
                buffer.clear()
                if block.strip():
                    yield _make_sdf_record(block, id_field, index)
                    index += 1
            else:
                buffer.append(line)

        # A final record with no terminator is malformed but recoverable.
        block = "".join(buffer)
        if block.strip():
            yield _make_sdf_record(block, id_field, index)


def _make_sdf_record(block: str, id_field: str | None, index: int) -> Record:
    compound_id = _sdf_id(block, id_field, index)
    return Record(compound_id, compound_id, block, "sdf", index)


def _sdf_id(block: str, id_field: str | None, index: int) -> str:
    """Pull an identifier out of an SDF record.

    Order: the requested field, then any recognised property name, then the title
    line, then a positional fallback.

    Properties are tried before the title deliberately. The title is where the
    format intends the name to live, but in practice bulk-generated SDFs leave it
    blank or fill it with the generating program's banner, whereas an explicit
    ``> <COMPOUND_ID>`` tag is there because someone meant it.
    """
    if id_field:
        if value := _sdf_property(block, id_field):
            return sanitize_id(value) or f"mol{index + 1}"

    for field in ID_FIELDS:
        if value := _sdf_property(block, field):
            if cleaned := sanitize_id(value):
                return cleaned

    lines = block.splitlines()
    if lines:
        title = lines[0].strip().lstrip(_BOM).strip()
        # Guard against the byte-order mark, blank titles, and the writer banners
        # that occupy this line in machine-generated files.
        if title and len(title) <= 120 and not title.startswith(("#", "$")):
            if cleaned := sanitize_id(title):
                return cleaned

    return f"mol{index + 1}"


_PROP_RE_CACHE: dict[str, re.Pattern[str]] = {}


def _sdf_property(block: str, field: str) -> str | None:
    """Read one ``> <FIELD>`` property value from an SDF record."""
    pattern = _PROP_RE_CACHE.get(field)
    if pattern is None:
        pattern = re.compile(rf"^>\s+<{re.escape(field)}>.*$", re.MULTILINE)
        _PROP_RE_CACHE[field] = pattern

    match = pattern.search(block)
    if not match:
        return None
    rest = block[match.end() :].lstrip("\n")
    value = rest.split("\n", 1)[0].strip()
    return value or None


def _iter_smi(path: Path) -> Iterator[Record]:
    """Stream a SMILES file: one ``SMILES [id]`` per line.

    A header naming the columns is skipped if present, since exports from
    spreadsheets routinely include one and it is not a molecule.
    """
    with _open_text(path) as fh:
        index = 0
        for lineno, raw_line in enumerate(fh):
            line = raw_line.strip().lstrip(_BOM).strip()
            if not line or line.startswith("#"):
                continue
            parts = line.split()
            smiles = parts[0]
            if lineno == 0 and smiles.lower() in ("smiles", "smi", "structure"):
                continue
            ligand_id = sanitize_id(parts[1]) if len(parts) > 1 else ""
            if not ligand_id:
                ligand_id = f"mol{index + 1}"
            yield Record(ligand_id, ligand_id, smiles, "smi", index)
            index += 1


def _iter_mol2(path: Path) -> Iterator[Record]:
    """Stream a multi-molecule MOL2 by splitting on ``@<TRIPOS>MOLECULE``."""
    with _open_text(path) as fh:
        buffer: list[str] = []
        index = 0
        for line in fh:
            if line.startswith("@<TRIPOS>MOLECULE") and buffer:
                block = "".join(buffer)
                buffer.clear()
                yield _make_mol2_record(block, index)
                index += 1
            buffer.append(line)

        if buffer:
            block = "".join(buffer)
            if block.strip():
                yield _make_mol2_record(block, index)


def _make_mol2_record(block: str, index: int) -> Record:
    compound_id = _mol2_id(block, index)
    return Record(compound_id, compound_id, block, "mol2", index)


def _mol2_id(block: str, index: int) -> str:
    """The MOL2 name is the line immediately after ``@<TRIPOS>MOLECULE``."""
    lines = block.splitlines()
    for i, line in enumerate(lines):
        if line.startswith("@<TRIPOS>MOLECULE") and i + 1 < len(lines):
            if cleaned := sanitize_id(lines[i + 1]):
                return cleaned
            break
    return f"mol{index + 1}"


def count_records(path: str | Path, fmt: str | None = None) -> int:
    """Count compounds without parsing them.

    Worth its own pass: knowing the total up front turns an indeterminate
    progress bar into a real one, and the scan is IO-bound rather than
    chemistry-bound, so it costs about a second even on a 174 MB library.
    """
    path = Path(path)
    fmt = fmt or detect_format(path)

    if fmt == "smi":
        with _open_text(path) as fh:
            return sum(1 for line in fh if line.strip() and not line.startswith("#"))

    if fmt == "mol2":
        # The marker opens each molecule, so markers and molecules correspond.
        with _open_text(path) as fh:
            return sum(1 for line in fh if line.startswith("@<TRIPOS>MOLECULE"))

    # In SDF the marker *terminates* a record, so a final record written without
    # one would go uncounted. Track whether anything followed the last marker.
    count = 0
    pending = False
    with _open_text(path) as fh:
        for line in fh:
            if line.startswith("$$$$"):
                count += 1
                pending = False
            elif line.strip():
                pending = True
    return count + (1 if pending else 0)


def survey(path: str | Path, fmt: str | None = None, id_field: str | None = None) -> dict:
    """Summarise a library without preparing it.

    Reports the split between records and distinct compounds, which is the number
    users need before committing to a screen: it is the difference between "how
    long will this take" (records) and "how many compounds am I testing"
    (distinct identifiers).
    """
    total = 0
    counts: defaultdict[str, int] = defaultdict(int)
    for record in iter_library(path, fmt=fmt, id_field=id_field, unique_ids=False):
        total += 1
        counts[record.compound_id] += 1

    repeated = {cid: n for cid, n in counts.items() if n > 1}
    return {
        "records": total,
        "distinct_compounds": len(counts),
        "compounds_with_variants": len(repeated),
        "max_variants": max(repeated.values(), default=1),
        "most_repeated": sorted(repeated.items(), key=lambda kv: -kv[1])[:5],
    }
