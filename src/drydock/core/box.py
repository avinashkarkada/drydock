"""The search box.

Docking searches a rectangular volume, and where that volume sits is one of the
few decisions in a screen that can invalidate everything downstream without ever
producing an error. A box in the wrong place returns plausible-looking affinities
for poses in the wrong pocket.

Drydock supports two ways to define it, and deliberately no more:

**Explicit** -- centre and size in Angstroms, as in a Vina ``config.txt``. Use
this when you already know the coordinates, typically from a co-crystal ligand.

**From residues** -- name the residues lining the site and let the box be
computed to enclose them, plus padding. Use this when you know the site
biochemically rather than numerically.

There is no automatic pocket finder. When the site is known, guessing at it adds
a failure mode without adding information.
"""

from __future__ import annotations

from collections.abc import Sequence
from dataclasses import dataclass
from typing import Any

# Vina warns above 27,000 A^3, and it is right to: search difficulty grows with
# volume, so a box much larger than the site wastes exhaustiveness on empty
# space and quietly degrades every result in the run.
LARGE_BOX_VOLUME = 27_000.0

# Padding beyond the selected atoms. A ligand needs room to place substituents
# past the residues that line the pocket; too little and poses are clipped at the
# boundary, which Vina does not report as an error.
DEFAULT_PADDING = 5.0


@dataclass(frozen=True, slots=True)
class Box:
    """A docking search volume, in Angstroms."""

    center: tuple[float, float, float]
    size: tuple[float, float, float]

    def __post_init__(self) -> None:
        if any(s <= 0 for s in self.size):
            raise ValueError(f"box dimensions must be positive, got {self.size}")

    @property
    def volume(self) -> float:
        return self.size[0] * self.size[1] * self.size[2]

    @property
    def is_large(self) -> bool:
        """True if big enough that search quality is likely to suffer."""
        return self.volume > LARGE_BOX_VOLUME

    @property
    def minimum(self) -> tuple[float, float, float]:
        return tuple(c - s / 2 for c, s in zip(self.center, self.size, strict=True))

    @property
    def maximum(self) -> tuple[float, float, float]:
        return tuple(c + s / 2 for c, s in zip(self.center, self.size, strict=True))

    def contains(self, point: Sequence[float]) -> bool:
        return all(
            lo <= p <= hi
            for lo, p, hi in zip(self.minimum, point, self.maximum, strict=True)
        )

    @classmethod
    def from_atoms(
        cls,
        coordinates: Sequence[Sequence[float]],
        padding: float = DEFAULT_PADDING,
        cubic: bool = False,
    ) -> Box:
        """Enclose a set of atoms, with padding.

        Args:
            coordinates: Atom positions to enclose.
            padding: Angstroms added on every side.
            cubic: Force equal dimensions. Occasionally wanted for AutoGrid maps,
                where non-cubic boxes are legal but awkward to reason about.
        """
        if not coordinates:
            raise ValueError("cannot build a box from no atoms")

        lows = [min(c[i] for c in coordinates) for i in range(3)]
        highs = [max(c[i] for c in coordinates) for i in range(3)]

        center = tuple(round((lo + hi) / 2, 3) for lo, hi in zip(lows, highs, strict=True))
        size = [round(hi - lo + 2 * padding, 3) for lo, hi in zip(lows, highs, strict=True)]

        if cubic:
            size = [max(size)] * 3

        return cls(center, tuple(size))

    @classmethod
    def from_config(cls, data: dict[str, Any]) -> Box:
        """Build from a parsed ``[box]`` config table."""
        if "center" in data and "size" in data:
            return cls(
                tuple(float(v) for v in data["center"]),
                tuple(float(v) for v in data["size"]),
            )
        try:
            center = (float(data["center_x"]), float(data["center_y"]), float(data["center_z"]))
            size = (float(data["size_x"]), float(data["size_y"]), float(data["size_z"]))
        except KeyError as exc:
            raise ValueError(f"box definition is missing {exc}") from exc
        return cls(center, size)

    def to_dict(self) -> dict[str, list[float]]:
        return {"center": list(self.center), "size": list(self.size)}

    def to_vina_config(self, receptor: str | None = None) -> str:
        """Render as a Vina ``config.txt``.

        Emitted in the layout AutoDock Vina and PaDEL-ADV both use, so the box
        can be handed to those tools directly.
        """
        lines = []
        if receptor:
            lines.append(f"receptor = {receptor}")
            lines.append("")
        lines += [
            f"center_x = {self.center[0]}",
            f"center_y = {self.center[1]}",
            f"center_z = {self.center[2]}",
            "",
            f"size_x = {self.size[0]}",
            f"size_y = {self.size[1]}",
            f"size_z = {self.size[2]}",
        ]
        return "\n".join(lines) + "\n"

    def warnings(self) -> list[str]:
        """Problems worth surfacing before a long run commits to this box."""
        issues = []
        if self.is_large:
            issues.append(
                f"box volume {self.volume:.0f} A^3 exceeds {LARGE_BOX_VOLUME:.0f} A^3; "
                "search quality degrades in large boxes -- consider tightening it "
                "or raising exhaustiveness"
            )
        if any(s < 10 for s in self.size):
            issues.append(
                f"box dimension {min(self.size):.1f} A is small; ligands larger than "
                "the box cannot be posed and will score poorly for the wrong reason"
            )
        return issues

    def __str__(self) -> str:
        cx, cy, cz = self.center
        sx, sy, sz = self.size
        return f"center ({cx}, {cy}, {cz})  size ({sx} x {sy} x {sz})  {self.volume:.0f} A^3"
