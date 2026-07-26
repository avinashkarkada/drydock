"""Table model for screening results.

The reason this is a ``QAbstractTableModel`` rather than a widget-per-row table:
a screen produces tens of thousands of rows, and Qt only ever asks the model for
the cells it is about to paint. Rows exist as plain Python objects and become
"widgets" only while visible. A 47,000-row `QTableWidget` would allocate 47,000
rows of widgets up front and make the interface unusable long before the run
finished.

The model is append-oriented, matching how results arrive: the watcher hands over
whatever appeared in the journal since the last poll, and the model appends.
"""

from __future__ import annotations

from typing import Any

from PySide6.QtCore import QAbstractTableModel, QModelIndex, Qt

from drydock.core.rundir import LigandResult

# Qt's model API requires an invalid QModelIndex as the default parent, but
# building one per call in a default argument is evaluated once at import and
# shared anyway. Naming that shared instance makes the intent explicit.
_NO_PARENT = QModelIndex()


class ResultsTableModel(QAbstractTableModel):
    """Ligand results, sortable, appended to as a run progresses."""

    COLUMNS: tuple[tuple[str, str], ...] = (
        ("#", "rank"),
        ("Ligand", "ligand_id"),
        ("Affinity", "best_affinity"),
        ("Modes", "n_modes"),
        ("Status", "status"),
        ("Time (s)", "elapsed_s"),
        ("Seed", "seed"),
    )

    def __init__(self, parent=None) -> None:
        super().__init__(parent)
        self._rows: list[LigandResult] = []
        self._sort_column: int = 2  # affinity
        self._sort_order: Qt.SortOrder = Qt.SortOrder.AscendingOrder

    # -- Qt interface ------------------------------------------------------

    def rowCount(self, parent: QModelIndex = _NO_PARENT) -> int:
        return 0 if parent.isValid() else len(self._rows)

    def columnCount(self, parent: QModelIndex = _NO_PARENT) -> int:
        return 0 if parent.isValid() else len(self.COLUMNS)

    def headerData(self, section: int, orientation: Qt.Orientation, role: int = 0) -> Any:
        if role != Qt.ItemDataRole.DisplayRole:
            return None
        if orientation == Qt.Orientation.Horizontal:
            return self.COLUMNS[section][0]
        return section + 1

    def data(self, index: QModelIndex, role: int = Qt.ItemDataRole.DisplayRole) -> Any:
        if not index.isValid():
            return None
        record = self._rows[index.row()]
        field = self.COLUMNS[index.column()][1]

        if role == Qt.ItemDataRole.DisplayRole:
            return self._display(record, field, index.row())

        if role == Qt.ItemDataRole.TextAlignmentRole:
            if field in ("rank", "best_affinity", "n_modes", "elapsed_s", "seed"):
                return int(Qt.AlignmentFlag.AlignRight | Qt.AlignmentFlag.AlignVCenter)
            return int(Qt.AlignmentFlag.AlignLeft | Qt.AlignmentFlag.AlignVCenter)

        if role == Qt.ItemDataRole.ToolTipRole and record.error:
            return record.error

        # Sorting a mixed-type column through the display string would order
        # "-10.2" before "-9.1". Expose the raw value for anything that needs to
        # compare rows.
        if role == Qt.ItemDataRole.UserRole:
            return self._raw(record, field, index.row())

        return None

    def sort(self, column: int, order: Qt.SortOrder = Qt.SortOrder.AscendingOrder) -> None:
        field = self.COLUMNS[column][1]
        if field == "rank":
            # Rank is positional, so sorting by it means restoring affinity order.
            field = "best_affinity"

        self.layoutAboutToBeChanged.emit()
        reverse = order == Qt.SortOrder.DescendingOrder
        self._rows.sort(key=lambda r: _sort_key(r, field), reverse=reverse)
        self._sort_column, self._sort_order = column, order
        self.layoutChanged.emit()

    # -- population --------------------------------------------------------

    def append(self, records: list[LigandResult]) -> None:
        """Add records, keeping the current sort.

        Appending then re-sorting is O(n log n) per poll, which is fine at this
        scale and avoids the bookkeeping of an insertion sort. If it ever stops
        being fine, this is where to look.
        """
        if not records:
            return
        first = len(self._rows)
        self.beginInsertRows(_NO_PARENT, first, first + len(records) - 1)
        self._rows.extend(records)
        self.endInsertRows()
        self.sort(self._sort_column, self._sort_order)

    def reset(self, records: list[LigandResult] | None = None) -> None:
        """Replace all rows. Used when attaching to a different run."""
        self.beginResetModel()
        self._rows = list(records or [])
        self.endResetModel()
        if self._rows:
            self.sort(self._sort_column, self._sort_order)

    def record_at(self, row: int) -> LigandResult | None:
        if 0 <= row < len(self._rows):
            return self._rows[row]
        return None

    @property
    def records(self) -> list[LigandResult]:
        return self._rows

    # -- formatting --------------------------------------------------------

    def _display(self, record: LigandResult, field: str, row: int) -> str:
        if field == "rank":
            return str(row + 1)
        if field == "ligand_id":
            return record.ligand_id
        if field == "best_affinity":
            affinity = record.best_affinity
            return "-" if affinity is None else f"{affinity:.1f}"
        if field == "n_modes":
            return str(len(record.modes))
        if field == "status":
            return record.status
        if field == "elapsed_s":
            return f"{record.elapsed_s:.1f}"
        if field == "seed":
            return "" if record.seed is None else str(record.seed)
        return ""

    def _raw(self, record: LigandResult, field: str, row: int) -> Any:
        if field == "rank":
            return row
        if field == "n_modes":
            return len(record.modes)
        return getattr(record, field, None)


def _sort_key(record: LigandResult, field: str) -> tuple[int, Any]:
    """Sort key that keeps missing values together at the end.

    Failed ligands have no affinity. Sorting them as ``inf`` would be wrong for a
    descending sort (they would lead), so absence is encoded in a leading flag
    that always sorts last in ascending order and first in descending -- either
    way, grouped and out of the way of real results.
    """
    if field == "n_modes":
        return (0, len(record.modes))

    value = getattr(record, field, None)
    if value is None:
        return (1, 0)
    return (0, value)
