"""Shared test configuration.

Qt is forced onto the offscreen platform before PySide6 is imported anywhere, so
the GUI tests run identically on a developer's desktop and on a CI runner with no
display. Without this the suite passes locally and fails in CI, which is exactly
the failure mode Drydock is trying not to have.
"""

from __future__ import annotations

import os

os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")

import pytest  # noqa: E402


@pytest.fixture(scope="session")
def qapp():
    """One QApplication for the whole session.

    Qt does not support creating a second one in a process, so this is
    session-scoped rather than per-test.
    """
    from PySide6.QtWidgets import QApplication

    app = QApplication.instance() or QApplication([])
    yield app
