"""Optional progress reporting for long-running pipeline steps.

Callers (for example a web app running HydroEO in a subprocess) register a
callback; pipeline code calls :func:`report`. With no callback registered,
:func:`report` does nothing, so library and CLI use is unaffected.

Examples
--------
>>> events = []
>>> with listen(events.append):
...     report("download", 1, 4, "tile")
>>> events[0].stage, events[0].done, events[0].total
('download', 1, 4)
"""

from __future__ import annotations

import contextlib
import logging
from collections.abc import Callable, Iterator
from dataclasses import dataclass

logger = logging.getLogger(__name__)


@dataclass(frozen=True)
class ProgressEvent:
    """One progress update: ``done`` of ``total`` units of ``stage``."""

    stage: str
    done: int
    total: int
    message: str = ""


ProgressCallback = Callable[[ProgressEvent], None]

_callback: ProgressCallback | None = None


def set_callback(callback: ProgressCallback | None) -> None:
    """Register ``callback`` for all later events, or ``None`` to remove it."""
    global _callback
    _callback = callback


@contextlib.contextmanager
def listen(callback: ProgressCallback) -> Iterator[None]:
    """Register ``callback`` inside a ``with`` block, restoring the previous one after."""
    global _callback
    previous, _callback = _callback, callback
    try:
        yield
    finally:
        _callback = previous


def report(stage: str, done: int, total: int, message: str = "") -> None:
    """Send a progress event to the registered callback, if any.

    A failing callback is logged and ignored: progress must never stop a run.
    """
    if _callback is None:
        return
    try:
        _callback(ProgressEvent(stage, done, total, message))
    except Exception:
        logger.exception("Progress callback failed")
