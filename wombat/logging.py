"""Structured logging setup for ConvergentDecidua."""

from __future__ import annotations

import logging

from rich.logging import RichHandler


def setup_logging(verbosity: int = 0) -> None:
    """Configure logging with Rich output.

    Parameters
    ----------
    verbosity : int
        0 = WARNING, 1 = INFO, 2+ = DEBUG.
    """
    level = logging.WARNING
    if verbosity == 1:
        level = logging.INFO
    elif verbosity >= 2:
        level = logging.DEBUG

    logging.basicConfig(
        level=level,
        format="%(message)s",
        datefmt="[%X]",
        handlers=[RichHandler(rich_tracebacks=True)],
    )
