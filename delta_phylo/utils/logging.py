"""
Logging utilities for delta-phylo.

Provides a consistent logging configuration for the package.
"""

from __future__ import annotations

import logging
import sys
from typing import Optional


_FORMAT = "%(asctime)s [%(levelname)s] %(name)s: %(message)s"
_DATE_FORMAT = "%Y-%m-%d %H:%M:%S"


def get_logger(
    name: str,
    level: int = logging.INFO,
    handler: Optional[logging.Handler] = None,
) -> logging.Logger:
    """Get a configured logger for the given name.

    This function creates (or retrieves) a logger with the specified name,
    optionally attaching a handler.  When called from CLI entry points,
    a StreamHandler writing to ``stderr`` is added automatically.

    Args:
        name: Logger name (typically ``__name__`` of the calling module).
        level: Logging level (e.g. ``logging.DEBUG``).
        handler: Optional custom handler.  If None and the logger has no
            handlers, a :class:`logging.StreamHandler` writing to stderr
            is added.

    Returns:
        Configured :class:`logging.Logger`.
    """
    logger = logging.getLogger(name)
    logger.setLevel(level)

    if not logger.handlers:
        if handler is None:
            handler = logging.StreamHandler(sys.stderr)
            handler.setFormatter(logging.Formatter(_FORMAT, datefmt=_DATE_FORMAT))
        logger.addHandler(handler)

    return logger


def configure_root_logger(level: int = logging.INFO) -> None:
    """Configure the root logger for the entire delta-phylo package.

    Should be called once at application startup (e.g. in the CLI entry point).

    Args:
        level: Root logging level.
    """
    root = logging.getLogger("delta_phylo")
    root.setLevel(level)
    if not root.handlers:
        handler = logging.StreamHandler(sys.stderr)
        handler.setFormatter(logging.Formatter(_FORMAT, datefmt=_DATE_FORMAT))
        root.addHandler(handler)
