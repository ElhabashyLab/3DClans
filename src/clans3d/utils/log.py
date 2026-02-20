"""
Centralized logging configuration for Clans-3D.

Call ``setup_logging()`` once at application start-up (e.g. in ``main()``)
before any other clans3d imports that emit log messages.

Every module should create its own logger via::

    import logging
    logger = logging.getLogger(__name__)
"""

import logging


def setup_logging(verbose: bool = False) -> None:
    """Configure the *clans3d* logger hierarchy.

    Args:
        verbose: If ``True`` the log level is set to ``DEBUG``,
                 otherwise ``INFO``.
    """
    level = logging.DEBUG if verbose else logging.INFO
    formatter = logging.Formatter(
        fmt="[%(levelname)-7s] %(message)s",
    )
    handler = logging.StreamHandler()
    handler.setFormatter(formatter)

    root = logging.getLogger("clans3d")
    # Avoid duplicate handlers when setup_logging is called more than once
    if not root.handlers:
        root.addHandler(handler)
    else:
        root.handlers[0].setFormatter(formatter)
    root.setLevel(level)
    root.propagate = False
