"""Logging configuration for HydroEO application."""

import logging


def setup_logging(level=logging.INFO):
    """Configure logging for HydroEO to display messages in terminal.

    Configures the HydroEO logger to output to console with a standard format.
    This should be called at application startup before running any HydroEO operations.

    Parameters
    ----------
    level : int
        Logging level (logging.DEBUG, logging.INFO, logging.WARNING, etc.)
        Default is logging.INFO.

    Examples
    --------
    >>> from HydroEO.logging_config import setup_logging
    >>> setup_logging(logging.DEBUG)  # Enable debug messages
    """
    logger = logging.getLogger("HydroEO")
    logger.setLevel(level)

    # Create console handler if one doesn't exist
    if not logger.handlers:
        handler = logging.StreamHandler()
        handler.setLevel(level)

        # Create formatter
        formatter = logging.Formatter(
            "%(asctime)s - %(name)s - %(levelname)s - %(message)s",
            datefmt="%Y-%m-%d %H:%M:%S",
        )
        handler.setFormatter(formatter)
        logger.addHandler(handler)
