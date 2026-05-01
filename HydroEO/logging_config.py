"""Logging configuration for HydroEO application."""

import logging
from datetime import datetime
from pathlib import Path


def _ensure_log_dir():
    """Create logs directory at project root if it doesn't exist.

    Returns
    -------
    Path
        Path to the logs directory.
    """
    project_root = Path(
        __file__
    ).parent.parent  # HydroEO/logging_config.py -> HydroEO/ -> project root
    logs_dir = project_root / "logs"
    logs_dir.mkdir(exist_ok=True)
    return logs_dir


def _get_log_filepath():
    """Generate a timestamped log file path.

    Returns
    -------
    Path
        Path to the log file with timestamp.
    """
    logs_dir = _ensure_log_dir()
    timestamp = datetime.now().strftime("%Y-%m-%d_%H-%M-%S")
    return logs_dir / f"HydroEO_{timestamp}.log"


def setup_logging(level=logging.INFO, enable_file_logging=True):
    """Configure logging for HydroEO with console and optional file output.

    Configures the HydroEO logger to output to console at the specified level.
    If file logging is enabled, also logs to a timestamped file at DEBUG level.
    This should be called at application startup before running any HydroEO operations.

    Parameters
    ----------
    level : int
        Logging level for console output (logging.DEBUG, logging.INFO, logging.WARNING, etc.)
        Default is logging.INFO.
    enable_file_logging : bool
        If True, also logs to a timestamped file at DEBUG level.
        Default is True.

    Examples
    --------
    >>> from HydroEO.logging_config import setup_logging
    >>> setup_logging(logging.DEBUG)  # Console at DEBUG, file at DEBUG
    >>> setup_logging(logging.INFO, enable_file_logging=False)  # Console only
    """
    logger = logging.getLogger("HydroEO")
    logger.setLevel(logging.DEBUG)  # Logger itself should be at DEBUG to capture all

    # Create formatter
    formatter = logging.Formatter(
        "%(asctime)s - %(name)s - %(levelname)s - %(message)s",
        datefmt="%Y-%m-%d %H:%M:%S",
    )

    # Create console handler if one doesn't exist
    if not logger.handlers:
        console_handler = logging.StreamHandler()
        console_handler.setLevel(level)
        console_handler.setFormatter(formatter)
        logger.addHandler(console_handler)

        # Create file handler if enabled
        if enable_file_logging:
            log_filepath = _get_log_filepath()
            file_handler = logging.FileHandler(log_filepath)
            file_handler.setLevel(logging.DEBUG)
            file_handler.setFormatter(formatter)
            logger.addHandler(file_handler)
