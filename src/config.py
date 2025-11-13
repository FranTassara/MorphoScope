"""
Configuration File for MorphoScope - Structural Plasticity Analyzer

This module centralizes all application configuration constants and settings,
making it easy to modify behavior without changing core code.

Usage
-----
Import at the start of your main application:
    from config import Config, setup_logging
    setup_logging()  # Initialize logging system
    
Then access constants:
    print(Config.APP_NAME)
    print(Config.VERSION)

Author: Francisco Tassara
Date: 2025-11-12
Based on: Petsakou, Sapsis & Blau, Cell 2015
"""

import sys
import logging
from pathlib import Path


# =============================================================================
# WINDOWS UTF-8 CONFIGURATION (must be before logging setup)
# =============================================================================
if sys.platform == "win32":
    import codecs
    # Force UTF-8 encoding for console output on Windows
    # This prevents UnicodeEncodeError with special characters (✓, µ, etc.)
    try:
        sys.stdout = codecs.getwriter('utf-8')(sys.stdout.buffer, 'strict')
        sys.stderr = codecs.getwriter('utf-8')(sys.stderr.buffer, 'strict')
    except AttributeError:
        # Fallback for environments without buffer attribute
        pass


class Config:
    """
    Global configuration constants for MorphoScope application.
    
    This class provides a centralized location for all application settings,
    including UI colors, file paths, logging configuration, and default values.
    
    Attributes
    ----------
    APP_NAME : str
        Application display name
    VERSION : str
        Current version (semantic versioning: MAJOR.MINOR.PATCH)
    LOG_FILE : str
        Name of the log file
    LOG_LEVEL : int
        Logging level (DEBUG, INFO, WARNING, ERROR, CRITICAL)
    LOG_FORMAT : str
        Format string for log messages
    DEFAULT_VOXEL_SIZE_Z : float
        Default Z-axis voxel size in micrometers
    CONFIG_FILE : str
        User configuration file name (JSON format)
    COLOR_ROI_SELECTED : str
        Hex color for selected ROI visualization
    COLOR_PROCESSED : str
        Color name for successfully processed images
    COLOR_ERROR : str
        Hex color for processing errors
    CSV_HEADERS : list
        Column headers for CSV output file
    
    Notes
    -----
    - Colors support both hex (#RRGGBB) and named colors
    - Voxel sizes are in micrometers (µm)
    - Log level can be changed to logging.DEBUG for troubleshooting
    
    Examples
    --------
    >>> from config import Config
    >>> print(Config.APP_NAME)
    MorphoScope - Structural Plasticity Analyzer
    >>> print(Config.VERSION)
    2.0.0
    """
    
    # =========================================================================
    # APPLICATION INFORMATION
    # =========================================================================
    APP_NAME = "MorphoScope - Structural Plasticity Analyzer"
    VERSION = "2.0.0"
    AUTHOR = "Francisco Tassara"
    INSTITUTION = "Your Institution"  # TODO: Update with your institution
    
    # =========================================================================
    # LOGGING CONFIGURATION
    # =========================================================================
    LOG_FILE = "morphoscope.log"
    LOG_LEVEL = logging.INFO  # Change to logging.DEBUG for detailed output
    LOG_FORMAT = '%(asctime)s - %(name)s - %(levelname)s - %(message)s'
    
    # Console logging format (simpler for terminal display)
    CONSOLE_LOG_FORMAT = '%(levelname)s: %(message)s'
    
    # =========================================================================
    # IMAGE PROCESSING DEFAULTS
    # =========================================================================
    # These are fallback values when metadata is unavailable
    DEFAULT_VOXEL_SIZE_X = 0.1  # µm (typically for high-res confocal)
    DEFAULT_VOXEL_SIZE_Y = 0.1  # µm
    DEFAULT_VOXEL_SIZE_Z = 1.0  # µm (typical Z-step size)
    
    # Threshold for surface contour identification (percentage of max intensity)
    SURFACE_THRESHOLD = 0.7  # 70% of maximum fluorescence
    
    # =========================================================================
    # FILE PATHS AND NAMES
    # =========================================================================
    CONFIG_FILE = "user_config.json"  # User preferences (saved between sessions)
    DEFAULT_OUTPUT_NAME = "results.csv"  # Default CSV filename
    
    # =========================================================================
    # UI COLORS AND STYLING
    # =========================================================================
    # ROI visualization
    COLOR_ROI_SELECTED = "#90EE90"  # Light green for selected ROI
    COLOR_ROI_ACTIVE = "lightgreen"  # ROI creation mode indicator
    
    # Image list feedback
    COLOR_PROCESSED = "lightgreen"  # Successfully processed image
    COLOR_ERROR = "#FB889A"  # Error during processing (light red)
    COLOR_PENDING = "white"  # Not yet processed
    
    # =========================================================================
    # CSV EXPORT CONFIGURATION
    # =========================================================================
    CSV_HEADERS = [
        'Image filename',
        'Spread x [pixel]',
        'Spread y [pixel]',
        'Spread z [pixel]',
        'Spread x*y [pixel²]',
        'Spread x*y*z [pixel³]',
        'Spread x [µm]',
        'Spread y [µm]',
        'Spread z [µm]',
        'Spread x*y [µm²]',
        'Spread x*y*z [µm³]',
        'Axonal Volume',
        'Fluorescence_px',
        'Fluorescence_um',
        'Observation'
    ]
    
    # =========================================================================
    # SUPPORTED FILE FORMATS
    # =========================================================================
    SUPPORTED_FORMATS = {
        'czi': 'Zeiss CZI files',
        'tif': 'TIFF stacks',
        'tiff': 'TIFF stacks',
        'lsm': 'Zeiss LSM files'
    }
    
    FILE_DIALOG_FILTER = "Image Files (*.tif *.tiff *.czi *.lsm)"
    
    # =========================================================================
    # PROCESSING PARAMETERS
    # =========================================================================
    # Maximum values for validation
    MAX_VOXEL_SIZE_XY = 10.0  # µm (larger values trigger warning)
    MAX_VOXEL_SIZE_Z = 50.0   # µm
    
    # Interpolation order for image rotation
    ROTATION_INTERPOLATION_ORDER = 1  # Bilinear (0=nearest, 1=bilinear, 3=cubic)


def setup_logging(log_to_file: bool = True, verbose: bool = False):
    """
    Configure logging system for the application.
    
    This function initializes the logging system with appropriate handlers
    and formatters. It should be called once at application startup.
    
    Parameters
    ----------
    log_to_file : bool, optional
        If True, log messages are written to a file (default: True)
    verbose : bool, optional
        If True, use DEBUG level for detailed output (default: False)
    
    Notes
    -----
    - Automatically configures UTF-8 encoding for Windows
    - Creates log file in the current working directory
    - Reduces noise from third-party libraries (matplotlib, PIL)
    - Uses different formats for file and console output
    
    Examples
    --------
    Basic usage:
    >>> setup_logging()
    
    Verbose mode for debugging:
    >>> setup_logging(verbose=True)
    
    Console-only logging:
    >>> setup_logging(log_to_file=False)
    """
    # Set logging level
    level = logging.DEBUG if verbose else Config.LOG_LEVEL
    
    # Prepare handlers
    handlers = []
    
    # File handler (if enabled)
    if log_to_file:
        file_handler = logging.FileHandler(
            Config.LOG_FILE, 
            encoding='utf-8'  # Explicit UTF-8 for log file
        )
        file_handler.setFormatter(logging.Formatter(Config.LOG_FORMAT))
        handlers.append(file_handler)
    
    # Console handler (always enabled)
    console_handler = logging.StreamHandler(sys.stdout)
    console_handler.setFormatter(logging.Formatter(Config.LOG_FORMAT))
    handlers.append(console_handler)
    
    # Configure root logger
    logging.basicConfig(
        level=level,
        handlers=handlers
    )
    
    # Reduce noise from third-party libraries
    logging.getLogger('matplotlib').setLevel(logging.WARNING)
    logging.getLogger('PIL').setLevel(logging.WARNING)
    logging.getLogger('tifffile').setLevel(logging.WARNING)
    
    # Log successful initialization
    logger = logging.getLogger(__name__)
    logger.info(f"{Config.APP_NAME} v{Config.VERSION}")
    logger.info(f"Logging initialized - Level: {logging.getLevelName(level)}")
    if log_to_file:
        logger.info(f"Log file: {Config.LOG_FILE}")


def get_version_info() -> dict:
    """
    Get application version information.
    
    Returns
    -------
    dict
        Dictionary containing version details:
        - 'version': Version string
        - 'app_name': Application name
        - 'author': Author name
        - 'institution': Institution name
    
    Examples
    --------
    >>> info = get_version_info()
    >>> print(f"{info['app_name']} v{info['version']}")
    """
    return {
        'version': Config.VERSION,
        'app_name': Config.APP_NAME,
        'author': Config.AUTHOR,
        'institution': Config.INSTITUTION
    }


def validate_voxel_sizes(voxel_x: float, voxel_y: float, voxel_z: float) -> tuple:
    """
    Validate voxel size parameters against configuration limits.
    
    Parameters
    ----------
    voxel_x, voxel_y, voxel_z : float
        Voxel dimensions in micrometers
    
    Returns
    -------
    tuple
        (is_valid: bool, warning_message: str)
        is_valid is True if all values are within acceptable ranges
        warning_message is empty string if valid, otherwise contains warning
    
    Examples
    --------
    >>> is_valid, msg = validate_voxel_sizes(0.1, 0.1, 1.0)
    >>> print(is_valid)  # True
    >>> is_valid, msg = validate_voxel_sizes(0.1, 0.1, 100.0)
    >>> print(msg)  # "Z voxel size exceeds recommended maximum..."
    """
    if voxel_x <= 0 or voxel_y <= 0 or voxel_z <= 0:
        return False, "All voxel sizes must be positive numbers"
    
    warnings = []
    
    if voxel_x > Config.MAX_VOXEL_SIZE_XY or voxel_y > Config.MAX_VOXEL_SIZE_XY:
        warnings.append(
            f"XY voxel size exceeds {Config.MAX_VOXEL_SIZE_XY}µm "
            f"(unusual for confocal microscopy)"
        )
    
    if voxel_z > Config.MAX_VOXEL_SIZE_Z:
        warnings.append(
            f"Z voxel size exceeds {Config.MAX_VOXEL_SIZE_Z}µm "
            f"(very large Z-step)"
        )
    
    if warnings:
        return False, " | ".join(warnings)
    
    return True, ""


# =============================================================================
# MODULE INITIALIZATION
# =============================================================================

# This runs when the module is imported
logger = logging.getLogger(__name__)

# Inform user if running on Windows (UTF-8 configuration applied)
if sys.platform == "win32":
    logger.debug("Windows platform detected - UTF-8 console encoding configured")