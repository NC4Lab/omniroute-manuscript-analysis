"""
Module: versioning.py

Purpose:
    Utilities for recording versioning metadata (git hash, processing timestamp)
    for all derived outputs and metadata objects in the pipeline.
"""

import subprocess
from datetime import datetime
import pickle
from pathlib import Path

from preprocessing.utils.omni_anal_logger import omni_anal_logger


def get_git_hash() -> str:
    """
    Get the current git commit hash of the codebase.

    Returns:
        str: Git SHA1 hash of current commit, or "unknown" if unavailable.
    """
    try:
        return subprocess.check_output(["git", "rev-parse", "HEAD"]).decode("utf-8").strip()
    except Exception:
        return "unknown"


def get_processed_timestamp() -> str:
    """
    Get the current timestamp in session-style format (e.g., '20250328_134136').

    Returns:
        str: Timestamp formatted as 'YYYYMMDD_HHMMSS'.
    """
    return datetime.now().strftime("%Y%m%d_%H%M%S")


def get_version_info() -> dict:
    """
    Get a dictionary containing version metadata.

    Returns:
        dict: Dictionary with keys 'git_hash' and 'processed_timestamp'.
    """
    return {
        "git_hash": get_git_hash(),
        "processed_timestamp": get_processed_timestamp(),
    }


def save_version_info(output_dir: Path, filename: str = "version_info.pkl") -> None:
    """
    Save version metadata to a file.

    Parameters:
        output_dir (Path): Directory to save the version info file.
        filename (str): Filename to use (default: 'version_info.pkl').
    """
    version_info = get_version_info()
    out_path = output_dir / filename
    with open(out_path, "wb") as f:
        pickle.dump(version_info, f)


def load_version_info(path: Path) -> dict:
    """
    Load version metadata from a pickle file.

    Parameters:
        path (Path): Full path to the version_info.pkl file.

    Returns:
        dict: Dictionary containing version metadata.
    """
    with open(path, "rb") as f:
        return pickle.load(f)


def print_version_info(version_info: dict) -> None:
    """
    Print git hash and processing timestamp in a readable format.

    Parameters:
        version_info (dict): Dictionary containing 'git_hash' and 'processed_timestamp'.
    """
    git_hash = version_info.get("git_hash", "unknown")
    timestamp = version_info.get("processed_timestamp", "unknown")
    
    omni_anal_logger.info(f"Version Info — Git hash: {git_hash} | Processed: {timestamp}")