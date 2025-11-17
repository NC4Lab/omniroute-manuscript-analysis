"""
Module: path.py

Purpose:
    Standardized utilities to resolve all key file and folder paths used in the project.
"""

from pathlib import Path

from preprocessing.utils.config import DATASET_ROOT, TRODES_DIR

# ------------------------ #
# Tools / Environment Directories
# ------------------------ #

def get_trodes_dir() -> Path:
    """
    Returns: 'C:/Program Files/Trodes/'
    """
    return Path(TRODES_DIR)

# ------------------------ #
# Data Directories
# ------------------------ #

def get_data_dir() -> Path:
    """
    Returns: 'data/'
    """
    return Path(DATASET_ROOT)

def get_rat_dir(rat_id: str) -> Path:
    """
    Parameters:
        rat_id (str)

    Returns: 'data/NC40001'
    """
    return get_data_dir() / rat_id

def get_raw_dir(rat_id: str, session_name: str) -> Path:
    """
    Parameters:
        rat_id (str)
        session_name (str)

    Returns: 'data/NC40001/20250328_134136/raw'
    """
    return get_session_path(rat_id, session_name) / "raw"

def get_processed_dir(rat_id: str, session_name: str) -> Path:
    """
    Parameters:
        rat_id (str)
        session_name (str)

    Returns: 'data/NC40001/20250328_134136/processed'
    """
    return get_session_path(rat_id, session_name) / "processed"

def get_dio_dir(rat_id: str, session_name: str) -> Path:
    """
    Parameters:
        rat_id (str)
        session_name (str)

    Returns: 'data/NC40001/20250328_134136/processed/20250328_134136.DIO'
    """
    return get_processed_dir(rat_id, session_name) / f"{session_name}.DIO"

def get_csc_dir(rat_id: str, session_name: str) -> Path:
    """
    Parameters:
        rat_id (str)
        session_name (str)

    Returns: 'data/NC40001/20250328_134136/processed/20250328_134136.CSC'
    """
    return get_processed_dir(rat_id, session_name) / f"{session_name}.CSC"

def get_spike_dir(rat_id: str, session_name: str) -> Path:
    """
    Parameters:
        rat_id (str)
        session_name (str)

    Returns: 'data/NC40001/20250328_134136/processed/20250328_134136.spikes'
    """
    return get_processed_dir(rat_id, session_name) / f"{session_name}.spikes"

def get_synced_ts_dir(rat_id: str, session_name: str) -> Path:
    """
    Parameters:
        rat_id (str)
        session_name (str)

    Returns: 'data/NC40001/20250328_134136/processed/20250328_134136.ts'
    """
    return get_processed_dir(rat_id, session_name) / f"{session_name}.ts"


# ------------------------ #
# Data Paths
# ------------------------ #

def get_experiment_metadata_csv_path() -> Path:
    """
    Returns: 'data/experiment_metadata.csv'
    """
    return get_data_dir() / "experiment_metadata.csv"

def get_ephys_channel_map_csv_path(rat_id: str) -> Path:
    """
    Parameters:
        rat_id (str)

    Returns: 'data/NC40001/ephys_channel_map_metadata.csv'
    """
    return get_rat_dir(rat_id) / "ephys_channel_map_metadata.csv"

def get_session_path(rat_id: str, session_name: str) -> Path:
    """
    Parameters:
        rat_id (str)
        session_name (str)

    Returns: 'data/NC40001/20250328_134136'
    """
    return get_rat_dir(rat_id) / session_name

def get_rec_path(rat_id: str, session_name: str) -> Path:
    """
    Parameters:
        rat_id (str)
        session_name (str)

    Returns: 'data/NC40001/20250328_134136/raw/Trodes/20250328_134136_merged.rec'
    """
    return get_raw_dir(rat_id, session_name) / "Trodes" / f"{session_name}_merged.rec"

def get_rosbag_path(rat_id: str, session_name: str) -> Path:
    """
    Parameters:
        rat_id (str)
        session_name (str)

    Returns: 'data/NC40001/20250328_134136/raw/ROS/ros_session_*.bag'
    """
    ros_dir = get_raw_dir(rat_id, session_name) / "ROS"
    bag_files = list(ros_dir.glob("*.bag"))
    if not bag_files:
        raise FileNotFoundError(f"No .bag file found in: {ros_dir}")
    return bag_files[0]

def get_session_metadata_path(rat_id: str, session_name: str) -> Path:
    """
    Parameters:
        rat_id (str)
        session_name (str)

    Returns: 'data/NC40001/20250328_134136/processed/session_metadata.pkl'
    """
    return get_processed_dir(rat_id, session_name) / "session_metadata.pkl"

def get_ephys_metadata_path(rat_id: str, session_name: str) -> Path:
    """
    Parameters:
        rat_id (str)
        session_name (str)

    Returns: 'data/NC40001/20250328_134136/processed/ephys_metadata.pkl'
    """
    return get_processed_dir(rat_id, session_name) / "ephys_metadata.pkl"
