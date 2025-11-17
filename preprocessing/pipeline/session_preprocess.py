"""
Module: session_preprocess.py

Purpose:
    This module performs initial preprocessing steps required to prepare a session for analysis.
    These steps establish the core context—file paths, session type, channel
    metadata, sampling rate, and timebase mapping—that are required for all downstream processing
    and analysis stages.
"""

from preprocessing.utils.metadata import SessionMetadata, EphysMetadata
from preprocessing.utils.path import (
    get_rec_path,
    get_dio_dir,
    get_rosbag_path,
)
from preprocessing.utils.io_trodes import extract_dio_from_rec, get_csc_sg_ts
from preprocessing.utils.ts_sync import compute_ts_sync_parameters, convert_sg_ts_to_ros_time, save_ts_sync_binary
from preprocessing.utils.omni_anal_logger import omni_anal_logger


def setup_session_metadata(rat_id: str, session_name: str, overwrite: bool = False) -> None:
    """
    Initialize and save SessionMetadata and, if applicable, EphysMetadata for a given session.

    Parameters:
        rat_id (str): Unique animal ID (e.g., "NC40001")
        session_name (str): Session folder name (e.g., "20250328_134136")
        overwrite (bool): If True, overwrite existing metadata pickle files

    Returns:
        None
    """
    omni_anal_logger.info(f"--- Setting up metadata for session: {rat_id}/{session_name} ---")

    # Create and save SessionMetadata object
    session_metadata = SessionMetadata(rat_id, session_name)
    session_metadata.load_or_initialize_pickle(overwrite=overwrite)
    session_metadata.save_pickle(overwrite=overwrite)
    omni_anal_logger.info("SessionMetadata initialized and saved")

    # If session is an ephys recording, create and save EphysMetadata
    if session_metadata.session_type == "ephys":
        ephys_metadata = EphysMetadata(rat_id, session_name)
        ephys_metadata.load_or_initialize_pickle(overwrite=overwrite)
        ephys_metadata.save_pickle(overwrite=overwrite)
        omni_anal_logger.info("EphysMetadata initialized and saved")
    else:
        omni_anal_logger.info("Session type is not 'ephys' — skipping EphysMetadata setup")

    omni_anal_logger.info(f"--- Metadata setup complete for session: {rat_id}/{session_name} ---")


def setup_dio_extraction(rat_id: str, session_name: str, overwrite: bool = False) -> None:
    """
    Extract DIO event channel files from a .rec file using SpikeGadgets exportdio.

    Parameters:
        rat_id (str): Unique animal ID
        session_name (str): Session folder name
        overwrite (bool): If True, force re-extraction even if DIO folder exists

    Returns:
        None
    """
    dio_dir = get_dio_dir(rat_id, session_name)
    rec_path = get_rec_path(rat_id, session_name)

    # Skip extraction if already done and overwrite is False
    if dio_dir.exists() and not overwrite:
        omni_anal_logger.warning(f"DIO directory already exists at {dio_dir} — skipping extraction.")
        return

    # Perform DIO extraction
    extract_dio_from_rec(rec_path=rec_path, dio_dir=dio_dir, overwrite=overwrite)
    omni_anal_logger.info("DIO extraction complete.")


def setup_timestamp_sync(rat_id: str, session_name: str, dio_channel: int = 2, overwrite: bool = False, save_ts_binary: bool = False) -> None:
    """
    Compute and save timestamp synchronization parameters between SpikeGadgets and ROS timebases.

    Parameters:
        rat_id (str): Unique animal ID
        session_name (str): Session folder name
        dio_channel (int): DIO channel to use for sync pulse detection (default = 2)
        overwrite (bool): If True, overwrite existing timestamp_mapping in metadata
        save_ts_binary (bool): If True, generate and save ROS-aligned timestamps for each CSC sample

    Returns:
        None
    """

    omni_anal_logger.info(f"--- Starting timestamp sync for session: {rat_id}/{session_name} ---")

    # Load session-level metadata and check session type
    session_metadata = SessionMetadata(rat_id, session_name)
    session_metadata.load_or_initialize_pickle(overwrite=False) # Do not overwrite existing metadata

    if session_metadata.session_type != "ephys":
        omni_anal_logger.info("Session type is not 'ephys' — skipping timestamp synchronization.")
        return

    # Load EphysMetadata to check for existing sync and sampling rate
    ephys_metadata = EphysMetadata(rat_id, session_name)
    ephys_metadata.load_or_initialize_pickle(overwrite=False) # Do not overwrite existing metadata

    if ephys_metadata.timestamp_mapping is not None and not overwrite:
        omni_anal_logger.warning("Timestamp sync already computed — skipping (use overwrite=True to force).")
        return

    # Prepare paths to DIO, .rec, and .bag files
    dio_dir = get_dio_dir(rat_id, session_name)
    rec_path = get_rec_path(rat_id, session_name)
    rosbag_path = get_rosbag_path(rat_id, session_name)

    # Ensure DIO extraction is present
    extract_dio_from_rec(rec_path=rec_path, dio_dir=dio_dir, overwrite=overwrite)

    # Compute and store sync parameters
    sync_mapping = compute_ts_sync_parameters(
        dio_path=dio_dir,
        dio_channel=dio_channel,
        sampling_rate_hz=ephys_metadata.sampling_rate_hz,
        rosbag_path=rosbag_path,
    )

    # Store sync parameters in EphysMetadata and save
    ephys_metadata.timestamp_mapping = sync_mapping
    ephys_metadata.save_pickle(overwrite=True) # Force overwrite to ensure sync parameters are saved

    omni_anal_logger.info("Timestamp sync parameters saved to EphysMetadata")
    omni_anal_logger.info(f"Sync polyfit: {sync_mapping['poly_coeffs']}, R² = {sync_mapping['r_squared']:.5f}")

    # Optionally compute and save full synced CSC timestamp vector
    if save_ts_binary:
        sg_ts = get_csc_sg_ts(ephys_metadata.num_samples, ephys_metadata.sampling_rate_hz)
        ros_ts = convert_sg_ts_to_ros_time(sg_ts, sync_mapping)
        save_ts_sync_binary(ros_ts, rat_id, session_name, overwrite=overwrite)

    omni_anal_logger.info(f"--- Timestamp sync complete for session: {rat_id}/{session_name} ---")

