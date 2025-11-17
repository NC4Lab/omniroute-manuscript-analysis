"""
Module: ts_sync.py

Purpose:
    Compute time synchronization parameters between Trodes (SpikeGadgets) and ROS systems
    based on shared digital sync pulses.
"""

from pathlib import Path
import numpy as np
from pandas import DataFrame

from preprocessing.utils.omni_anal_logger import omni_anal_logger
from preprocessing.utils.io_trodes import load_dio_binary
from preprocessing.utils.io_rosbag import load_ros_sync_ts
from preprocessing.utils.path import get_synced_ts_dir
from preprocessing.utils.versioning import save_version_info

def compute_ts_sync_parameters(
    dio_path: Path,
    dio_channel: int,
    sampling_rate_hz: float,
    rosbag_path: Path
) -> dict[str, any]:
    """
    Compute polynomial sync mapping between SpikeGadgets and ROS timestamps.

    Parameters:
        dio_path (Path): Directory containing .dio files extracted from the .rec file.
        dio_channel (int): DIO channel to use for sync pulse detection.
        sampling_rate_hz (float): Sampling rate of the Trodes recording in Hz.
        rosbag_path (Path): Path to the ROS .bag file containing sync topic.

    Returns:
        dict: {
            "poly_coeffs": list of polynomial coefficients,
            "r_squared": float indicating fit quality
        }
    """
    # Load DIO sync pulses from Trodes (based on specified channel)
    dio_df = load_dio_binary(dio_path, channel=dio_channel).dio
    trodes_ts = dio_df[dio_df["state"] == True].index.to_numpy() / sampling_rate_hz
    omni_anal_logger.info(f"Found {len(trodes_ts)} Trodes sync pulses from Din{dio_channel}")

    # Load sync pulses from ROS bag
    ros_ts = load_ros_sync_ts(rosbag_path)
    omni_anal_logger.info(f"Found {len(ros_ts)} ROS sync pulses")

    if len(trodes_ts) < 2 or len(ros_ts) < 2:
        raise ValueError("Not enough sync pulses found in one or both data sources.")

    # Align timestamps using Needleman-Wunsch and compute polyfit
    p, x_match, y_match = align_timestamps_nw(trodes_ts, ros_ts, new=True)

    # Compute R² of fit
    predicted = np.polyval(p, x_match)
    ss_res = np.sum((y_match - predicted) ** 2)
    ss_tot = np.sum((y_match - np.mean(y_match)) ** 2)
    r_squared = 1 - ss_res / ss_tot if ss_tot > 0 else 0

    omni_anal_logger.info(f"Sync polyfit: y = {p[0]:.6f} * x + {p[1]:.6f}, R² = {r_squared:.5f}")

    # Validate alignment
    sync_check_passed = _validate_sync_alignment(
        trodes_ts=trodes_ts,
        ros_ts=ros_ts,
        dio_df=dio_df,
        sampling_rate_hz=sampling_rate_hz,
        x_match=x_match,
        y_match=y_match,
        poly_coeffs=np.array(p)
    )

    return {
        "poly_coeffs": p.tolist(),
        "r_squared": float(r_squared)
    }


def _validate_sync_alignment(
    trodes_ts: np.ndarray,
    ros_ts: np.ndarray,
    dio_df: DataFrame,
    sampling_rate_hz: float,
    x_match: np.ndarray,
    y_match: np.ndarray,
    poly_coeffs: np.ndarray
) -> bool:
    """
    Internal validation for timestamp sync.
    Checks whether mapped Trodes timestamps align with high-state DIO windows at ROS sync pulse times.
    
    Prints a summary of residuals, fit quality, and pulse alignment.
    """

    # Compute residuals on matched pairs
    predicted = np.polyval(poly_coeffs, x_match)
    residuals = y_match - predicted
    rms_error = np.sqrt(np.mean(residuals ** 2))
    max_residual = np.max(np.abs(residuals))
    mean_residual = np.mean(np.abs(residuals))
    slope = poly_coeffs[0]
    intercept = poly_coeffs[1]

    # Map SG timestamps to ROS timebase
    transformed_ts = np.polyval(poly_coeffs, trodes_ts)

    # For each ROS timestamp, find nearest transformed SG timestamp
    # Then map back to original SG time
    high_hits = []
    for ros_sync in ros_ts:
        idx = np.argmin(np.abs(transformed_ts - ros_sync))
        sg_ts = trodes_ts[idx]
        sg_sample = int(round(sg_ts * sampling_rate_hz))
        if sg_sample in dio_df.index:
            state = dio_df.loc[sg_sample, "state"]
            high_hits.append(state == 1)
        else:
            high_hits.append(False)

    # Compute high-state fraction
    high_hit_fraction = np.mean(high_hits)
    pass_check = high_hit_fraction > 0.9

    # Summary print
    omni_anal_logger.info("--- Sync Validation Summary ---")
    omni_anal_logger.info(f"  Matched pulse pairs        : {len(x_match)}")
    omni_anal_logger.info(f"  Total ROS sync pulses      : {len(ros_ts)}")
    omni_anal_logger.info(f"  RMS error (matched)        : {rms_error:.6f} s")
    omni_anal_logger.info(f"  Mean residual              : {mean_residual:.6f} s")
    omni_anal_logger.info(f"  Max residual               : {max_residual:.6f} s")
    omni_anal_logger.info(f"  Fit slope                  : {slope:.6f}")
    omni_anal_logger.info(f"  Fit intercept              : {intercept:.2f}")
    omni_anal_logger.info(f"  ROS syncs hit high DIO     : {sum(high_hits)} / {len(ros_ts)}")
    omni_anal_logger.info(f"  High-state hit fraction    : {high_hit_fraction:.3f}")
    omni_anal_logger.info(f"  PASS: {pass_check}")
    omni_anal_logger.info("--- End Sync Validation ---")

    # Raise error if validation fails
    if not pass_check:
        raise ValueError(
            f"Sync validation failed: high-state DIO hit fraction = {high_hit_fraction:.3f} (must be > 0.9)"
        )

def align_timestamps_nw(x, y, new, match=1, mismatch=1, gap=1, thresh=0.1):
    """
    Align sets of digital timestamps using the Needleman-Wunsch algorithm.

    Compares between timestamp difference vectors. If the differences go over
    a threshold, it is not a match. Constructs a scoring matrix and a direction
    matrix for the NW algorithm, and then computes a maximum score path
    through the scoring matrix. The path is decoded to form the optimal
    alignment. Good differences are used to form match vectors and a polynomial is
    fitted to determine drift of timestamps relative to each other.

    References: 
        http://en.wikipedia.org/wiki/Needleman%E2%80%93Wunsch_algorithm
        http://www.avatar.se/molbioinfo2001/dynprog/dynamic.html
        http://www.hrbc-genomics.net/training/bcd/Curric/PrwAli/node3.html

    Adapted from:
        Kamil Slowikowski (github: slowkow) Python NW implementation
        https://gist.github.com/slowkow/06c6dba9180d013dfd82bec217d22eb5

    and previously developed MATLAB version:
        Manu Madhav, 2015
        https://www.mathworks.com/matlabcentral/fileexchange/52819-align_timestamps

    Manu Madhav  
    25-Apr-2023
    """

    omni_anal_logger.info(f"Aligning timestamps of length {len(x)} and {len(y)}")

    # Compute difference vectors
    dx = np.diff(x)
    dy = np.diff(y)
    nx = len(dx)
    ny = len(dy)

    # Initialize scoring matrix and direction matrix
    F = np.zeros((nx + 1, ny + 1))
    F[:, 0] = np.linspace(0, -nx * gap, nx + 1)
    F[0, :] = np.linspace(0, -ny * gap, ny + 1)

    P = np.zeros((nx + 1, ny + 1))
    P[:, 0] = 3
    P[0, :] = 4

    # Fill scoring matrix based on similarity of differences
    t = np.zeros(3)
    for i in range(nx):
        for j in range(ny):
            if abs(dx[i] - dy[j]) <= thresh:
                t[0] = F[i, j] + match
            else:
                t[0] = F[i, j] - mismatch
            t[1] = F[i, j + 1] - gap
            t[2] = F[i + 1, j] - gap

            tmax = np.max(t)
            F[i + 1, j + 1] = tmax

            if t[0] == tmax:
                P[i + 1, j + 1] += 2
            if t[1] == tmax:
                P[i + 1, j + 1] += 3
            if t[2] == tmax:
                P[i + 1, j + 1] += 4

    # Trace back to find optimal alignment
    i = nx
    j = ny
    x_idx = np.array([])
    y_idx = np.array([])

    while i > 0 or j > 0:
        if P[i, j] in [2, 5, 6, 9]:
            x_idx = np.append(x_idx, i - 2)
            y_idx = np.append(y_idx, j - 2)
            i -= 1
            j -= 1
        elif P[i, j] in [3, 5, 7, 9]:
            x_idx = np.append(x_idx, i - 2)
            y_idx = np.append(y_idx, np.nan)
            i -= 1
        elif P[i, j] in [4, 6, 7, 9]:
            x_idx = np.append(x_idx, np.nan)
            y_idx = np.append(y_idx, j - 2)
            j -= 1

    # Reverse alignment indices
    x_idx = np.flip(x_idx)
    y_idx = np.flip(y_idx)

    # Identify valid matched differences
    good_diffs_idx = np.logical_not(np.logical_or(np.isnan(x_idx), np.isnan(y_idx)))

    x_match_idx = np.union1d(x_idx[good_diffs_idx], x_idx[good_diffs_idx] + 1).astype(int)[1:]
    y_match_idx = np.union1d(y_idx[good_diffs_idx], y_idx[good_diffs_idx] + 1).astype(int)[1:]

    x_match = x[x_match_idx]
    y_match = y[y_match_idx]

    omni_anal_logger.info(f"Found {len(x_match)} matching timestamps")

    # Fit polynomial
    if not new:
        # Fit residual polynomial: y_ts = np.polyval(p, x_ts) + x_ts
        p = np.polyfit(x_match, y_match - x_match, 1)
    else:
        # Fit direct polynomial: y_ts = np.polyval(p, x_ts)
        p = np.polyfit(x_match, y_match, 1)

    return p, x_match, y_match

def save_ts_sync_binary(ros_ts_array: np.ndarray, rat_id: str, session_name: str, overwrite: bool = False) -> None:
    """
    Save ROS-aligned timestamps for each SpikeGadgets samples to the session's synced timestamp directory.

    Parameters:
        ros_ts_array (np.ndarray): 1D array of ROS timestamps (one per CSC sample).
        rat_id (str): Unique animal ID (e.g., "NC40001").
        session_name (str): Session folder name (e.g., "20250328_134136").
        overwrite (bool): If False and file exists, skip saving.
    """

    ts_dir = get_synced_ts_dir(rat_id, session_name)
    ts_path = ts_dir / "ros_times_from_csc.dat"

    if ts_path.exists() and not overwrite:
        omni_anal_logger.warning(f"ROS timestamp .dat file already exists at {ts_path} — skipping save.")
        return

    ts_dir.mkdir(parents=True, exist_ok=True)

    # Save timestamp array as NumPy binary (with .dat suffix)
    with omni_anal_logger.time_block("Saving CSC-aligned ROS timestamps to .dat"):
        np.save(ts_path, ros_ts_array)

    # Save version info
    save_version_info(ts_dir)
    omni_anal_logger.info(f"Saved ROS-aligned timestamps to: {ts_path}")

def load_ts_sync_binary(ts_dir: Path) -> np.ndarray:
    """
    Load ROS-aligned CSC timestamp vector from a binary .dat file.

    Parameters:
        ts_dir (Path): Directory containing the timestamp file (e.g., processed/<session>.ts/)

    Returns:
        np.ndarray: 1D array of ROS timestamps (one per CSC sample).

    Raises:
        FileNotFoundError: If the expected .dat file is not found in the directory.
    """
    ts_path = ts_dir / "ros_times_from_csc.dat"
    if not ts_path.exists():
        raise FileNotFoundError(f"Timestamp binary not found at: {ts_path}")

    with omni_anal_logger.time_block("Loading synced ROS timestamps from binary"):
        return np.load(ts_path)
    
def convert_sg_ts_to_ros_time(
    sg_ts: np.ndarray,
    sync_mapping: dict[str, any]
) -> np.ndarray:
    """
    Convert SpikeGadgets timestamps to ROS time using fitted sync parameters.

    Parameters:
        sg_ts (np.ndarray): Array of timestamps in SpikeGadgets timebase (seconds).
        sync_mapping (dict): Contains 'poly_coeffs' from compute_ts_sync_parameters().

    Returns:
        np.ndarray: Array of timestamps in ROS timebase (seconds).
    """
    if "poly_coeffs" not in sync_mapping:
        raise ValueError("sync_mapping must contain 'poly_coeffs'")

    poly = np.array(sync_mapping["poly_coeffs"])
    ros_ts = np.polyval(poly, sg_ts)

    omni_anal_logger.info(f"Converted {len(sg_ts)} SG timestamps to ROS timebase")
    return ros_ts