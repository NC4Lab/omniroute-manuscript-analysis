"""
Module: io_trodes.py

Purpose:
    Utilities for reading continuous data from SpikeGadgets .rec files using SpikeInterface.
"""

from pathlib import Path
import os
import numpy as np
import subprocess
import spikeinterface.full as si
from spikeinterface.core import BaseRecording, NumpyRecording
from typing import Optional

from preprocessing.utils.omni_anal_logger import omni_anal_logger
from preprocessing.utils.versioning import save_version_info
from preprocessing.utils.binary_utils import TrodesDIOBinaryLoader
from preprocessing.utils.path import get_trodes_dir

def load_sample_rate_from_rec(rec_path: Path) -> None:
    """
    Load sampling rate from the .rec file and store it in the EphysMetadata.

    Parameters:
        rec_path (Path): Path to the .rec file.
    """
    rec = si.read_spikegadgets(rec_path)
    return rec.get_sampling_frequency()

def load_num_samples_from_rec(rec_path: Path) -> int:
    """
    Load total number of samples from the .rec file using SpikeInterface.

    Parameters:
        rec_path (Path): Path to the .rec file.

    Returns:
        int: Number of total samples in the recording.
    """
    rec = si.read_spikegadgets(rec_path)
    return rec.get_num_samples()

def load_gain_to_uV_from_rec(rec_path: Path) -> float:
    """
    Load the per-channel gain from the .rec file using SpikeInterface.

    Assumes that all channels have the same gain (which is typically true for SpikeGadgets CSC).
    Only the first entry of the 'gain_to_uV' array is returned.

    Parameters:
        rec_path (Path): Path to the .rec file.

    Returns:
        float: Gain value in microvolts per bit (e.g., 0.195).
    """
    rec = si.read_spikegadgets(rec_path)
    return float(rec.get_property("gain_to_uV")[0])

def get_dio_sg_ts(dio_df: np.ndarray, sampling_rate_hz: float) -> np.ndarray:
    """
    Compute SpikeGadgets timebase timestamps for each entry in the DIO trace.

    Parameters:
        dio_df (np.ndarray): DIO DataFrame with index as sample numbers (from load_dio_binary(...).dio).
        sampling_rate_hz (float): Sampling rate of the Trodes system in Hz.

    Returns:
        np.ndarray: Array of timestamps in SpikeGadgets timebase (in seconds).
    """
    return dio_df.index.to_numpy() / sampling_rate_hz


def get_csc_sg_ts(n_samples: int, sampling_rate_hz: float) -> np.ndarray:
    """
    Compute SpikeGadgets timebase timestamps for each sample in CSC data.

    Parameters:
        n_samples (int): Number of CSC samples.
        sampling_rate_hz (float): Sampling rate of the CSC recording in Hz.

    Returns:
        np.ndarray: Array of timestamps in SpikeGadgets timebase (in seconds).
    """
    return np.arange(n_samples) / sampling_rate_hz

def load_dio_binary(dio_dir: Path, channel: int) -> np.ndarray:
    """
    Load a binary DIO trace from a specific channel using TrodesDIOBinaryLoader.

    Parameters:
        dio_dir (Path): Path to the extracted DIO directory.
        channel (int): Channel number to load (e.g., 1 for Din1).

    Returns:
        np.ndarray: 1D array of 0/1 values.
    """
    # Look for a file ending in Din{channel}.dat
    pattern = f"Din{channel}.dat"

    matched_file = None
    for file in dio_dir.iterdir():
        if file.name.endswith(pattern):
            matched_file = file
            break

    if matched_file is None:
        raise FileNotFoundError(f"No Din{channel}.dat file found in {dio_dir}")

    # Load the .dat file and log the time taken
    with omni_anal_logger.time_block(f"Loading DIO channel {channel}"):
        return TrodesDIOBinaryLoader(matched_file)

def extract_dio_from_rec(rec_path: Path, dio_dir: Path, overwrite: bool = False) -> None:
    """
    Extract DIO channels from a .rec file using the SpikeGadgets exportdio tool.
    Saves them using SpikeGadgets .dat format.

    Each file is auto named by SpikeGadgets, which I fucking hate!
    (e.g., 20250328_134136.dio_Controller_Din1.dat).

    Parameters:
        rec_path (Path): Path to the .rec file.
        dio_dir (Path): Destination directory to save DIO files.
        overwrite (bool): If True, overwrite existing files. Default is False.
    """
    # Check if DIO directory already exists and overwrite is disabled
    if dio_dir.exists() and not overwrite:
        omni_anal_logger.warning(f"DIO already extracted at {dio_dir}")
        return

    # Resolve full path to the SpikeGadgets exportdio executable in the Trodes directory
    exportdio_exe_path = Path(get_trodes_dir()) / "exportdio.exe"
    if not exportdio_exe_path.exists():
        raise FileNotFoundError(f"exportdio.exe not found at {exportdio_exe_path}")

    # Build the exportdio command
    cmd = [
        str(exportdio_exe_path),
        "-rec", str(rec_path),
        "-outputdirectory", str(dio_dir.parent),
        "-output", dio_dir.stem,  # removes .DIO suffix
    ]

    # Run the exportdio tool as a subprocess
    with omni_anal_logger.time_block("Extracting DIO using exportdio"):
        result = subprocess.run(cmd, capture_output=True, text=True)

    # Log the subprocess output for inspection
    omni_anal_logger.info(f"exportdio stdout:\n{result.stdout}")
    if result.stderr:
        omni_anal_logger.warning(f"exportdio stderr:\n{result.stderr}")

    # Raise an error if exportdio failed
    if result.returncode != 0:
        raise RuntimeError(
            f"exportdio.exe failed with return code {result.returncode}.\n"
            f"Command: {' '.join(cmd)}\n"
            f"Stderr:\n{result.stderr}"
        )

    # Save version metadata alongside the extracted DIO files
    save_version_info(dio_dir)

def load_csc_from_rec(
    rec_path: Path,
    channel_trodes_id: list[int],
    channel_headstage_hardware_id: list[int],
    trodes_id_include: Optional[list[int]] = None
) -> BaseRecording:
    """
    Load CSC data from a .rec file using SpikeInterface, map hardware IDs to Trodes IDs,
    and optionally filter by a provided list of Trodes IDs.

    Parameters:
        rec_path (Path): Path to the .rec file.
        channel_trodes_id (list[int]): List of Trodes IDs (logical labels from CSV).
        channel_headstage_hardware_id (list[int]): Corresponding hardware IDs for slicing.
        trodes_id_include (Optional[list[int]]): Subset of Trodes IDs to include (if provided).

    Returns:
        BaseRecording: A mutable NumpyRecording with Trodes IDs and hardware ID properties.
    """
    if len(channel_trodes_id) != len(channel_headstage_hardware_id):
        raise ValueError("channel_trodes_id and channel_headstage_hardware_id must be the same length.")

    # Load original SpikeGadgets recording
    base_rec = si.read_spikegadgets(rec_path)

    # Raw extractor may return string channel IDs — convert to int for safe comparison
    raw_channel_ids = base_rec.get_channel_ids()
    id_map = {int(cid): cid for cid in raw_channel_ids}  # Maps int_id -> original (likely str)

    # Validate that all requested hardware IDs are available
    missing = [hid for hid in channel_headstage_hardware_id if hid not in id_map]
    if missing:
        omni_anal_logger.info(f"Available hardware IDs in .rec file: {sorted(id_map.keys())}")
        raise ValueError(f"The following hardware IDs are missing in .rec: {missing}")

    # Convert hardware IDs to the raw extractor format (likely strings) for slicing
    slice_ids = [id_map[hid] for hid in channel_headstage_hardware_id]
    base_rec = base_rec.channel_slice(channel_ids=slice_ids)

    # Convert to memory-backed NumpyRecording
    recording = NumpyRecording(
        traces_list=[base_rec.get_traces()],
        sampling_frequency=base_rec.get_sampling_frequency(),
        channel_ids=channel_trodes_id,  # Replace extractor-native IDs with Trodes IDs
    )
    recording.set_property("gain_to_uV", base_rec.get_property("gain_to_uV"))
    recording.set_property("offset_to_uV", base_rec.get_property("offset_to_uV"))
    recording.set_property("physical_unit", base_rec.get_property("physical_unit"))
    recording.set_property("gain_to_physical_unit", base_rec.get_property("gain_to_physical_unit"))
    recording.set_property("offset_to_physical_unit", base_rec.get_property("offset_to_physical_unit"))
    recording.set_property("channel_names", base_rec.get_property("channel_names"))
    recording.set_property("hardware_id", channel_headstage_hardware_id) # Store original hardware ID

    # Optionally filter to a subset of Trodes IDs
    if trodes_id_include is not None:
        missing = set(trodes_id_include) - set(channel_trodes_id)
        if missing:
            raise ValueError(f"Requested Trodes IDs not in known map: {missing}")
        recording = recording.channel_slice(channel_ids=trodes_id_include)

    return recording

def save_csc_binary(csc_dir: Path, recording: BaseRecording, overwrite: bool = False) -> None:
    """
    Save CSC trace data from a BaseRecording as separate raw .dat files
    (one per channel) in the specified CSC directory.

    Each file is named using the session base name + .csc_chan{N}.dat
    (e.g., 20250328_134136.csc_chan1.dat).

    Parameters:
        csc_dir (Path): Output directory (should end with .CSC).
        recording (BaseRecording): The pre-sliced recording object with channel_ids set to Trodes IDs.
        overwrite (bool): If False and file exists, skip saving that channel.
    """
    if not csc_dir.name.endswith(".CSC"):
        raise ValueError(f"CSC directory must end in .CSC, got: {csc_dir.name}")

    # Derive session base name from folder name (strip ".CSC")
    session_base_name = csc_dir.name.replace(".CSC", "")

    # Ensure output directory exists
    csc_dir.mkdir(parents=True, exist_ok=True)

    with omni_anal_logger.time_block("Saving CSC traces to individual .dat files"):
        traces = recording.get_traces()  # shape: (n_samples, n_channels)
        channel_ids = recording.get_channel_ids()

        for i, chan_id in enumerate(channel_ids):
            file_name = f"{session_base_name}.csc_chan{i+1}.dat"  # 1-based index
            file_path = csc_dir / file_name

            if file_path.exists() and not overwrite:
                omni_anal_logger.warning(f"File already exists at {file_path} — skipping.")
                continue

            traces[:, i].astype(np.float32).tofile(file_path)
            omni_anal_logger.info(f"Saved CSC channel {chan_id} to: {file_path}")

    # Save processing version info
    save_version_info(csc_dir)
    omni_anal_logger.info(f"Saved CSC version info to: {csc_dir}")

def load_csc_binary(
    csc_dir: Path,
    channel_trodes_id: list[int],
    sampling_rate_hz: float,
    trodes_id_include: Optional[list[int]] = None
) -> BaseRecording:
    """
    Load CSC traces from individual .dat files in a CSC directory and return as a BaseRecording.

    Parameters:
        csc_dir (Path): Directory containing per-channel .dat files.
        channel_trodes_id (list[int]): Trodes IDs corresponding to the channels in order.
        sampling_rate_hz (float): Sampling rate to embed in the BaseRecording.
        trodes_id_include (Optional[list[int]]): If provided, filters the channels to include only these Trodes IDs.

    Returns:
        NumpyRecording: A SpikeInterface BaseRecording object with channel_ids set to Trodes IDs.
    """
    if not csc_dir.exists() or not csc_dir.is_dir():
        raise FileNotFoundError(f"CSC directory not found: {csc_dir}")

    # Filter channels if inclusion list is provided
    if trodes_id_include is not None:
        id_mask = [tid in trodes_id_include for tid in channel_trodes_id]
        selected_ids = [tid for tid, keep in zip(channel_trodes_id, id_mask) if keep]
        file_indices = [i for i, keep in enumerate(id_mask) if keep]
    else:
        selected_ids = channel_trodes_id
        file_indices = list(range(len(channel_trodes_id)))

    # Load traces from each selected channel
    traces_list = []
    session_base = csc_dir.name.replace(".CSC", "")

    for i in file_indices:
        chan_idx = i + 1  # 1-based file naming convention
        file_path = csc_dir / f"{session_base}.csc_chan{chan_idx}.dat"

        if not file_path.exists():
            raise FileNotFoundError(f"Missing CSC channel file: {file_path}")

        trace = np.fromfile(file_path, dtype=np.float32)
        traces_list.append(trace)

    # Stack traces into shape (n_samples, n_channels)
    traces_array = np.stack(traces_list, axis=1)

    # Create and return BaseRecording
    recording = NumpyRecording(
        traces_list=[traces_array],
        sampling_frequency=sampling_rate_hz,
        channel_ids=selected_ids
    )
    return recording
