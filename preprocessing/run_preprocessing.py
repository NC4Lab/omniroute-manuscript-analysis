"""
Script: run_preprocessing.py

Purpose:
    Initialize and load experiment-level metadata for batch session preprocessing.
    Export session CSC data to MATLAB .mat files for downstream analysis.

Run script:
    python -m preprocessing.run_preprocessing
"""

import sys
import subprocess
from preprocessing.utils.metadata import ExperimentMetadata
from preprocessing.utils.omni_anal_logger import omni_anal_logger
from preprocessing.pipeline.session_preprocess import (
    setup_session_metadata,
    setup_dio_extraction,
    setup_timestamp_sync,
)

#-----------------------------------------------
# Parameters 
#-----------------------------------------------

# Flag to overwrite existing metadata files
experiment_csv_overwrite = False

# Flag to overwrite existing session/ephys metadata files
metadata_overwrite = True

#Flag to overwrite existing DIO files
dio_overwrite = True

# DIO channel index used for sync pulses (as per your rig conventions)
dio_channel = 1

# Flag to overwrite existing timestamp sync files
ts_sync_overwrite = True

# Flag to save timestamp sync binary files
save_ts_sync_binary = True

# Flag to overwrite existing MATLAB CSC export
matlab_export_overwrite = True

#-----------------------------------------------
# Setup ExperimentMetadata 
#-----------------------------------------------
omni_anal_logger.info("=== Initializing ExperimentMetadata ===")

# Initialize ExperimentMetadata instance
exp_meta = ExperimentMetadata()

# Generate experiment_metadata.csv if needed
ExperimentMetadata.initialize_experiment_metadata_csv(
    overwrite=experiment_csv_overwrite
)

# Load experiment metadata CSV into memory and group by rat/session for batch processing
exp_meta.load_experiment_metadata_csv()
omni_anal_logger.info(f"Loaded batch session structure:\n{exp_meta.batch_rat_list}")

#raise SystemExit

#-----------------------------------------------
# Batch session preprocessing
#-----------------------------------------------
for rat_block in exp_meta.batch_rat_list:
    rat_id = rat_block["rat_id"]
    for session_name in rat_block["batch_session_list"]:
        omni_anal_logger.info(f"=== Processing session: {rat_id}/{session_name} ===")

        # Step 1: Generate and save session/ephys metadata
        setup_session_metadata(rat_id, session_name, overwrite=metadata_overwrite)

        # Step 2: Extract DIO if needed
        setup_dio_extraction(rat_id, session_name, overwrite=dio_overwrite)

        # Step 3: Compute timestamp sync if needed
        setup_timestamp_sync(
            rat_id,
            session_name,
            dio_channel=dio_channel,
            overwrite=ts_sync_overwrite,
            save_ts_binary=save_ts_sync_binary,
        )

        # Step 4: Export CSC to MATLAB .mat
        cmd = [
            sys.executable,
            "-m",
            "preprocessing.pipeline.export_matlab_csc",
            "--rat_id",
            rat_id,
            "--session_name",
            session_name,
        ]
        if matlab_export_overwrite:
            cmd.append("--overwrite")

        omni_anal_logger.info(f"Exporting MATLAB CSC: {' '.join(cmd)}")
        subprocess.run(cmd, check=True)
