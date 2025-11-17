from pathlib import Path
from dotenv import load_dotenv
import os

# Start from this file's directory and walk upward to find .env
env_path = None
for parent in Path(__file__).resolve().parents:
    candidate = parent / ".env"
    if candidate.exists():
        env_path = candidate
        break

if not env_path:
    raise FileNotFoundError(
        f".env file not found in any parent directory of {__file__}"
    )

# Load environment variables
load_dotenv(dotenv_path=env_path)

# Read required vars
TRODES_DIR = os.getenv("TRODES_DIR")
DATASET_ROOT = os.getenv("DATASET_ROOT")

# Validate presence
if not TRODES_DIR or not DATASET_ROOT:
    raise EnvironmentError(
        "TRODES_DIR and DATASET_ROOT must be set in the .env file."
    )
