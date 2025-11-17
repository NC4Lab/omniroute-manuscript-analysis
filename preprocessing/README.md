# Preprocessing Environment and Usage

This directory contains the Python preprocessing pipeline for the Omniroute Maze manuscript analysis. It parses raw data (e.g., Trodes/ROS output), aligns streams, and generates intermediate files consumed by the MATLAB analysis code and saves them to `methods_manuscript_datasets` located under `DATASET_ROOT`.

---

## Environment Variables

This repository expects a `.env` file in the project root with the following variables:

```text
TRODES_DIR=...
DATASET_ROOT=...
PYTHONPATH=.
```

- `TRODES_DIR` – path to the local Trodes installation (used for file conversion/utilities)
- `DATASET_ROOT` – root directory containing the Omniroute Maze dataset downloaded from OSF (used for both preprocessing and analysis)
- `PYTHONPATH` – should be set to `.` so local modules are importable when running via `python -m`

---

## Conda Environment

Create and activate the `omniroute_analysis` conda environment from the repository root:

```bash
conda env create -f environment.yml
conda activate omniroute_analysis
```

If `conda` is not recognized in your terminal, initialize it (one time) and restart your shell:

```bash
conda init powershell      # or: conda init bash / zsh
```

---

SpikeInterface Setup
--------------------

The `spikeinterface` package is required and is already included in `environment.yml`.  
We recommend using `spikeinterface` version `0.102.0`, which is the version specified in `environment.yml`, for compatibility with this pipeline.

---

Trodes Setup
------------

Install the Trodes software suite from SpikeGadgets:

* [https://www.spikegadgets.com/trodes](https://www.spikegadgets.com/trodes)
    

Ensure `TRODES_DIR` is set appropriately in the project’s root `.env` file (see main README).

---

Interpreter Selection in VS Code (Optional)
-------------------------------------------

If you are using VS Code, select the conda environment so scripts run with the correct interpreter:

1. Open the Command Palette (`Ctrl+Shift+P`)
    
2. Run **Python: Select Interpreter**
    
3. Choose the interpreter corresponding to `omniroute_analysis`  
    (path typically ends with `.../envs/omniroute_analysis/python.exe`)
    

---

Running Preprocessing Scripts
-----------------------------

Run preprocessing from the **repository root** using the `-m` module flag so imports resolve correctly:

```bash
python -m preprocessing.run_preprocessing
```

In general, prefer:

```bash
python -m package.module
```

over:

```bash
python path/to/script.py
```

Using `-m` ensures the repository root is treated as the top-level package, which avoids common relative-import issues.
