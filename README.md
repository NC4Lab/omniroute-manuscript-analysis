# Omniroute Manuscript Analysis  
Preprocessing and analysis pipeline for the Omniroute behavioral validation manuscript.

<p align="center">
  <img src="assets/omniroute_render.png" alt="Omniroute Render" width="650">
</p>

## Paper
This repository contains the full preprocessing and analysis code used for:

**Lester AW, Mombeini AG, Madhav MS (2025).  
*The Omniroute maze: a novel rodent navigation apparatus that integrates dynamically configurable routes, sensory cues, and automated reward delivery.*  
In review at eLife.**

All figures and quantitative results in the manuscript are generated from this codebase.

---

## Data Availability

The pipeline supports two modes:
1. **Full pipeline:** raw → ros/ephys preprocessing → analysis  
2. **Analysis-only:** skip preprocessing by using the preprocessed data provided on OSF  

All raw and preprocessed data required to reproduce the analyses and figures are in `methods_manuscript_datasets` hosted on the Omniroute OSF project:

**OSF project:** https://osf.io/xfpsv/

### Reconstructing `methods_manuscript_datasets`
Before running any preprocessing or analysis, you must reconstruct the `methods_manuscript_datasets` directory from the multi-part archive provided on OSF.

1. Install 7-Zip or another tool that supports multi-part `.7z` archives.  
   For Windows, 7-Zip is available at https://7-zip.org.
2. Create a directory with read and write access to store the dataset, for example `/path/to/data/`. 
3. Download **all** parts of the archive from OSF into that directory  
   (files named `methods_manuscript_datasets.7z.001`, `methods_manuscript_datasets.7z.002`, and so on).
4. Use your archive tool (e.g., 7-Zip) to extract **only** `methods_manuscript_datasets.7z.001` into `/path/to/data/`.  
   For example, on Windows with 7-Zip: right-click → **7-Zip → Extract Here**.  
   The tool will automatically read `.002`, `.003`, etc. and reconstruct  
   `/path/to/data/methods_manuscript_datasets/`.
5. Set the `DATASET_ROOT` variable in the root `.env` file to the reconstructed `methods_manuscript_datasets` directory.  
   Example: `DATASET_ROOT=/path/to/data/methods_manuscript_datasets`

---

## Repository Structure (Top Level)

```

omniroute-manuscript-analysis/  
│  
├── preprocessing/              # Python pipeline for extracting and aligning raw data
│   ├── run_preprocessing.py    # Main entry point for preprocessing
│   └── README.md               # Detailed preprocessing setup
│ 
├── analysis/           # MATLAB scripts generating all manuscript figures and stats  
│ └── run_analysis.m    # Main entry point for MATLAB analysis  
│ 
├── results/  
│ ├── figures/          # Outputs: .png/.pdf/.fig  
│ ├── tables/           # Outputs: .mat/.csv and intermediate summaries  
│ └── logs/             # Preprocessing and analysis logs  
│  
├── environment.yml     # Conda environment for Python preprocessing  
└── README.md

```

---

## Python Requirements (Preprocessing)

Preprocessing is performed in **Python 3** using the conda environment provided:

- Python ≥ 3.11  
- numpy, pandas, scipy  
- rosbags  
- python-dotenv  
- spikeinterface[full] (0.102.0)  

Create and activate the environment:

```bash
conda env create -f environment.yml
conda activate omniroute_analysis
```

**Note:**  
See `preprocessing/README.md` for full details.

---

MATLAB Requirements (Analysis)
------------------------------

Analyses and figure generation were tested using:

* **MATLAB R2022b**
    
* Statistics and Machine Learning Toolbox
    
* Image Processing Toolbox
    
* Signal Processing Toolbox
    
* Curve Fitting Toolbox
    

No third-party MATLAB dependencies are required.

---

Preprocessing (Python)
----------------------

The preprocessing pipeline parses raw ROS bag data, synchronizes sensor streams, extracts task-relevant events, and generates intermediate continuous spike channel (CSC) files consumed by MATLAB for the `ephys` analysis.

**Entry point:**

```bash
python -m preprocessing.run_preprocessing
```

After the run, processed continuous spike channel (CSC) data will be written to `methods_manuscript_datasets`.
    

**Note:**  
This step is optional if using the preprocessed dataset from OSF.

---

Analysis (MATLAB)
-----------------

The MATLAB pipeline generates all figures and quantitative analyses reported in the manuscript.

**Entry point:**  
Open MATLAB and run:

```matlab
analysis/run_analysis.m
```

This script:

* Loads required data from `methods_manuscript_datasets` located under `DATASET_ROOT`
* Performs all relevant analyses
* Generates and saves all manuscript figures
* Writes summary tables and logs to the `results/` directory
    

**Outputs will be placed in:**

* `results/figures/`
* `results/tables/`
* `results/logs/analysis_*.txt`
    

**Note:**  
Descriptive and summary statistics are stored in the associated log txt files.

---

End-to-End Reproduction Instructions
------------------------------------

1. **Clone the repository**  

```bash
git clone https://github.com/NC4Lab/omniroute-manuscript-analysis.git
cd omniroute-manuscript-analysis
```
# Checkout the frozen, manuscript-specific branch used for the published analyses
```bash
git checkout manuscript-omniroute-lester
```

2. **Download `methods_manuscript_datasets` from OSF**  
    Place the folder in any directory with read/write access, and set the `DATASET_ROOT` variable in the root `.env` file.
    
3. **(Optional) Run preprocessing**
    

```bash
conda activate omniroute_analysis
python -m preprocessing.run_preprocessing
```

4. **Run the MATLAB analysis**
    

Open MATLAB and execute:

```matlab
analysis/run_analysis.m
```

This will regenerate all paper figures and tables.

---

Licensing
---------

Analysis code in this repository is released under the MIT License (see LICENSE).
    

---
