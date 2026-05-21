# Sequence Embedding and Machine Learning-Guided Antibody Library Analysis

This repository contains the data analysis, image processing, and predictive modeling code for the manuscript:
> **"Landscape-scale navigation unlocks antibody CDR structural logic for AI-guided rescue and therapeutic optimization"** *(Under Review)*

---

## Overview

This repository provides a complete computational framework for antibody library optimization, spanning from raw single-molecule image quantification to AI-driven sequence design. By integrating single-molecule biophysics with machine learning–guided predictive modeling, we establish a rational strategy to analyze fitness landscapes and discover top-tier variants without the need for exhaustive experimental screening.

The repository is organized into two complementary components:

1. **MATLAB Image Analysis GUI** (`matlab_scripts/`) — processes single-molecule fluorescence images to extract intensity profiles and quantify antigen-binding occupancy.
2. **Python ML Pipeline** (`python_pipeline/`) — performs physicochemical sequence embedding and trains predictive models to navigate the combinatorial fitness landscape.

---

## Repository Structure

```
PPI-Landscape-2025/
├── matlab_scripts/
│   └── spotAnalysisApp_V35_1_Final.m        # MATLAB GUI for single-molecule fluorescence image analysis
├── python_pipeline/
│   ├── analysis_pipeline.py                 # Command-line / IDE executable version
│   └── Analysis_Pipeline_Cleaned.ipynb      # Jupyter / Colab notebook version (same logic)
├── data/                                    # Sample datasets for demo runs
│   ├── demo_image_stack.tif                 # Example single-molecule fluorescence image (MATLAB demo)
│   └── demo_library_measurements.csv        # Example processed library data (Python demo)
├── LICENSE                                  # MIT License
├── README.md                                # Project overview and instructions
└── requirements.txt                         # Python dependencies (pinned versions)
```

---

## System Requirements

### Software dependencies

**Python pipeline** (see `requirements.txt`):
- Python 3 (tested on the default Google Colab Python runtime)
- numpy
- pandas
- matplotlib
- seaborn
- scikit-learn
- umap-learn == 0.5.11  *(with densMAP support)*
- numba == 0.60.0

`umap-learn` and `numba` are pinned to the versions used during development to ensure reproducibility of the densMAP embedding step. Other dependencies follow whichever versions are resolved by `pip` at install time; this matches the behavior expected on a fresh Google Colab runtime.

**MATLAB pipeline**:
- MATLAB R2024b or later
- Image Processing Toolbox

### Tested environment

The Python pipeline has been developed and tested **exclusively on Google Colab** (the default free-tier runtime, Linux backend with the standard pre-installed Python 3 environment). We have not independently verified the pipeline on other operating systems.

Because the dependencies (`numpy`, `pandas`, `scikit-learn`, `umap-learn`, `numba`, `matplotlib`, `seaborn`) are all cross-platform Python packages, the pipeline is expected to run on any standard Python 3 environment (Linux, macOS, Windows) with `requirements.txt` installed. Users wishing to run the code locally are welcome to do so, but should be aware that the exact pinned versions of `umap-learn` and `numba` may interact differently with their local Python and BLAS libraries. **For maximum reproducibility, we recommend running the Python pipeline on Google Colab.**

The MATLAB pipeline has been tested on:
- MATLAB R2024b or later with the Image Processing Toolbox

### Hardware requirements

**Non-standard hardware is not required.** The Python pipeline was developed and tested on the **default Google Colab free-tier runtime** (no GPU/TPU allocation; standard CPU instance with ~12 GB RAM). All steps — physicochemical vectorization, densMAP embedding, KMeans clustering, and Ridge regression — run on CPU within this default environment.

For users running locally:
- A standard desktop or laptop with **≥ 8 GB RAM** is sufficient for the demo dataset.
- The MATLAB image analysis is single-threaded; batch processing benefits from multi-core CPUs but does not require them.
- No GPU is needed for any component of this pipeline.

---

## Installation Guide

### Python pipeline

We recommend an isolated environment (Anaconda or `venv`). To ensure reproducibility and prevent API drift, install the exact pinned dependencies:

```bash
git clone https://github.com/tyyoonlab-snu/PPI-Landscape-2025.git
cd PPI-Landscape-2025
pip install -r requirements.txt
```

**Typical install time on Google Colab**: ~1–2 minutes (most dependencies are pre-installed on the Colab runtime; only `umap-learn` and pinned versions of `numba` need to be installed). On a local machine with a fresh Python environment, expect ~3–5 minutes depending on internet connection.

### MATLAB pipeline

No separate installation step is required beyond having MATLAB R2024b+ with the Image Processing Toolbox installed. Simply clone the repository and navigate to `matlab_scripts/` within MATLAB.

**Typical setup time**: < 1 minute (assuming MATLAB and the Image Processing Toolbox are already installed).

---

## Demo

### Python pipeline demo

A small example dataset is provided in `data/demo_library_measurements.csv`. To run the demo:

#### Option A — Jupyter / Colab notebook *(recommended)*

```bash
jupyter notebook python_pipeline/Analysis_Pipeline_Cleaned.ipynb
```

Or open the notebook directly in Google Colab via the **Open in Colab** button on GitHub.

#### Option B — Python script

```bash
python python_pipeline/analysis_pipeline.py
```

**Expected output**:
- 2D densMAP embedding scatter plot showing the projected combinatorial sequence space, colored by KMeans clusters
- Predictive recovery curve showing top-1% champion recovery rate as a function of training set size (≤ double mutants)
- Tabular summary of cluster-wise mean affinity / productivity statistics (printed to console or notebook cell)
- A `results/` directory containing PNG figures and CSV tables of the embedding coordinates and model predictions

**Expected run time on demo data** (default Google Colab CPU runtime):
- Full notebook execution: ~5–10 minutes
- Most time is spent on densMAP embedding (~3–5 min) and Ridge model cross-validation (~2 min)

### MATLAB pipeline demo

A demo `.tif` image stack is provided in `data/demo_image_stack.tif`.

1. Open MATLAB and navigate to `matlab_scripts/`.
2. Launch the application:
   ```matlab
   spotAnalysisApp_V35_1_Final
   ```
3. In the GUI, load `../data/demo_image_stack.tif` and run with default parameters.

**Expected output**:
- An `.xlsx` report containing:
  - Aggregated per-well statistics (mean intensity, spot counts, monomer/trimer/trimer+ classification)
  - Raw intensity distributions for downstream analysis
  - QC summary (saturated spots, specific spot counts, NSB correction values)

**Expected run time on demo image**:
- ~2–3 minutes per image stack on a standard desktop

---

## Instructions for Use

### Running the Python pipeline on your own data

1. Prepare your library measurement file as a CSV following the schema of `data/demo_library_measurements.csv`. Required columns:
   - `sequence` — variant sequence or mutation identifier
   - `affinity` — measured binding affinity (or occupancy)
   - `productivity` — measured expression / productivity
   - `n_mutations` — number of mutations from wild-type
2. Edit the **User inputs** cell (notebook) or top-of-file constants (`analysis_pipeline.py`) to point to your CSV.
3. Adjust filtering thresholds (e.g., minimum replicate count, expression cutoff) as appropriate for your dataset.
4. Run the four pipeline stages sequentially. Each stage's output is the input for the next, so re-running from the top is safe.

### Running the MATLAB pipeline on your own data

1. Place your `.tif` image stacks in a working directory.
2. Launch `spotAnalysisApp_V35_1_Final` in MATLAB.
3. Configure detection parameters (Gaussian smoothing radius, adaptive sensitivity, morphological opening, area filter) using the GUI sliders.
4. Specify control wells for NSB correction.
5. Run batch processing; results are exported as `.xlsx` to a user-specified output directory.

### Reproduction of manuscript results *(optional)*

To reproduce the quantitative results reported in the manuscript:
1. Obtain the full processed dataset from the corresponding author (see **Data Availability** below).
2. Place it in `data/full_library_measurements.csv`.
3. Run `Analysis_Pipeline_Cleaned.ipynb` end-to-end with default parameters and fixed random seeds (already set in Section 1 of the notebook).

---

## Pipeline Structure (Python)

Both `Analysis_Pipeline_Cleaned.ipynb` and `analysis_pipeline.py` are organized into four sequential sections:

1. **Libraries & Imports** — installs pinned dependencies and fixes random seeds.
2. **Data Upload & Preprocessing** — applies Data Quality + biological filtering and aggregates replicates.
3. **UMAP Embedding & KMeans Clustering** — projects the combinatorial space to 2D and assigns clusters.
4. **Predictive Recovery from Low-Order (≤ Double) Mutants** — trains the additive Ridge model and reports top 1% recovery curves.

### Key methodological features

- **Biological preprocessing & filtering** — curates experimental measurements by strictly enforcing biological constraints (e.g., masking occupancy/affinity values for variants that fail to express).
- **Physicochemical vectorization & densMAP embedding** — transforms combinatorial mutations into high-dimensional numerical vectors based on position-specific physicochemical properties (Kyte–Doolittle hydropathy, residue volume, and isoelectric point). The sequence space is then projected onto a topology-preserving 2D landscape with `densMAP` (UMAP) and systematically categorized via KMeans clustering.
- **ML-guided predictive recovery** — bypasses brute-force screening by training an additive Ridge regression model exclusively on low-order (≤ double) mutants. By capturing fundamental pairwise epistatic interactions, the model drastically condenses the search space and accurately recovers the global top 1% champions in both affinity and productivity.

### MATLAB image analysis key features

- **Automated spot detection** — identifies GFP spots using adjustable parameters (Gaussian smoothing, adaptive sensitivity, morphological opening, and area filtering).
- **Intensity profiling & oligomeric classification** — extracts mean intensities, performs image-level background subtraction, and calibrates against an internal/external monomeric reference to classify spots into distinct states (Monomer, Trimer, Trimer+).
- **NSB correction & quality control** — automatically corrects for Non-Specific Binding (NSB) using control wells and provides a comprehensive QC summary.
- **Batch export** — generates detailed `.xlsx` reports containing both aggregated statistics and raw intensity distributions.

---

## Data Availability

Representative processed datasets for demo purposes are provided in `data/`. Raw image stacks and the full experimental dataset are available from the corresponding author upon reasonable request.

---

## Citation

If you use this code or the analysis framework in your work, please cite:

```
(Citation to be updated upon publication)
```

---

## License

This project is released under the MIT License. See [`LICENSE`](LICENSE) for details.

---

## Contact

For questions about the code or the manuscript, please open an issue on this repository or contact the corresponding author.
