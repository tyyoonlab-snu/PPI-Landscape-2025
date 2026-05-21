# Single-Molecule Image Analysis and ML-Guided Navigation of Antibody Library Landscapes

This repository contains the data analysis, image processing, and predictive modeling code for the manuscript:
> **"Landscape-scale navigation unlocks antibody CDR structural logic for AI-guided rescue and therapeutic optimization"** *(Under Review)*

---

## Overview

This repository hosts **two independent computational tools** developed for the accompanying manuscript on antibody library analysis. The two tools are released together for convenience but are **not coupled**: they address different problems, take different inputs, and can be used entirely separately.

### Tool 1 — MATLAB Single-Molecule Image Analysis GUI (`matlab_scripts/`)

A batch-processing graphical application for quantifying single-molecule fluorescence images (`.tif`). The tool automatically detects fluorescent spots, extracts per-spot intensities, performs image-level background subtraction and Non-Specific Binding (NSB) correction, and classifies each spot's oligomeric state (Monomer / Trimer / Trimer+) by calibrating against a monomeric reference.

In the manuscript, this tool was used to quantify fluorescent spots and verify that **TNF-α bound to Adalimumab and its variants existed in the expected trimeric state**, providing the experimental ground truth for downstream interpretation. The tool itself is general-purpose and can be applied to any single-molecule fluorescence experiment requiring spot detection and oligomeric-state classification.

### Tool 2 — Python ML-Guided Library Navigation Pipeline (`python_pipeline/`)

A Python pipeline that uses experimental measurements of **single and double mutants** to train an additive Ridge regression model. The model learns the contribution (weight) of each single–amino-acid substitution from this low-order data alone, and is then used to predict the fitness of the full combinatorial library *in silico*.

The key result enabled by this tool is that **global top-tier variants can be recovered by screening only a drastically reduced subset of the combinatorial library**, rather than exhaustively measuring every combination. The pipeline also projects the combinatorial sequence space onto a topology-preserving 2D landscape via densMAP (UMAP) and categorizes variants by KMeans clustering to visualize the fitness landscape.

This tool takes a tabular CSV/Excel file of variant measurements as input — its measurements need not have been produced by Tool 1.

---

## Repository Structure

```
PPI-Landscape-2025/
├── matlab_scripts/                          # Tool 1: independent MATLAB image analysis GUI
│   └── spotAnalysisApp_V35_1_Final.m        # GUI app for single-molecule fluorescence image analysis
├── python_pipeline/                         # Tool 2: independent Python ML pipeline
│   ├── analysis_pipeline.py                 # Command-line / IDE executable version
│   └── Analysis_Pipeline_Cleaned.ipynb      # Jupyter / Colab notebook version (same logic)
├── data/                                    # Sample datasets for demo runs (per-tool)
│   ├── demo_image_root/                     # MATLAB demo: root folder of fluorescence images
│   │   ├── Sample_Adalimumab_WT/            #   each subfolder = one well/condition, containing .tif files
│   │   │   ├── WT_07_AVG.tif
│   │   │   ├── WT_08_AVG.tif
│   │   │   └── ...
│   │   ├── Sample_Variant/
│   │   │   └── ...
│   │   └── Control_NSB/                     #   subfolder designated as NSB control
│   │       └── ...
│   └── demo_library_measurements.csv        # Python demo: tabular variant measurements
├── LICENSE                                  # MIT License
├── README.md                                # Project overview and instructions
└── requirements.txt                         # Python dependencies (applies to Tool 2 only)
```

The two tools are released together for distribution convenience but operate independently. Users interested in only one of them need not install or run the other.

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
- MATLAB R2024b (tested) with the Image Processing Toolbox
- Likely compatible with earlier versions (R2021a+) since no R2024b-specific features are used, but this has not been independently verified
- Image Processing Toolbox

### Tested environment

The Python pipeline has been developed and tested **exclusively on Google Colab** (the default free-tier runtime, Linux backend with the standard pre-installed Python 3 environment). We have not independently verified the pipeline on other operating systems.

Because the dependencies (`numpy`, `pandas`, `scikit-learn`, `umap-learn`, `numba`, `matplotlib`, `seaborn`) are all cross-platform Python packages, the pipeline is expected to run on any standard Python 3 environment (Linux, macOS, Windows) with `requirements.txt` installed. Users wishing to run the code locally are welcome to do so, but should be aware that the exact pinned versions of `umap-learn` and `numba` may interact differently with their local Python and BLAS libraries. **For maximum reproducibility, we recommend running the Python pipeline on Google Colab.**

The MATLAB pipeline has been developed and tested on:
- **MATLAB R2024b** with the Image Processing Toolbox

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

**Typical install time on Google Colab**: ~30–60 seconds. Most dependencies (`numpy`, `pandas`, `scikit-learn`, `matplotlib`, `seaborn`) are pre-installed in the Colab runtime, so `pip install -r requirements.txt` effectively only needs to install `umap-learn==0.5.11` and ensure `numba==0.60.0` is present. On a fresh local Python environment, expect ~2–4 minutes for full dependency resolution and installation.

### MATLAB pipeline

No separate installation step is required beyond having MATLAB R2024b (or a compatible version) with the Image Processing Toolbox installed. Simply clone the repository and navigate to `matlab_scripts/` within MATLAB.

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

**Expected run time on demo data** (default Google Colab CPU runtime, dataset size–dependent):

| Pipeline stage | Time (approximate) |
| --- | --- |
| Data loading & preprocessing | < 5 seconds |
| densMAP embedding (`n_neighbors=50`, `densmap=True`) | ~30 sec – 3 min |
| KMeans clustering (`k=5`, on 2D coordinates) | < 5 seconds |
| Ridge regression + recovery curve (100 ratios × 2 metrics) | ~10–30 seconds |
| Figure rendering & saving | < 5 seconds |
| **Total (end-to-end)** | **~1–5 minutes** |

The densMAP step is the dominant cost and scales with the number of unique variants. For reference:
- ~1,000 variants → end-to-end ~1 minute
- ~5,000 variants → end-to-end ~3 minutes
- ~20,000 variants → end-to-end ~5–10 minutes (densMAP alone may take 5–8 minutes)

These estimates are based on the algorithmic complexity of each stage on the default Colab CPU runtime; actual values may vary by ±50% depending on Colab backend allocation.

### MATLAB pipeline demo

A demo dataset is provided as a **root folder** at `data/demo_image_root/`, containing several subfolders. Each subfolder represents one well or condition and holds multiple `.tif` images. This structure is required by the GUI: the application expects a root folder whose immediate subfolders correspond to individual samples/wells.

1. Open MATLAB and navigate to `matlab_scripts/`.
2. Launch the application:
   ```matlab
   spotAnalysisApp_V35_1_Final
   ```
3. In the GUI, click **Select Root** and choose the `data/demo_image_root/` folder. The subfolders inside it will automatically populate the **Process** / **NSB Control** / **Display** list panels.
4. Select the subfolders to process and the subfolder to use as NSB control, then run the analysis with default parameters.

**Expected output**:
- An `.xlsx` report containing:
  - Aggregated per-well statistics (mean intensity, spot counts, monomer/trimer/trimer+ classification)
  - Raw intensity distributions for downstream analysis
  - QC summary (saturated spots, specific spot counts, NSB correction values)

**Expected run time** (per image / per dataset):
- Single `.tif` image stack analysis (spot detection + intensity extraction + classification): ~5–30 seconds per image, depending on image size and spot density
- Full batch processing of one experiment (e.g., 96-well plate, ~96 images): ~2–5 minutes
- NSB correction and `.xlsx` report generation: < 30 seconds

These estimates are based on the cost of the underlying image-processing operations (`imgaussfilt` → adaptive `imbinarize` → `imopen` → `regionprops`) which are linear in image area and run on a single CPU thread.

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

1. Organize your `.tif` image files into a **root folder containing one subfolder per well/condition**, with the `.tif` images for that well/condition placed inside the corresponding subfolder. Subfolder names starting with `.` or `_`, and any folder named `_processed`, are automatically ignored by the GUI.
   ```
   my_experiment/                ← root folder you will select in the GUI
   ├── Well_A1/
   │   ├── img_001.tif
   │   └── img_002.tif
   ├── Well_A2/
   │   └── ...
   └── NSB_Control/
       └── ...
   ```
2. Launch `spotAnalysisApp_V35_1_Final` in MATLAB.
3. Click **Select Root** and choose the root folder of your experiment. Subfolders will appear in the **Process** / **NSB Control** / **Display** list panels.
4. Configure detection parameters (Gaussian smoothing radius, adaptive sensitivity, morphological opening, area filter) using the GUI sliders. The live preview pane assists with parameter tuning before batch processing.
5. Select which subfolders to use as NSB controls and (optionally) as the internal monomer anchor for calibration.
6. Run batch processing; results are exported as `.xlsx` to a user-specified output directory.

### Reproduction of manuscript results *(optional)*

The two tools contribute to different parts of the manuscript and are reproduced independently.

**Python ML pipeline results** (combinatorial library navigation and top-variant recovery curves):
1. Obtain the full processed dataset from the corresponding author (see **Data Availability** below).
2. Place it in `data/full_library_measurements.csv`.
3. Run `Analysis_Pipeline_Cleaned.ipynb` end-to-end with default parameters and fixed random seeds (already set in Section 1 of the notebook).

**MATLAB image analysis results** (TNF-α trimer-state quantification on Adalimumab variants):
1. Obtain the raw `.tif` image stacks organized in the required root/subfolder structure from the corresponding author (see **Data Availability** below).
2. Launch `spotAnalysisApp_V35_1_Final` and use **Select Root** to load the root directory.
3. Use the parameter values reported in the manuscript Methods section to reproduce the reported spot statistics.

---

## Pipeline Structure (Python)

Both `Analysis_Pipeline_Cleaned.ipynb` and `analysis_pipeline.py` are organized into four sequential sections:

1. **Libraries & Imports** — installs pinned dependencies and fixes random seeds.
2. **Data Upload & Preprocessing** — applies Data Quality + biological filtering and aggregates replicates.
3. **UMAP Embedding & KMeans Clustering** — projects the combinatorial space to 2D and assigns clusters.
4. **Predictive Recovery from Low-Order (≤ Double) Mutants** — trains the additive Ridge model and reports top 1% recovery curves.

### Key methodological features

- **Biological preprocessing & filtering** — curates experimental measurements by strictly enforcing biological constraints (e.g., masking occupancy/affinity values for variants that fail to express, so non-expressing variants do not contaminate the binding signal).
- **Physicochemical vectorization & densMAP embedding** — transforms combinatorial mutations into high-dimensional numerical vectors based on position-specific physicochemical properties (Kyte–Doolittle hydropathy, residue volume, and isoelectric point). The sequence space is then projected onto a topology-preserving 2D landscape with `densMAP` (UMAP) and systematically categorized via KMeans clustering, allowing visual inspection of the fitness landscape.
- **ML-guided top-variant recovery from reduced screening** — the core insight is that a Ridge regression model trained **only on single and double mutants** can learn per-position amino-acid weights that capture the fundamental contributions and pairwise epistatic interactions driving fitness. Using these learned weights, the model predicts the fitness of the full combinatorial library *in silico* and ranks variants accordingly. This enables effective recovery of the global top 1% variants in both affinity and productivity **without exhaustive measurement of the combinatorial library** — only the low-order (≤ double) mutant subset needs to be experimentally screened.

### MATLAB image analysis key features

- **Automated spot detection** — identifies fluorescent spots using adjustable parameters (Gaussian smoothing, adaptive sensitivity, morphological opening, and area filtering).
- **Intensity profiling & oligomeric classification** — extracts mean intensities, performs image-level background subtraction, and calibrates against an internal/external monomeric reference to classify spots into distinct oligomeric states (Monomer, Trimer, Trimer+). In the accompanying manuscript, this was used to confirm that TNF-α bound to Adalimumab and its variants existed predominantly as trimers.
- **NSB correction & quality control** — automatically corrects for Non-Specific Binding (NSB) using user-designated control wells and provides a comprehensive QC summary (saturated spots, specific spot counts, spot-area distributions).
- **Batch export** — generates detailed `.xlsx` reports containing both aggregated per-well statistics and raw intensity distributions for downstream analysis.

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
