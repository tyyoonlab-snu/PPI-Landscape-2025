# Single-Molecule Image Analysis and ML-Guided Navigation of Antibody Library Landscapes

This repository contains the data analysis, image processing, and predictive modeling code for the manuscript:
> **"Landscape-scale navigation unlocks antibody CDR structural logic for AI-guided rescue and therapeutic optimization"** *(Under Review)*

---

## Overview

This repository hosts **two independent computational tools** developed for the accompanying manuscript on antibody library analysis. The two tools are released together for convenience but are **not coupled**: they address different problems, take different inputs, and can be used entirely separately.

### Tool 1 — MATLAB Single-Molecule Image Analysis GUI (`matlab_scripts/`)

A batch-processing graphical application for quantifying single-molecule fluorescence images (`.tif`). The tool automatically detects fluorescent spots, extracts per-spot intensities, performs image-level background subtraction and Non-Specific Binding (NSB) correction, and classifies each spot's oligomeric state (Monomer / Trimer / Trimer+) by calibrating against a monomeric reference.

In the manuscript, this tool was used to quantify fluorescent spots and verify that **TNF-α bound to Adalimumab and its variants existed in the expected trimeric state**, providing the experimental ground truth for downstream interpretation. The tool itself is general-purpose and can be applied to any single-molecule fluorescence experiment requiring spot detection and oligomeric-state classification.

### Tool 2 — Python Library Landscape Analysis Pipeline (`python_pipeline/`)

A Python pipeline that turns tabular variant measurements into the fitness landscapes reported in the manuscript, and navigates them.

The pipeline has three outputs:

1. **A topology-preserving 2D landscape.** Combinatorial variants are encoded by position-specific physicochemical properties and projected with densMAP (UMAP).
2. **Topographical regions of that landscape.** Local fitness statistics computed over each variant's nearest neighbours are used to delineate peak clusters, rugged interfaces and broad valleys, and continuous surfaces are rendered for visualization.
3. **ML-guided recovery of top variants.** An additive Ridge regression model is trained on **single and double mutants only**, learning the contribution (weight) of each single–amino-acid substitution from this low-order data alone, and is then used to predict the fitness of the full combinatorial library *in silico*.

The key result enabled by (3) is that **global top-tier variants can be recovered by screening only a drastically reduced subset of the combinatorial library**, rather than exhaustively measuring every combination.

This tool takes a tabular CSV/Excel file of variant measurements as input — its measurements need not have been produced by Tool 1.

---

## Repository Structure

```
PPI-Landscape-2025/
├── matlab_scripts/                            # Tool 1: independent MATLAB image analysis GUI
│   └── spotAnalysisApp_V35_1_Final.m          # GUI app for single-molecule fluorescence image analysis
├── python_pipeline/                           # Tool 2: independent Python analysis pipeline
│   ├── Analysis_Pipeline_Cleaned.ipynb        # Reference implementation (Jupyter / Colab)
│   ├── analysis_pipeline.py                   # Command-line / IDE equivalent
│   ├── topographical_classification.py        # Local-statistics delineation of landscape regions
│   ├── fitness_landscape_visualization.py     # Continuous 3D surface and 2D top-view rendering
│   ├── l108g_control_landscape.py             # L108G-fixed control landscape (build, project, render)
│   └── l108g_macro_comparison.py              # Functional vs L108G-control distribution comparison
├── data/                                      # Datasets for demo runs (per-tool)
│   ├── demo_image_root/                       # MATLAB demo: root folder of fluorescence images
│   │   ├── Sample_Adalimumab_WT/              #   each subfolder = one well/condition, containing .tif files
│   │   │   ├── WT_07_AVG.tif
│   │   │   ├── WT_08_AVG.tif
│   │   │   └── ...
│   │   ├── Sample_Variant/
│   │   │   └── ...
│   │   └── Control_NSB/                       #   subfolder designated as NSB control
│   │       └── ...
│   ├── demo_library_measurements.csv          # Python demo: main combinatorial library (9,588 variants)
│   └── demo_L108G_library_measurements.csv    # Python demo: L108G-fixed control library (1,910 variants)
├── LICENSE                                    # MIT License
├── README.md                                  # Project overview and instructions
└── requirements.txt                           # Python dependencies (applies to Tool 2 only)
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
- scipy
- scikit-learn
- plotly
- umap-learn == 0.5.11  *(with densMAP support)*
- numba == 0.60.0

`umap-learn` and `numba` are pinned to the versions used during development to ensure reproducibility of the densMAP embedding step. Other dependencies follow whichever versions are resolved by `pip` at install time; this matches the behavior expected on a fresh Google Colab runtime.

**MATLAB pipeline**:
- MATLAB R2024b (tested) with the Image Processing Toolbox
- Likely compatible with earlier versions (R2021a+) since no R2024b-specific features are used, but this has not been independently verified

### Tested environment

The Python pipeline has been developed and tested **exclusively on Google Colab** (the default free-tier runtime, Linux backend with the standard pre-installed Python 3 environment). We have not independently verified the pipeline on other operating systems.

Because the dependencies are all cross-platform Python packages, the pipeline is expected to run on any standard Python 3 environment (Linux, macOS, Windows) with `requirements.txt` installed. Users wishing to run the code locally are welcome to do so, but should be aware that the exact pinned versions of `umap-learn` and `numba` may interact differently with their local Python and BLAS libraries. **For maximum reproducibility, we recommend running the Python pipeline on Google Colab.**

The MATLAB pipeline has been developed and tested on:
- **MATLAB R2024b** with the Image Processing Toolbox

### Hardware requirements

**Non-standard hardware is not required.** The Python pipeline was developed and tested on the **default Google Colab free-tier runtime** (no GPU/TPU allocation; standard CPU instance with ~12 GB RAM). All steps — physicochemical vectorization, densMAP embedding, clustering, topographical classification, surface rendering, and Ridge regression — run on CPU within this default environment.

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

**Typical install time on Google Colab**: ~30–60 seconds. Most dependencies (`numpy`, `pandas`, `scipy`, `scikit-learn`, `matplotlib`, `seaborn`, `plotly`) are pre-installed in the Colab runtime, so `pip install -r requirements.txt` effectively only needs to install `umap-learn==0.5.11` and ensure `numba==0.60.0` is present. On a fresh local Python environment, expect ~2–4 minutes for full dependency resolution and installation.

### MATLAB pipeline

No separate installation step is required beyond having MATLAB R2024b (or a compatible version) with the Image Processing Toolbox installed. Simply clone the repository and navigate to `matlab_scripts/` within MATLAB.

**Typical setup time**: < 1 minute (assuming MATLAB and the Image Processing Toolbox are already installed).

---

## Demo

### Python pipeline demo

The processed measurements for the adalimumab HCDR2/HCDR3 combinatorial library (9,588 variants) are provided in `data/demo_library_measurements.csv`. To run the demo:

#### Option A — Jupyter / Colab notebook *(recommended)*

```bash
jupyter notebook python_pipeline/Analysis_Pipeline_Cleaned.ipynb
```

Or open the notebook directly in Google Colab via the **Open in Colab** button on GitHub. Run the sections in order; the upload cell in Section 2 accepts `data/demo_library_measurements.csv`.

#### Option B — Python scripts

```bash
python python_pipeline/analysis_pipeline.py \
    --input data/demo_library_measurements.csv --outdir output

python python_pipeline/topographical_classification.py \
    --input output/1_UMAP_2D_Coordinates.csv \
    --output output/2_Topographical_Classification.csv

python python_pipeline/fitness_landscape_visualization.py \
    --input output/2_Topographical_Classification.csv --outdir output
```

To reproduce the L108G-fixed control analysis (Fig. 3p–s, Extended Data Fig. 7):

```bash
python python_pipeline/l108g_control_landscape.py \
    --input data/demo_L108G_library_measurements.csv --outdir output_l108g

python python_pipeline/l108g_macro_comparison.py \
    --main data/demo_library_measurements.csv \
    --l108g data/demo_L108G_library_measurements.csv --outdir output_l108g
```

**Expected output** (written to `output/`):

| File | Contents |
| --- | --- |
| `1_UMAP_2D_Coordinates.csv` | densMAP coordinates and cluster assignment per variant |
| `2_Topographical_Classification.csv` | Local fitness statistics and topographical class per variant |
| `Landscape3D_<metric>.html` | Interactive 3D fitness surface |
| `TopView_<metric>.svg`, `TopView_Combined.svg` | 2D top-view contour maps |
| `SourceData_Surface_<metric>.csv` | Interpolated matrix underlying each rendered surface |
| `PanelB_Recovery.svg` | Top-1% recovery curve as a function of screening depth |

A tabular summary of class-wise affinity and productivity statistics is printed to the console or notebook cell.

**Expected run time on demo data** (default Google Colab CPU runtime, dataset size–dependent):

| Pipeline stage | Time (approximate) |
| --- | --- |
| Data loading & preprocessing | < 5 seconds |
| densMAP embedding (`n_neighbors=50`, `densmap=True`) | ~30 sec – 3 min |
| KMeans clustering (silhouette-selected *k*, on 2D coordinates) | ~5–20 seconds |
| Topographical classification (30-NN local statistics + DBSCAN) | ~5–15 seconds |
| Landscape rendering (150 × 150 grid, per metric) | ~10–20 seconds |
| Ridge regression + recovery curve (100 ratios × 2 metrics) | ~10–30 seconds |
| Figure rendering & saving | < 5 seconds |
| **Total (end-to-end)** | **~2–6 minutes** |

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

1. Prepare your library measurement file as a CSV or Excel file following the schema of `data/demo_library_measurements.csv`.

   **Required:**
   - `Mutation_Description` — variant identifier, given as underscore-separated single-mutation codes (e.g. `T52V_S55H_S100K`); the parental sequence is labelled `Parental`. If this column is absent, the first column of the file is used.
   - At least one measurement column from `Normalized_Occupancy` (binding, normalized to parental), `Normalized_Expression` (productivity, normalized to parental), or `MPNN Z-Score`.

   **Optional:**
   - `Data_Quality` — per-measurement quality label. Rows labelled `Valid`, `Valid (Non-binding)` or `Valid; Valid (Non-binding)` are retained for all analyses; rows labelled `Invalid (Low expression)` are retained for expression analysis only, and their occupancy values are masked, since binding cannot be inferred from a non-expressing variant. If the column is absent, no quality filtering is applied.
   - `*_Std` columns — per-variant standard deviations across replicates.
   - Replicate rows sharing the same `Mutation_Description` are aggregated to mean and standard deviation automatically.

2. Update the combinatorial design specification. The physicochemical encoding is defined by the `design_rules` dictionary in Section 3 of the notebook (`analysis_pipeline.py`: `perform_umap_clustering`), which lists each mutated position, its parental residue and the permitted substitutions. **This must match your library**, since the embedding is constructed from it rather than parsed from the data file.

3. Adjust filtering thresholds and classification parameters as appropriate for your dataset. The topographical thresholds are exposed as module-level constants at the top of `topographical_classification.py`.

4. Run the pipeline sections in order. Each section's output is the input for the next, so re-running from the top is safe.

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

**Python pipeline results** (fitness landscapes, topographical regions, and top-variant recovery curves):
1. Use `data/demo_library_measurements.csv`, which contains the processed measurements for the full combinatorial library reported in the manuscript.
2. Run `Analysis_Pipeline_Cleaned.ipynb` end-to-end with default parameters and fixed random seeds (already set in Section 1 of the notebook).

**MATLAB image analysis results** (TNF-α trimer-state quantification on Adalimumab variants):
1. Obtain the raw `.tif` image stacks organized in the required root/subfolder structure from the corresponding author (see **Data Availability** below).
2. Launch `spotAnalysisApp_V35_1_Final` and use **Select Root** to load the root directory.
3. Use the parameter values reported in the manuscript Methods section to reproduce the reported spot statistics.

---

## Pipeline Structure (Python)

The pipeline is organized into six sequential sections. `Analysis_Pipeline_Cleaned.ipynb` is the reference implementation; the `.py` modules expose the same steps for command-line use.

1. **Libraries & Imports** — installs pinned dependencies and fixes random seeds (`random_state = 42` throughout).
2. **Data Upload & Preprocessing** — applies quality and biological filtering and aggregates replicates.
3. **densMAP Embedding & KMeans Clustering** — projects the combinatorial space to 2D and assigns exploratory clusters.
4. **Topographical Classification** — delineates peak clusters, rugged interfaces and broad valleys from local fitness statistics.
5. **Fitness Landscape Rendering** — builds the continuous 3D surfaces and 2D top-view contour maps.
6. **Predictive Recovery from Low-Order (≤ Double) Mutants** — trains the additive Ridge model and reports top-1% recovery curves.

> **Note on Section 3.** The KMeans step provides an exploratory partition of the embedding and is retained for completeness. The topographical regions reported in the manuscript are assigned in Section 4, not by KMeans.

#### L108G-fixed control analysis (Fig. 3p–s, Extended Data Fig. 7)

Two additional scripts support the L108G control experiment, which tests whether pre-filtering the library to individually functional substitutions is necessary. They are independent of the six-section pipeline above and operate on their own dataset (`data/demo_L108G_library_measurements.csv`, 1,910 measured variants on the L108G background):

- **`l108g_control_landscape.py`** builds the 1,920-variant L108G-fixed combinatorial manifold, projects the measured variants onto it, and renders the affinity and expression landscapes.
- **`l108g_macro_comparison.py`** compares the distributions of affinity and expression, and the complete-loss fractions, between the main functional-variant library and the L108G-fixed control.

### Key methodological features

- **Biological preprocessing & filtering** — curates experimental measurements by strictly enforcing biological constraints (e.g., masking occupancy/affinity values for variants that fail to express, so non-expressing variants do not contaminate the binding signal).

- **Physicochemical vectorization & densMAP embedding** — transforms combinatorial mutations into high-dimensional numerical vectors based on position-specific physicochemical properties (Kyte–Doolittle hydropathy, residue volume in Å³, and isoelectric point), with parental residues placed at the origin of the property space. The standardized vectors are projected onto a topology-preserving 2D landscape with `densMAP` (UMAP; `n_neighbors=50`, `min_dist=0.1`, `spread=2.0`).

- **Topographical classification** — for each variant, the local mean and local standard deviation of log₂-normalized relative occupancy are computed over its 30 nearest neighbours in the embedding. Regions of high local variability are assigned as rugged interfaces; among the remainder, the highest and lowest local means define peak clusters and broad valleys. A density-based filter (DBSCAN) retains only spatially contiguous members of each class, so that the reported regions correspond to coherent areas of the landscape rather than scattered points. Variants with no detectable binding are assigned a floor occupancy of 0.05 rather than being discarded, so that non-binding regions still contribute to their neighbourhood statistics.

- **Continuous landscape rendering** — measured metrics are interpolated onto a regular 150 × 150 grid by distance-weighted k-nearest-neighbour regression (*k* = 15) followed by Gaussian smoothing (σ = 1.2 grid units), and expressed in standard-deviation units relative to the parental clone.

  > **Interpolated surfaces.** The grid extends 1.5 densMAP units beyond the sampled range and no convex-hull or density masking is applied, so parts of a rendered surface lying outside the sampled sequence space are inferred rather than measured. All quantitative claims about landscape ruggedness and navigability in the manuscript are computed on the discrete measured variants, not on these surfaces. The interpolated matrix underlying each surface is written to `SourceData_Surface_<metric>.csv`.

- **ML-guided top-variant recovery from reduced screening** — the core insight is that a Ridge regression model trained **only on single and double mutants** can learn per-position amino-acid weights that capture the fundamental contributions and pairwise epistatic interactions driving fitness. Using these learned weights, the model predicts the fitness of the full combinatorial library *in silico* and ranks variants accordingly. This enables effective recovery of the global top 1% variants in both affinity and productivity **without exhaustive measurement of the combinatorial library** — only the low-order (≤ double) mutant subset needs to be experimentally screened.

### MATLAB image analysis key features

- **Automated spot detection** — identifies fluorescent spots using adjustable parameters (Gaussian smoothing, adaptive sensitivity, morphological opening, and area filtering).
- **Intensity profiling & oligomeric classification** — extracts mean intensities, performs image-level background subtraction, and calibrates against an internal/external monomeric reference to classify spots into distinct oligomeric states (Monomer, Trimer, Trimer+). In the accompanying manuscript, this was used to confirm that TNF-α bound to Adalimumab and its variants existed predominantly as trimers.
- **NSB correction & quality control** — automatically corrects for Non-Specific Binding (NSB) using user-designated control wells and provides a comprehensive QC summary (saturated spots, specific spot counts, spot-area distributions).
- **Batch export** — generates detailed `.xlsx` reports containing both aggregated per-well statistics and raw intensity distributions for downstream analysis.

---

## Data Availability

`data/demo_library_measurements.csv` contains the processed measurements for the adalimumab HCDR2/HCDR3 combinatorial library analysed in the manuscript, including densMAP coordinates, normalized occupancy and expression metrics with replicate standard deviations, ProteinMPNN compatibility scores, local fitness statistics, and topographical class assignments.

`data/demo_L108G_library_measurements.csv` contains the processed measurements for the L108G-fixed control library (1,910 measured variants on the L108G background), with the same set of metric columns.

`data/demo_image_root/` provides example single-molecule image stacks illustrating the folder structure expected by the MATLAB analysis GUI.

Raw image stacks for the full dataset are available from the corresponding author upon reasonable request.

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
