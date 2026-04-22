# Sequence Embedding and Machine Learning-Guided Antibody Library Analysis

This repository contains the data analysis, image processing, and predictive modeling codes for the manuscript:
> "Landscape-scale navigation unlocks antibody CDR structural logic for AI-guided rescue and therapeutic optimization" (Under Review)

## Overview
This repository provides a complete computational framework for antibody library optimization, spanning from raw single-molecule image quantification to AI-driven sequence design. By integrating single-molecule biophysics with machine learning-guided predictive modeling, we establish a rational strategy to analyze fitness landscapes and discover top-tier variants without the need for exhaustive experimental screening.

The repository is structured into two main components:
1. MATLAB Image Analysis GUI (`/matlab_scripts`): Processes single-molecule fluorescence images to extract intensity profiles and quantify antigen-binding occupancy.
2. Python ML Pipeline (`/python_pipeline`): Performs physicochemical sequence embedding and trains predictive models to navigate the combinatorial fitness landscape.

## Repository Structure
```text
PPI-Landscape-2025/
├── matlab_scripts/        # MATLAB GUI for single-molecule fluorescence image analysis
├── python_pipeline/       # Python scripts for densMAP embedding and ML prediction
├── data/                  # Sample datasets (if available)
├── README.md              # Project overview and instructions
└── requirements.txt       # Python dependencies


### Part 1: Single-Molecule Image Analysis (MATLAB)

The spotAnalysisApp is a batch-processing GUI designed to extract quantitative binding data from .tif fluorescence images (e.g., GFP spots).
Key Features
- Automated Spot Detection: Identifies GFP spots using adjustable parameters (Gaussian smoothing, adaptive sensitivity, morphological opening, and area filtering).
- Intensity Profiling & Oligomeric Classification: Extracts mean intensities, performs image-level background subtraction, and calibrates against an internal/external monomeric reference to classify spots into distinct states (Monomer, Trimer, Trimer+).
- NSB Correction & Quality Control: Automatically corrects for Non-Specific Binding (NSB) using control wells and provides a comprehensive QC summary (e.g., saturated spots, specific spot counts).
- Batch Export: Generates detailed .xlsx reports containing both aggregated statistics and raw intensity distributions for downstream analysis.

### How to Run (MATLAB)
1. Open MATLAB (Requires the Image Processing Toolbox).
2. Navigate to the matlab_scripts/ directory.
3. Run the application by typing spotAnalysisApp_V35_1_Final in the command window.


### Part 2: Sequence Embedding & ML Prediction (Python)
Using the experimental data generated from the MATLAB pipeline, this component optimizes the combinatorial sequence space.

### Key Features
- Biological Preprocessing & Filtering: Curates experimental measurements by strictly enforcing biological constraints (e.g., masking occupancy/affinity values for variants that fail to express).
- Physicochemical Vectorization & densMAP Embedding: Transforms combinatorial mutations into high-dimensional numerical vectors based on position-specific physicochemical properties (Kyte–Doolittle hydropathy, residue volume, and isoelectric point). This sequence space is then projected into a topology-preserving 2D landscape using densMAP (UMAP) and systematically categorized via KMeans clustering.
- ML-Guided Predictive Recovery: Bypasses brute-force screening by training an additive Ridge regression model exclusively on low-order (≤ double) mutants. By effectively capturing fundamental epistatic interactions (pairwise mutations), the model drastically condenses the search space and accurately recovers the global top 1% champions in both affinity and productivity.

### Environment Setup (Python)
We recommend using a virtual environment (e.g., Anaconda or venv). To ensure reproducibility and prevent API drift, install the exact pinned dependencies:
git clone https://github.com/tyyoonlab-snu/PPI-Landscape-2025.git
cd PPI-Landscape-2025
pip install -r requirements.txt

### How to Run (Python)
Ensure your experimental data (.csv or .xlsx) is placed in the data/ folder.

cd python_pipeline
python analysis_pipeline.py --input ../data/your_dataset.csv --outdir output/