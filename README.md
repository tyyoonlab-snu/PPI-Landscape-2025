# Sequence Embedding and Machine Learning-Guided Antibody Library Analysis

This repository contains the data analysis and predictive modeling code for the manuscript:
> Landscape-scale navigation unlocks antibody CDR structural logic for AI-guided rescue and therapeutic optimization (Under Review)

## Overview
This pipeline demonstrates a comprehensive computational framework to efficiently navigate and optimize a combinatorial antibody sequence space. By integrating physicochemical property-based sequence vectorization with machine learning-guided predictive modeling, we establish a rational strategy to analyze fitness landscapes and discover top-tier variants without the need for exhaustive experimental screening.

### Key Features
1. Biological Preprocessing & Filtering: Curates experimental measurements by strictly enforcing biological constraints (e.g., masking occupancy/affinity values for variants that fail to express).
2. Physicochemical Vectorization & densMAP Embedding: Transforms combinatorial mutations into high-dimensional numerical vectors based on position-specific physicochemical properties (Kyte–Doolittle hydropathy, residue volume, and isoelectric point). This sequence space is then projected into a topology-preserving 2D landscape using 'densMAP' (UMAP) and systematically categorized via KMeans clustering.
3. ML-Guided Predictive Recovery: Bypasses brute-force screening by training an additive Ridge regression model exclusively on low-order (≤ double) mutants. By effectively capturing fundamental epistatic interactions (pairwise mutations), the model drastically condenses the search space and accurately recovers the global top 1% champions in both affinity and productivity.

## Environment Setup
We recommend using a virtual environment (e.g., Anaconda or 'venv'). To ensure reproducibility and prevent API drift, install the exact pinned dependencies:

```bash
git clone https://github.com/tyyoonlab-snu/PPI-Landscape-2025.git
cd PPI-Landscape-2025
pip install -r requirements.txt