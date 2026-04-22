"""
Code Availability: ML-Guided Antibody Library Design & Analysis
This script performs data preprocessing, UMAP embedding, KMeans clustering,
and predictive recovery analysis (Ridge regression) of antibody variants.
"""

import os
import re
import random
import warnings
import itertools
import argparse

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

import umap
from sklearn.preprocessing import StandardScaler, MultiLabelBinarizer
from sklearn.cluster import KMeans
from sklearn.linear_model import Ridge
from sklearn.metrics import silhouette_score

warnings.filterwarnings("ignore")

# ==========================================
# [Reproducibility] Fix global random seeds
# ==========================================
RANDOM_STATE = 42
random.seed(RANDOM_STATE)
np.random.seed(RANDOM_STATE)


def load_and_preprocess_data(filepath):
    print(f"📂 Loading data from: {filepath}")
    
    if filepath.endswith(".csv"):
        df_raw = pd.read_csv(filepath)
    else:
        df_raw = pd.read_excel(filepath)

    df_raw.columns = df_raw.columns.str.strip()
    target_id_col = "Mutation_Description"
    if target_id_col not in df_raw.columns:
        df_raw.rename(columns={df_raw.columns[0]: target_id_col}, inplace=True)
    df_raw[target_id_col] = df_raw[target_id_col].astype(str).str.strip()

    # Data Quality filtering
    VALID_QUALITY_LABELS = ["Valid", "Valid (Non-binding)", "Valid; Valid (Non-binding)"]
    EXPRESSION_ONLY_LABELS = ["Invalid (Low expression)"]
    quality_col = next((col for col in df_raw.columns if "data_quality" in col.lower()), None)

    if quality_col:
        df_raw[quality_col] = df_raw[quality_col].astype(str).str.strip()
        df_raw = df_raw[df_raw[quality_col].isin(VALID_QUALITY_LABELS + EXPRESSION_ONLY_LABELS)].copy()

    # Numeric coercion & Biological filter
    targets = [c for c in ["Normalized_Occupancy", "Normalized_Expression"] if c in df_raw.columns]
    for col in targets:
        df_raw[col] = pd.to_numeric(df_raw[col], errors="coerce")

    if "Normalized_Occupancy" in targets:
        if quality_col:
            df_raw.loc[df_raw[quality_col].isin(EXPRESSION_ONLY_LABELS), "Normalized_Occupancy"] = np.nan
        if "Normalized_Expression" in targets:
            df_raw.loc[df_raw["Normalized_Expression"].isna(), "Normalized_Occupancy"] = np.nan

    # Aggregate replicates
    df_agg = df_raw.groupby(target_id_col)[targets].agg(["mean", "std"])
    df_agg.columns = [f"{a}_{b}" for a, b in df_agg.columns]
    df_processed = df_agg.reset_index().rename(columns={target_id_col: "ID"})

    for col in targets:
        mean_col = f"{col}_mean"
        if mean_col in df_processed.columns:
            df_processed[f"Log_{col}"] = np.log2(df_processed[mean_col].clip(lower=0.05))

    print(f"✅ Preprocessing complete: {len(df_processed)} unique variants.")
    return df_processed


def perform_umap_clustering(df_processed):
    print("🏃 Running densMAP and KMeans clustering...")
    
    # Simple feature extraction (Example mapping, adjust aa_props as in your original code)
    # [Note: Full aa_props and design_rules dicts are abbreviated here for readability]
    aa_props = {"WT": [0.0, 0.0, 0.0], "V": [4.2, 140.0, 5.96], "Y": [-1.3, 193.6, 5.66]} 
    
    # For a general GitHub repo, you might construct features directly from the ID column
    # to avoid hardcoding design rules, but here we use a placeholder random feature 
    # to maintain the pipeline structure if rules aren't strictly hardcoded.
    # In your real upload, paste the FULL 'aa_props' and 'design_rules' here.
    
    X_full = np.random.rand(len(df_processed), 27) # Placeholder for actual features
    X_scaled = StandardScaler().fit_transform(X_full)

    reducer = umap.UMAP(n_neighbors=50, min_dist=0.1, spread=2.0, densmap=True, random_state=RANDOM_STATE)
    embedding = reducer.fit_transform(X_scaled)

    df_processed["UMAP1"] = embedding[:, 0]
    df_processed["UMAP2"] = embedding[:, 1]

    umap_coords = df_processed[["UMAP1", "UMAP2"]]
    kmeans = KMeans(n_clusters=5, random_state=RANDOM_STATE, n_init=10)
    df_processed["Cluster"] = kmeans.fit_predict(umap_coords)

    print("✅ UMAP + Clustering complete.")
    return df_processed


def analyze_predictive_recovery(df_plot, out_dir="output"):
    print("📈 Running Panel B Predictive Recovery Analysis...")
    os.makedirs(out_dir, exist_ok=True)

    def get_mut_list(m_str):
        if pd.isna(m_str) or str(m_str).strip().upper() in ("WT", "PARENTAL", "NONE", ""): return []
        return [m.strip() for m in re.split(r"[,\s_+/]+", str(m_str)) if m.strip()]

    df_plot["mut_list"] = df_plot["ID"].apply(get_mut_list)
    df_plot["mut_count"] = df_plot["mut_list"].apply(len)

    mlb = MultiLabelBinarizer()
    X_all_sparse = mlb.fit_transform(df_plot["mut_list"])

    results_list = []
    target_top_pct = 1
    ratios = np.arange(0.01, 1.01, 0.01)

    for target_col, metric_name in [("Log_Normalized_Occupancy", "Affinity"), ("Log_Normalized_Expression", "Expression")]:
        if target_col not in df_plot.columns: continue

        valid_idx = df_plot[target_col].notna()
        df_valid, X_valid, y_valid = df_plot[valid_idx].copy(), X_all_sparse[valid_idx], df_plot.loc[valid_idx, target_col].values
        
        global_thresh = df_valid[target_col].quantile(0.99)
        actual_top_indices = set(df_valid[df_valid[target_col] >= global_thresh].index)
        
        # Train on <= Double mutants
        train_mask = df_valid["mut_count"] <= 2
        model = Ridge(alpha=1.0).fit(X_valid[train_mask], y_valid[train_mask])
        df_valid["pred"] = model.predict(X_valid)
        
        df_ranked = df_valid.sort_values("pred", ascending=False)
        
        for r in ratios:
            top_n = max(1, int(round(len(df_ranked) * r)))
            recovered = len(set(df_ranked.head(top_n).index).intersection(actual_top_indices))
            results_list.append({"Metric": metric_name, "Ratio": r * 100, "Recovery_Rate": (recovered / len(actual_top_indices) * 100)})

    # Plotting code
    df_res = pd.DataFrame(results_list)
    plt.figure(figsize=(8, 6))
    sns.lineplot(data=df_res, x="Ratio", y="Recovery_Rate", hue="Metric", linewidth=3)
    plt.axvline(10, color="gray", linestyle=":"); plt.axhline(80, color="gray", linestyle=":")
    plt.title(f"Predictive Recovery from ≤Double Mutants")
    
    svg_path = os.path.join(out_dir, "PanelB_Recovery.svg")
    plt.savefig(svg_path, format="svg")
    print(f"✅ Saved figure to {svg_path}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Antibody Library ML Analysis")
    parser.add_argument("--input", type=str, required=True, help="Path to input CSV/Excel data file")
    parser.add_argument("--outdir", type=str, default="output", help="Directory to save outputs")
    args = parser.parse_args()

    df = load_and_preprocess_data(args.input)
    df = perform_umap_clustering(df)
    analyze_predictive_recovery(df, out_dir=args.outdir)
    print("🎉 All analyses completed successfully.")