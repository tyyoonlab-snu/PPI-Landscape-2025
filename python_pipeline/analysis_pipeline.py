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
    print(f" Loading data from: {filepath}")
    
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

    print(f" Preprocessing complete: {len(df_processed)} unique variants.")
    return df_processed


# Amino-acid features: [Kyte-Doolittle hydropathy, residue volume (Å^3), pI].
# Parental residues are placed at the origin so that each variant is encoded
# by its deviation from the parental sequence.
AA_PROPS = {
    "A": [1.8, 88.6, 6.00], "R": [-4.5, 173.4, 10.76], "N": [-3.5, 114.1, 5.41],
    "D": [-3.5, 111.1, 2.77], "C": [2.5, 108.5, 5.07], "E": [-3.5, 138.4, 3.22],
    "Q": [-3.5, 143.8, 5.65], "G": [-0.4, 60.1, 5.97], "H": [-3.2, 153.2, 7.59],
    "I": [4.5, 166.7, 6.02], "L": [3.8, 166.7, 5.98], "K": [-3.9, 168.6, 9.74],
    "M": [1.9, 162.9, 5.74], "F": [2.8, 189.9, 5.48], "P": [-1.6, 112.7, 6.30],
    "S": [-0.8, 89.0, 5.68], "T": [-0.7, 116.1, 5.60], "W": [-0.9, 227.8, 5.89],
    "Y": [-1.3, 193.6, 5.66], "V": [4.2, 140.0, 5.96], "WT": [0.0, 0.0, 0.00],
}

# Combinatorial design rules: {position: parental residue + permitted mutants}.
# This defines the sequence space of the adalimumab HCDR2/HCDR3 library and
# must match the library being analysed.
DESIGN_RULES = {
    52: {"wt": "T", "muts": ["V"]},
    55: {"wt": "S", "muts": ["H", "G", "A", "W"]},
    57: {"wt": "H", "muts": ["Y", "W"]},
    58: {"wt": "I", "muts": ["H", "D", "T"]},
    99: {"wt": "V", "muts": ["T"]},
    100: {"wt": "S", "muts": ["K", "T", "A", "L"]},
    103: {"wt": "S", "muts": ["P"]},
    106: {"wt": "S", "muts": ["G"]},
    108: {"wt": "L", "muts": ["S"]},
}


def build_feature_matrix(design_rules=DESIGN_RULES, aa_props=AA_PROPS):
    """Enumerate the full combinatorial library and encode each variant.

    Returns the list of variant IDs (underscore-separated mutation codes, with
    the parental sequence labelled "Parental") and the matching physicochemical
    feature matrix, one 3-value block per designed position.
    """
    positions = sorted(design_rules.keys())
    combos = itertools.product(
        *[["WT"] + [f"{design_rules[p]['wt']}{p}{m}" for m in design_rules[p]["muts"]]
          for p in positions]
    )

    ids, vectors = [], []
    for combo in combos:
        active = [m for m in combo if m != "WT"]
        ids.append("Parental" if not active else "_".join(active))
        vec = []
        for pos, item in zip(positions, combo):
            residue = design_rules[pos]["wt"] if item == "WT" else item[-1]
            vec.extend(aa_props[residue])
        vectors.append(vec)

    return ids, np.asarray(vectors, dtype=float)


def perform_umap_clustering(df_processed):
    """densMAP projection of the combinatorial library, with an exploratory
    silhouette-selected KMeans partition, merged onto the experimental data.

    The embedding is built from the physicochemical encoding of the design
    space (`build_feature_matrix`), not from random features. Topographical
    regions used in the manuscript are assigned separately, in
    `topographical_classification.py`.
    """
    print("Running densMAP embedding...")

    ids, features = build_feature_matrix()
    X_scaled = StandardScaler().fit_transform(features)

    reducer = umap.UMAP(n_neighbors=50, min_dist=0.1, spread=2.0,
                        densmap=True, random_state=RANDOM_STATE)
    embedding = reducer.fit_transform(X_scaled)

    df_embed = pd.DataFrame({
        "ID": ids,
        "UMAP1": embedding[:, 0],
        "UMAP2": embedding[:, 1],
    }).merge(df_processed, on="ID", how="left")

    print("Selecting k by silhouette score...")
    coords = df_embed[["UMAP1", "UMAP2"]]
    best_k, best_score = 5, -1.0
    for k in range(3, 11):
        labels = KMeans(n_clusters=k, random_state=RANDOM_STATE,
                        n_init=10).fit_predict(coords)
        score = silhouette_score(coords, labels)
        if score > best_score:
            best_score, best_k = score, k

    df_embed["Cluster"] = KMeans(n_clusters=best_k, random_state=RANDOM_STATE,
                                 n_init=10).fit_predict(coords)

    print(f"UMAP + clustering complete (best k = {best_k}, "
          f"silhouette = {best_score:.3f}).")
    return df_embed


def analyze_predictive_recovery(df_plot, out_dir="output"):
    print(" Running Panel B Predictive Recovery Analysis...")
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
    print(f" Saved figure to {svg_path}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Antibody Library ML Analysis")
    parser.add_argument("--input", type=str, required=True, help="Path to input CSV/Excel data file")
    parser.add_argument("--outdir", type=str, default="output", help="Directory to save outputs")
    args = parser.parse_args()

    os.makedirs(args.outdir, exist_ok=True)

    df = load_and_preprocess_data(args.input)
    df = perform_umap_clustering(df)

    # Save the embedding so the classification and visualization modules can
    # consume it as their `--input`.
    coords_path = os.path.join(args.outdir, "1_UMAP_2D_Coordinates.csv")
    df.to_csv(coords_path, index=False)
    print(f"Saved embedding to {coords_path}")

    analyze_predictive_recovery(df, out_dir=args.outdir)
    print(" All analyses completed successfully.")