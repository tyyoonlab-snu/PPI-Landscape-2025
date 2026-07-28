"""
Topographical classification of the combinatorial fitness landscape.

Assigns every variant in the densMAP embedding to one of four topographical
classes -- Peak Cluster, Rugged Interface, Broad Valley, Normal/Other -- from
local fitness statistics computed over its k nearest neighbours in the 2D
embedding, followed by a density-based contiguity filter.

Pipeline position
-----------------
Consumes the 2D embedding produced by Section 3 (`1_UMAP_2D_Coordinates.csv`)
and reproduces the `Local_Mean_Aff`, `Local_Std_Aff` and
`Topolgy_Classification` columns distributed with
`data/demo_library_measurements.csv`.

Corresponding manuscript items
------------------------------
Fig. 3d-k, Extended Data Fig. 6e-l.
"""

import argparse
import warnings

import numpy as np
import pandas as pd
from sklearn.cluster import DBSCAN
from sklearn.neighbors import NearestNeighbors

warnings.filterwarnings("ignore")

# ==========================================
# [Reproducibility] Fixed random seed
# ==========================================
RANDOM_STATE = 42

# ==========================================
# Classification parameters
# ==========================================
# Number of nearest neighbours (in the 2D densMAP plane) used to compute the
# local mean and local standard deviation of affinity for each variant.
K_NEIGHBORS = 30

# Percentile thresholds applied to the local statistics.
#   Rugged interface : top    RUGGED_PERCENTILE % of local standard deviation
#   Peak cluster     : top   (100 - PEAK_PERCENTILE) % of local mean
#   Broad valley     : bottom VALLEY_PERCENTILE % of local mean
RUGGED_PERCENTILE = 89
PEAK_PERCENTILE = 91
VALLEY_PERCENTILE = 20

# Density-based contiguity filter. Variants that satisfy a threshold but are
# spatially isolated are reassigned to "Normal/Other", so that the reported
# classes correspond to contiguous regions of the landscape.
DBSCAN_EPS = 0.8
DBSCAN_MIN_SAMPLES = 10

# Occupancy floor applied before the log2 transform. Variants with no
# detectable antigen binding are assigned this fractional occupancy rather
# than being discarded, so that non-binding regions still contribute to the
# local statistics of their neighbourhood (log2(0.05) = -4.32).
OCCUPANCY_FLOOR = 0.05

CLASS_PEAK = "Peak Cluster"
CLASS_RUGGED = "Rugged Interface"
CLASS_VALLEY = "Broad Valley"
CLASS_OTHER = "Normal/Other"


def compute_local_statistics(df, coord_cols=("UMAP1", "UMAP2"),
                             value_col="Log_Normalized_Occupancy",
                             k=K_NEIGHBORS):
    """Local mean / standard deviation of `value_col` over the k nearest
    neighbours of each variant in the densMAP plane.

    The query point is included among its own k neighbours, matching the
    convention used to generate the deposited landscape data. Variants whose
    fitness value is missing (non-expressing clones) are skipped when the
    neighbourhood statistics are accumulated but still receive a value.
    """
    coords = df.loc[:, list(coord_cols)].to_numpy(dtype=float)
    values = df[value_col].to_numpy(dtype=float)

    # Variants without embedding coordinates cannot be placed on the landscape
    # and are excluded from both the neighbour search and the neighbourhoods.
    has_coords = np.isfinite(coords).all(axis=1)

    local_mean = np.full(len(df), np.nan)
    local_std = np.full(len(df), np.nan)
    if has_coords.sum() < k:
        return local_mean, local_std

    nn = NearestNeighbors(n_neighbors=k).fit(coords[has_coords])
    _, neighbour_idx = nn.kneighbors(coords[has_coords])

    neighbour_values = values[has_coords][neighbour_idx]
    local_mean[has_coords] = np.nanmean(neighbour_values, axis=1)
    local_std[has_coords] = np.nanstd(neighbour_values, axis=1, ddof=0)

    return local_mean, local_std


def apply_density_filter(df, mask, eps=DBSCAN_EPS,
                         min_samples=DBSCAN_MIN_SAMPLES,
                         coord_cols=("UMAP1", "UMAP2")):
    """Keep only spatially contiguous members of a candidate class.

    DBSCAN is run on the densMAP coordinates of the candidate variants; points
    labelled as noise (-1) are dropped from the class.
    """
    if mask.sum() == 0:
        return mask

    coords = df.loc[mask, list(coord_cols)].to_numpy(dtype=float)
    labels = DBSCAN(eps=eps, min_samples=min_samples).fit_predict(coords)

    retained = mask.copy()
    retained[mask] = labels != -1
    return retained


def classify_topography(df,
                        value_col="Log_Normalized_Occupancy",
                        k=K_NEIGHBORS,
                        rugged_percentile=RUGGED_PERCENTILE,
                        peak_percentile=PEAK_PERCENTILE,
                        valley_percentile=VALLEY_PERCENTILE,
                        eps=DBSCAN_EPS,
                        min_samples=DBSCAN_MIN_SAMPLES):
    """Assign topographical classes to every variant.

    Precedence is Rugged Interface -> Peak Cluster / Broad Valley -> Other:
    a variant sitting in a region of high local variability is reported as a
    rugged interface even if its local mean is extreme.
    """
    df = df.copy()

    local_mean, local_std = compute_local_statistics(
        df, value_col=value_col, k=k
    )
    df["Local_Mean_Aff"] = local_mean
    df["Local_Std_Aff"] = local_std

    std_cut = np.nanpercentile(local_std, rugged_percentile)
    peak_cut = np.nanpercentile(local_mean, peak_percentile)
    valley_cut = np.nanpercentile(local_mean, valley_percentile)

    print(f"   Rugged interface : Local_Std_Aff  >= {std_cut:.3f} "
          f"({rugged_percentile}th pct)")
    print(f"   Peak cluster     : Local_Mean_Aff >= {peak_cut:.3f} "
          f"({peak_percentile}th pct)")
    print(f"   Broad valley     : Local_Mean_Aff <= {valley_cut:.3f} "
          f"({valley_percentile}th pct)")

    scored = np.isfinite(local_mean) & np.isfinite(local_std)
    is_rugged = scored & (local_std >= std_cut)
    is_peak = scored & (~is_rugged) & (local_mean >= peak_cut)
    is_valley = scored & (~is_rugged) & (local_mean <= valley_cut)

    # Retain only contiguous regions within each candidate class.
    is_rugged = apply_density_filter(df, is_rugged, eps, min_samples)
    is_peak = apply_density_filter(df, is_peak, eps, min_samples)
    is_valley = apply_density_filter(df, is_valley, eps, min_samples)

    classification = pd.Series(CLASS_OTHER, index=df.index, dtype=object)
    classification[~scored] = np.nan
    classification[is_valley] = CLASS_VALLEY
    classification[is_peak] = CLASS_PEAK
    classification[is_rugged] = CLASS_RUGGED

    df["Topolgy_Classification"] = classification
    return df


def summarize(df):
    """Print class sizes and the fitness range spanned by each class."""
    summary = (
        df.groupby("Topolgy_Classification")[["Local_Mean_Aff", "Local_Std_Aff"]]
        .agg(["count", "min", "max", "mean"])
        .round(3)
    )
    print("\n=== Topographical class summary ===")
    print(summary)
    return summary


def validate_against_reference(df, reference_path):
    """Optional consistency check against deposited landscape data.

    Matches rows on whichever identifier column both tables share
    ("Mutation_Description" or "ID"), so the check works whether the input
    came from the notebook or from analysis_pipeline.py.
    """
    ref = pd.read_csv(reference_path)
    if "Topolgy_Classification" not in ref.columns:
        print("Reference file carries no classification column; skipping.")
        return None

    id_col = next((c for c in ("Mutation_Description", "ID")
                   if c in df.columns and c in ref.columns), None)
    if id_col is None:
        print("No shared identifier column with the reference; skipping check.")
        return None

    merged = df.merge(
        ref[[id_col, "Topolgy_Classification"]],
        on=id_col, how="inner", suffixes=("", "_ref"),
    )
    agreement = (
        merged["Topolgy_Classification"] == merged["Topolgy_Classification_ref"]
    ).mean()
    print(f"\nAgreement with reference classification: {agreement * 100:.1f}% "
          f"(n = {len(merged)})")
    return agreement


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Topographical classification of the fitness landscape"
    )
    parser.add_argument("--input", required=True,
                        help="CSV with UMAP1, UMAP2 and the fitness column")
    parser.add_argument("--value-col", default="Log_Normalized_Occupancy",
                        help="Fitness column used for local statistics")
    parser.add_argument("--output", default="2_Topographical_Classification.csv")
    parser.add_argument("--reference", default=None,
                        help="Optional reference CSV for a consistency check")
    args = parser.parse_args()

    data = pd.read_csv(args.input)

    if args.value_col not in data.columns and "Normalized_Occupancy" in data.columns:
        print(f"'{args.value_col}' not found; deriving it from "
              f"'Normalized_Occupancy' (log2, floored at {OCCUPANCY_FLOOR}).")
        data[args.value_col] = np.log2(
            data["Normalized_Occupancy"].replace(0, OCCUPANCY_FLOOR)
        )

    print(f"Classifying {len(data)} variants "
          f"(k = {K_NEIGHBORS}, seed = {RANDOM_STATE})...")
    result = classify_topography(data, value_col=args.value_col)
    summarize(result)

    if args.reference:
        validate_against_reference(result, args.reference)

    result.to_csv(args.output, index=False)
    print(f"\nSaved classification to {args.output}")
