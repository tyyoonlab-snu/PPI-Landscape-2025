"""
L108G-fixed control landscape.

Constructs a separate combinatorial landscape in which position 108 is held at
the deleterious L108G substitution, projects the measured L108G-background
variants onto it, and compares this control library against the main
functional-variant library.

This library tests whether pre-filtering the combinatorial library to
individually functional substitutions is necessary: including a non-functional
anchor (L108G) is expected to expand the low-fitness "void" regions and reduce
the probability of positive epistasis.

Pipeline position
-----------------
Independent of the main pipeline. Takes the L108G measurement table as input;
the macroscopic comparison additionally takes the 9,600-variant main library
table.

Corresponding manuscript items
------------------------------
Fig. 3p-s, Extended Data Fig. 7.
"""

import argparse
import itertools
import os
import re
import warnings

import matplotlib.colors as mcolors
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.ndimage import gaussian_filter
from sklearn.neighbors import KNeighborsRegressor
from sklearn.preprocessing import StandardScaler

try:
    import umap
except ImportError:  # embedding is optional if precomputed coordinates exist
    umap = None

warnings.filterwarnings("ignore")

plt.rcParams["svg.fonttype"] = "none"
plt.rcParams["pdf.fonttype"] = 42
plt.rcParams["ps.fonttype"] = 42

RANDOM_STATE = 42
OCCUPANCY_FLOOR = 0.05          # matches the main pipeline
GRID_RESOLUTION = 150
GRID_MARGIN = 1.5
KNN_NEIGHBORS = 15
GAUSSIAN_SIGMA = 1.2

# Amino-acid features: [Kyte-Doolittle hydropathy, residue volume (A^3), pI].
AA_PROPS = {
    "A": [1.8, 88.6, 6.00], "R": [-4.5, 173.4, 10.76], "N": [-3.5, 114.1, 5.41],
    "D": [-3.5, 111.1, 2.77], "C": [2.5, 108.5, 5.07], "E": [-3.5, 138.4, 3.22],
    "Q": [-3.5, 143.8, 5.65], "G": [-0.4, 60.1, 5.97], "H": [-3.2, 153.2, 7.59],
    "I": [4.5, 166.7, 6.02], "L": [3.8, 166.7, 5.98], "K": [-3.9, 168.6, 9.74],
    "M": [1.9, 162.9, 5.74], "F": [2.8, 189.9, 5.48], "P": [-1.6, 112.7, 6.30],
    "S": [-0.8, 89.0, 5.68], "T": [-0.7, 116.1, 5.60], "W": [-0.9, 227.8, 5.89],
    "Y": [-1.3, 193.6, 5.66], "V": [4.2, 140.0, 5.96], "WT": [0.0, 0.0, 0.00],
}

# Design rules for the L108G-fixed library. Position 108 is fixed to G; the
# remaining positions carry the functional substitutions used in the manuscript
# (S55W is excluded relative to the main library, giving 1,920 combinations).
DESIGN_RULES = {
    52: {"wt": "T", "muts": ["V"]},
    55: {"wt": "S", "muts": ["H", "G", "A"]},
    57: {"wt": "H", "muts": ["Y", "W"]},
    58: {"wt": "I", "muts": ["H", "D", "T"]},
    99: {"wt": "V", "muts": ["T"]},
    100: {"wt": "S", "muts": ["K", "T", "A", "L"]},
    103: {"wt": "S", "muts": ["P"]},
    108: {"wt": "L", "muts": ["G"]},   # fixed to the L108G background
}

ID_CANDIDATES = ["Mutation_Description", "ID", "Oligo Name", "Variant",
                 "Sequence", "Name"]


def find_id_column(df):
    for col in ID_CANDIDATES:
        if col in df.columns:
            return col
    for col in df.columns:
        if df[col].dtype == object and df[col].astype(str).str.contains("_").any():
            return col
    return df.columns[0]


def build_base_manifold():
    """Enumerate the 1,920 L108G-fixed combinations and embed them.

    Position 108 is fixed to G, so every variant lies on the L108G background.
    The parental sequence "L108G" (no additional substitutions) is the baseline.
    """
    if umap is None:
        raise ImportError("umap-learn is required to build the base manifold.")

    positions = sorted(DESIGN_RULES.keys())
    lists = []
    for p in positions:
        if p == 108:
            lists.append(["L108G"])   # fixed
        else:
            lists.append(["WT"] + [f"{DESIGN_RULES[p]['wt']}{p}{m}"
                                   for m in DESIGN_RULES[p]["muts"]])

    ids, vectors = [], []
    for combo in itertools.product(*lists):
        active = [m for m in combo if m != "WT"]
        ids.append("_".join(active))   # baseline reduces to "L108G"
        vec = []
        for pos, item in zip(positions, combo):
            if pos == 108:
                vec.extend(AA_PROPS["G"])
            else:
                residue = DESIGN_RULES[pos]["wt"] if item == "WT" else item[-1]
                vec.extend(AA_PROPS[residue])
        vectors.append(vec)

    scaler = StandardScaler().fit(np.asarray(vectors, dtype=float))
    reducer = umap.UMAP(n_neighbors=30, min_dist=0.4, n_components=2,
                        metric="euclidean", random_state=RANDOM_STATE)
    embedding = reducer.fit_transform(scaler.transform(
        np.asarray(vectors, dtype=float)))

    base = pd.DataFrame({"ID": ids,
                         "UMAP1_Base": embedding[:, 0],
                         "UMAP2_Base": embedding[:, 1]})
    base["Match_ID"] = base["ID"].str.upper()
    return base, scaler, reducer, positions


def project_measurements(df, scaler, reducer, positions):
    """Encode measured L108G-background variants and project them onto the
    fixed manifold, using the same scaler and reducer as the base library."""
    id_col = find_id_column(df)
    df = df.copy()
    df["ID"] = df[id_col].astype(str)
    df["Match_ID"] = df["ID"].str.upper()

    vectors = []
    for seq in df["ID"]:
        muts = {int(re.findall(r"\d+", m)[0]): m[-1]
                for m in str(seq).upper().split("_") if re.findall(r"\d+", m)}
        vec = []
        for pos in positions:
            if pos == 108:
                vec.extend(AA_PROPS["G"])
            else:
                vec.extend(AA_PROPS.get(muts.get(pos, DESIGN_RULES[pos]["wt"]),
                                        [0, 0, 0]))
        vectors.append(vec)

    coords = reducer.transform(scaler.transform(np.asarray(vectors, dtype=float)))
    df["UMAP1_Proj"] = coords[:, 0]
    df["UMAP2_Proj"] = coords[:, 1]

    # Additional substitutions beyond the fixed L108G background.
    df["Additional_Mutations"] = df["ID"].apply(
        lambda s: len([m for m in str(s).upper().split("_") if m and m != "L108G"])
    )
    return df


def ensure_log_column(df, keyword, fallback):
    """Return the name of a log2 fitness column, creating it (floored) if
    only the raw column is present."""
    for col in df.columns:
        if col.startswith("Log_") and keyword in col.lower():
            return col
    raw = next((c for c in df.columns
                if keyword in c.lower() and not c.startswith("Log_")
                and not c.lower().endswith("_std")), None)
    if raw is None:
        return None
    log_col = f"Log_{raw}"
    df[log_col] = np.log2(pd.to_numeric(df[raw], errors="coerce")
                          .clip(lower=OCCUPANCY_FLOOR))
    return log_col


def render_surface(df, log_col, base, out_dir, save_format="svg"):
    """Interpolated 2D top-view of a projected metric, baselined at L108G."""
    scored = df.dropna(subset=[log_col, "UMAP1_Proj", "UMAP2_Proj"])
    if len(scored) < 10:
        print(f"   Skipping '{log_col}': too few points.")
        return
    short = "Affinity" if "occupancy" in log_col.lower() else "Expression"

    x = np.linspace(base["UMAP1_Base"].min() - GRID_MARGIN,
                    base["UMAP1_Base"].max() + GRID_MARGIN, GRID_RESOLUTION)
    y = np.linspace(base["UMAP2_Base"].min() - GRID_MARGIN,
                    base["UMAP2_Base"].max() + GRID_MARGIN, GRID_RESOLUTION)
    grid_x, grid_y = np.meshgrid(x, y)

    knn = KNeighborsRegressor(n_neighbors=KNN_NEIGHBORS, weights="distance").fit(
        scored[["UMAP1_Proj", "UMAP2_Proj"]].to_numpy(float),
        scored[log_col].to_numpy(float))
    surface = gaussian_filter(
        knn.predict(np.c_[grid_x.ravel(), grid_y.ravel()]).reshape(grid_x.shape),
        sigma=GAUSSIAN_SIGMA)

    limit = max(abs(np.nanmin(surface)), abs(np.nanmax(surface)), 1.0)
    cmap = (plt.cm.RdBu_r if "occupancy" in log_col.lower()
            else mcolors.LinearSegmentedColormap.from_list(
                "TealOrange", ["#008080", "#FFFFFF", "#FF8C00"]))

    baseline = base[base["Match_ID"] == "L108G"]
    bx, by = (baseline[["UMAP1_Base", "UMAP2_Base"]].to_numpy()[0]
              if not baseline.empty else (0.0, 0.0))

    fig, ax = plt.subplots(figsize=(8, 7))
    filled = ax.contourf(grid_x, grid_y, surface, levels=100, cmap=cmap,
                         norm=mcolors.TwoSlopeNorm(vmin=-limit, vcenter=0, vmax=limit),
                         extend="both")
    ax.contour(grid_x, grid_y, surface, levels=15, colors="black",
               linewidths=0.3, alpha=0.3)
    ax.scatter(bx, by, s=300, c="gold", marker="*", edgecolors="black",
               linewidths=1.5, zorder=10)
    ax.text(bx + 0.3, by + 0.3, "L108G", fontsize=11, fontweight="bold",
            bbox=dict(facecolor="white", alpha=0.7, edgecolor="none", pad=1))
    plt.colorbar(filled, ax=ax, fraction=0.046, pad=0.04).set_label("Log2 score")
    ax.set_title(f"L108G-fixed {short} landscape (N={len(scored):,})")
    ax.set_xlabel("UMAP1"); ax.set_ylabel("UMAP2"); ax.set_aspect("equal")

    os.makedirs(out_dir, exist_ok=True)
    matrix_path = os.path.join(out_dir, f"SourceData_L108G_{short}.csv")
    pd.DataFrame(surface, index=np.round(y, 4),
                 columns=np.round(x, 4)).to_csv(matrix_path)
    fig_path = os.path.join(out_dir, f"L108G_{short}_landscape.{save_format}")
    plt.savefig(fig_path, dpi=300, bbox_inches="tight", format=save_format)
    plt.close()
    print(f"   Rendered L108G {short} landscape ({fig_path}); "
          f"source data {matrix_path}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="L108G-fixed control landscape analysis"
    )
    parser.add_argument("--input", required=True,
                        help="L108G measurement table (CSV or XLSX)")
    parser.add_argument("--outdir", default="output_l108g")
    parser.add_argument("--format", default="svg", choices=["svg", "png"])
    args = parser.parse_args()

    data = (pd.read_excel(args.input) if args.input.lower().endswith(("xlsx", "xls"))
            else pd.read_csv(args.input))
    print(f"Loaded {len(data)} L108G-background variants.")

    base, scaler, reducer, positions = build_base_manifold()
    print(f"Built L108G-fixed base manifold ({len(base)} combinations).")

    data = project_measurements(data, scaler, reducer, positions)

    aff_col = ensure_log_column(data, "occupancy", "Log_Normalized_Occupancy")
    exp_col = ensure_log_column(data, "expression", "Log_Normalized_Expression")

    for col in (aff_col, exp_col):
        if col:
            render_surface(data, col, base, args.outdir, save_format=args.format)

    os.makedirs(args.outdir, exist_ok=True)
    out_csv = os.path.join(args.outdir, "L108G_projected_coordinates.csv")
    keep = ["ID", "UMAP1_Proj", "UMAP2_Proj", "Additional_Mutations"] + \
           [c for c in (aff_col, exp_col) if c]
    data[keep].to_csv(out_csv, index=False)
    print(f"Saved projected coordinates to {out_csv}")
    print("Done.")
