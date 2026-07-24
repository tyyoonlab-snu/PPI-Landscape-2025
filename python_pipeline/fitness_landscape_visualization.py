"""
Continuous fitness landscape rendering (3D surface and 2D top-view).

Builds the continuous surfaces shown in the manuscript from the discrete
densMAP embedding. Measured fitness values are interpolated onto a regular
grid by distance-weighted k-nearest-neighbour regression and smoothed, then
expressed in standard-deviation units relative to the parental clone.

Interpolation is deliberately performed over the full rectangular domain
without convex-hull or density masking, so regions of the rendered surface
that fall outside the sampled sequence space are inferred rather than
measured. All quantitative statements about landscape ruggedness and
navigability in the manuscript are derived from the discrete measured
variants, not from these surfaces.

Pipeline position
-----------------
Consumes the merged table produced by Section 3 (densMAP coordinates plus
experimental metrics).

Corresponding manuscript items
------------------------------
Fig. 3c-i, Fig. 4a,b, Extended Data Fig. 8a,b.
"""

import argparse
import os
import warnings

import matplotlib.colors as mcolors
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import plotly.graph_objects as go
from scipy.ndimage import gaussian_filter
from sklearn.neighbors import KNeighborsRegressor

warnings.filterwarnings("ignore")

# Keep text editable in vector output (Nature Portfolio figure requirement).
plt.rcParams["svg.fonttype"] = "none"
plt.rcParams["pdf.fonttype"] = 42
plt.rcParams["ps.fonttype"] = 42

# ==========================================
# Surface construction parameters
# ==========================================
GRID_RESOLUTION = 150      # points per axis
GRID_MARGIN = 1.5          # densMAP units added beyond the data range
KNN_NEIGHBORS = 15         # neighbours used for distance-weighted interpolation
GAUSSIAN_SIGMA = 1.2       # smoothing width, in grid units

# Colour-scale limits, in standard deviations from the parental clone.
EXPERIMENTAL_CLIP = 3.0    # symmetric cap for affinity / expression metrics
MPNN_LOWER_BOUND = -6.0    # lower cap for ProteinMPNN compatibility scores
MPNN_UPPER_BOUND = 1.0


def _is_mpnn(target_name):
    return "MPNN" in target_name or "Z-Score" in target_name


def get_colormap(target_name):
    """Matplotlib colormap matching the manuscript colour conventions."""
    if _is_mpnn(target_name):
        return plt.cm.PRGn
    if "Expression" in target_name:
        return mcolors.LinearSegmentedColormap.from_list(
            "TealOrange", ["#008080", "#FFFFFF", "#FF8C00"]
        )
    return plt.cm.RdBu_r


def get_plotly_colorscale(target_name):
    if _is_mpnn(target_name):
        return "PRGn"
    if "Expression" in target_name:
        return [[0.0, "#008080"], [0.5, "#FFFFFF"], [1.0, "#FF8C00"]]
    return "RdBu_r"


def build_surface(df, target_col, is_mpnn=False):
    """Interpolate a fitness metric onto a regular grid over the embedding.

    Returns the grid coordinates, the surface expressed in standard-deviation
    units relative to the parental clone, and the parental coordinates.
    """
    scored = df.dropna(subset=[target_col, "UMAP1", "UMAP2"]).copy()
    values = scored[target_col].to_numpy(dtype=float)
    if is_mpnn:
        # Sign is inverted so that positive values denote improved
        # sequence-structure compatibility.
        values = -values

    x = np.linspace(df["UMAP1"].min() - GRID_MARGIN,
                    df["UMAP1"].max() + GRID_MARGIN, GRID_RESOLUTION)
    y = np.linspace(df["UMAP2"].min() - GRID_MARGIN,
                    df["UMAP2"].max() + GRID_MARGIN, GRID_RESOLUTION)
    grid_x, grid_y = np.meshgrid(x, y)

    knn = KNeighborsRegressor(
        n_neighbors=KNN_NEIGHBORS, weights="distance"
    ).fit(scored[["UMAP1", "UMAP2"]].to_numpy(dtype=float), values)

    surface = knn.predict(np.c_[grid_x.ravel(), grid_y.ravel()])
    surface = gaussian_filter(surface.reshape(grid_x.shape), sigma=GAUSSIAN_SIGMA)

    # Express the surface in standard deviations relative to the parental clone.
    scale = np.nanstd(surface)
    if scale < 1e-6:
        scale = 1.0

    id_col = next((c for c in ("ID", "Mutation_Description") if c in df.columns),
                  None)
    parental = (df[df[id_col].astype(str).str.contains("Parental|WT", case=False,
                                                       na=False)]
                if id_col else df.iloc[0:0])
    parental_xy = None
    if not parental.empty:
        parental_value = parental[target_col].fillna(0.0).to_numpy(dtype=float)[0]
        if is_mpnn:
            parental_value = -parental_value
        surface_sd = (surface - parental_value) / scale
        parental_xy = (parental["UMAP1"].to_numpy()[0],
                       parental["UMAP2"].to_numpy()[0])
    else:
        surface_sd = (surface - np.nanmean(surface)) / scale

    return grid_x, grid_y, surface_sd, parental_xy, x, y


def colour_limits(surface_sd, is_mpnn=False):
    """Colour-scale limits, capped so that outliers do not dominate."""
    if is_mpnn:
        lower = max(np.floor(np.nanpercentile(surface_sd, 2)), MPNN_LOWER_BOUND)
        return float(lower), MPNN_UPPER_BOUND

    robust = np.nanpercentile(np.abs(surface_sd), 98)
    upper = EXPERIMENTAL_CLIP if robust >= EXPERIMENTAL_CLIP else float(np.ceil(robust))
    return -upper, upper


def render_target(df, target_name, out_dir="output", save_format="svg", ax=None):
    """Render one metric as an interactive 3D surface and a 2D top-view."""
    target_col = (f"Log_{target_name}"
                  if f"Log_{target_name}" in df.columns else target_name)
    if target_col not in df.columns:
        print(f"   Skipping '{target_name}': column not found.")
        return

    os.makedirs(out_dir, exist_ok=True)
    is_mpnn = _is_mpnn(target_name)

    grid_x, grid_y, surface_sd, parental_xy, x, y = build_surface(
        df, target_col, is_mpnn
    )
    c_min, c_max = colour_limits(surface_sd, is_mpnn)
    ticks = np.unique(
        np.linspace(int(np.ceil(c_min)), int(np.floor(c_max)), 7).round().astype(int)
    )

    # Source Data: the interpolated matrix underlying the rendered surface.
    matrix_path = os.path.join(out_dir, f"SourceData_Surface_{target_name}.csv")
    pd.DataFrame(surface_sd, index=np.round(y, 4),
                 columns=np.round(x, 4)).to_csv(matrix_path)

    # ---- 3D surface (interactive) ----
    figure = go.Figure(data=[go.Surface(
        z=surface_sd, x=grid_x, y=grid_y,
        colorscale=get_plotly_colorscale(target_name),
        cmid=0, cmin=c_min, cmax=c_max, opacity=0.9,
        colorbar=dict(title="SD from parental", thickness=20, len=0.5, x=1.1,
                      tickmode="array", tickvals=ticks,
                      ticktext=[str(v) for v in ticks]),
    )])
    if parental_xy is not None:
        figure.update_layout(scene=dict(annotations=[dict(
            x=parental_xy[0], y=parental_xy[1], z=0.0, text="Parental",
            font=dict(color="black", size=13), showarrow=False, yshift=10,
        )]))
    figure.update_layout(
        title=dict(text=f"<b>{target_name}</b><br>"
                        f"<sup>Continuous interpolated surface</sup>", x=0.5),
        width=1000, height=800,
        scene=dict(xaxis_title="UMAP1", yaxis_title="UMAP2",
                   zaxis=dict(title="SD from parental",
                              range=[np.nanmin(surface_sd) - 0.5,
                                     np.nanmax(surface_sd) + 0.5]),
                   aspectratio=dict(x=1, y=1, z=0.85)),
    )
    figure.write_html(os.path.join(out_dir, f"Landscape3D_{target_name}.html"),
                      include_plotlyjs="cdn")

    # ---- 2D top-view ----
    standalone = ax is None
    if standalone:
        _, ax = plt.subplots(figsize=(8, 7))

    norm = (mcolors.TwoSlopeNorm(vmin=c_min, vcenter=0, vmax=c_max)
            if c_min < 0 < c_max and abs(c_min) != abs(c_max)
            else mcolors.Normalize(vmin=c_min, vmax=c_max))

    filled = ax.contourf(grid_x, grid_y, surface_sd,
                         levels=np.linspace(c_min, c_max, 100),
                         cmap=get_colormap(target_name), norm=norm, extend="both")
    ax.contour(grid_x, grid_y, surface_sd,
               levels=np.linspace(c_min, c_max, 15),
               colors="black", linewidths=0.5, alpha=0.5)

    if parental_xy is not None:
        ax.scatter(*parental_xy, s=150, c="gold", marker="D",
                   edgecolors="black", linewidths=2, zorder=10)
        ax.text(parental_xy[0], parental_xy[1] + 0.3, "Parental", fontsize=11,
                ha="center", va="bottom",
                bbox=dict(facecolor="white", alpha=0.5, edgecolor="none", pad=1))

    bar = plt.colorbar(filled, ax=ax, ticks=ticks, fraction=0.046, pad=0.04)
    bar.set_label("Standard deviations from parental", fontsize=11)

    ax.set_title(target_name, fontsize=13, pad=15)
    ax.set_xlabel("UMAP1", fontsize=11)
    ax.set_ylabel("UMAP2", fontsize=11)
    ax.set_xlim(grid_x.min(), grid_x.max())
    ax.set_ylim(grid_y.min(), grid_y.max())
    ax.set_aspect("equal")

    if standalone:
        plt.tight_layout()
        plt.savefig(os.path.join(out_dir, f"TopView_{target_name}.{save_format}"),
                    dpi=300, bbox_inches="tight", format=save_format)
        plt.close()

    print(f"   Rendered '{target_name}' (source data: {matrix_path})")


def render_all(df, targets, out_dir="output", save_format="svg"):
    """Render every requested metric and a combined top-view panel."""
    os.makedirs(out_dir, exist_ok=True)
    figure, axes = plt.subplots(1, len(targets), figsize=(7 * len(targets), 6))
    if len(targets) == 1:
        axes = [axes]

    for ax, target in zip(axes, targets):
        render_target(df, target, out_dir=out_dir,
                      save_format=save_format, ax=ax)

    plt.tight_layout()
    combined = os.path.join(out_dir, f"TopView_Combined.{save_format}")
    plt.savefig(combined, dpi=300, bbox_inches="tight", format=save_format)
    plt.close()
    print(f"Saved combined top-view to {combined}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Render continuous fitness landscapes"
    )
    parser.add_argument("--input", required=True,
                        help="CSV with UMAP1, UMAP2, ID and metric columns")
    parser.add_argument("--targets", nargs="+",
                        default=["Normalized_Occupancy", "Normalized_Expression"],
                        help="Metric columns to render")
    parser.add_argument("--outdir", default="output")
    parser.add_argument("--format", default="svg", choices=["svg", "png"])
    args = parser.parse_args()

    data = pd.read_csv(args.input)
    print(f"Rendering {len(args.targets)} landscape(s) from {len(data)} variants...")
    render_all(data, args.targets, out_dir=args.outdir, save_format=args.format)
    print("Done.")
