"""
Macroscopic comparison of the main functional-variant library against the
L108G-fixed control library.

Compares the distributions of affinity and expression, and the fraction of
variants with complete loss of each property, between the main 9,600-variant
library and the 1,920-variant L108G-fixed control. Including the non-functional
L108G anchor is expected to shift the distributions downward and increase the
complete-loss fraction, supporting the strategy of recombining only
individually functional substitutions.

Pipeline position
-----------------
Standalone. Takes two measurement tables: the main library and the L108G
library.

Corresponding manuscript items
------------------------------
Fig. 3r,s.
"""

import argparse
import warnings

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from scipy import stats

warnings.filterwarnings("ignore")

plt.rcParams["svg.fonttype"] = "none"
plt.rcParams["pdf.fonttype"] = 42
plt.rcParams["ps.fonttype"] = 42

OCCUPANCY_FLOOR = 0.05
ID_CANDIDATES = ["Mutation_Description", "ID", "Oligo Name", "Variant",
                 "Sequence", "Name"]


def match_id(df):
    id_col = next((c for c in ID_CANDIDATES if c in df.columns), df.columns[0])
    df = df.copy()
    df["Match_ID"] = df[id_col].astype(str).str.upper().str.replace("WT", "PARENTAL")
    return df


def log_column(df, keyword):
    """Return (raw_col, log_col) for a metric, creating the floored log2
    column if only the raw column is present."""
    raw = next((c for c in df.columns
                if keyword in c.lower() and not c.startswith("Log_")
                and not c.lower().endswith("_std")), None)
    if raw is None:
        return None, None
    log_col = f"Log_{raw}"
    if log_col not in df.columns:
        df[log_col] = np.log2(pd.to_numeric(df[raw], errors="coerce")
                              .clip(lower=OCCUPANCY_FLOOR))
    return raw, log_col


def stars(p):
    return ("****" if p <= 1e-4 else "***" if p <= 1e-3
            else "**" if p <= 1e-2 else "*" if p <= 5e-2 else "ns")


def complete_loss(df, raw_col):
    """Fraction (%) of variants with complete loss, defined as raw value (minus
    its replicate s.d. when available) at or below zero, with a 95% CI."""
    if raw_col not in df.columns:
        return 0.0, 0.0
    std_col = f"{raw_col}_std"
    values = pd.to_numeric(df[raw_col], errors="coerce")
    if std_col in df.columns:
        mask = (values <= 0) | ((values - df[std_col].fillna(0)) <= 0)
    else:
        mask = values <= 0
    p = mask.mean()
    ci = 1.96 * np.sqrt(p * (1 - p) / len(df)) * 100
    return p * 100, ci


def compare(main_df, l108g_df, keyword, short_name, out_dir, save_format="svg"):
    raw_m, log_m = log_column(main_df, keyword)
    raw_l, log_l = log_column(l108g_df, keyword)
    if not log_m or not log_l:
        print(f"   Skipping {short_name}: metric not found in both tables.")
        return

    main_scores = main_df[log_m].dropna()
    l108g_scores = l108g_df[log_l].dropna()
    label_m = f"Functional library (N={len(main_scores):,})"
    label_l = f"L108G control (N={len(l108g_scores):,})"

    parental = main_df[main_df["Match_ID"] == "PARENTAL"]
    parental_val = parental[log_m].dropna().to_numpy()
    parental_val = parental_val[0] if len(parental_val) else 0.0

    _, p_val = stats.mannwhitneyu(main_scores, l108g_scores, alternative="two-sided")
    loss_m, err_m = complete_loss(main_df, raw_m)
    loss_l, err_l = complete_loss(l108g_df, raw_l)

    violin = pd.concat([
        pd.DataFrame({"Library": label_m, "Score": main_scores}),
        pd.DataFrame({"Library": label_l, "Score": l108g_scores}),
    ], ignore_index=True)

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 6))
    colors = ["#B0BEC5", "#D62728"]

    sns.violinplot(data=violin, x="Library", y="Score", palette=colors,
                   ax=ax1, inner="quartile", linewidth=1.5, cut=0)
    ax1.axhline(parental_val, color="black", linestyle="--", linewidth=2,
                label=f"Parental ({parental_val:.2f})")
    y_max, y_min = violin["Score"].max(), violin["Score"].min()
    span = y_max - y_min
    bar_y = y_max + span * 0.05
    ax1.plot([0, 0, 1, 1], [bar_y, bar_y + span * 0.02, bar_y + span * 0.02, bar_y],
             lw=1.5, c="black")
    ax1.text(0.5, bar_y + span * 0.02, stars(p_val), ha="center", va="bottom",
             fontsize=16, fontweight="bold")
    ax1.set_title(f"{short_name} distribution")
    ax1.set_ylabel("Log2 score"); ax1.set_xlabel(""); ax1.legend(loc="upper right")

    bars = ax2.bar([label_m, label_l], [loss_m, loss_l], yerr=[err_m, err_l],
                   capsize=8, color=colors, edgecolor="black", linewidth=1.5)
    for bar, value, err in zip(bars, [loss_m, loss_l], [err_m, err_l]):
        ax2.text(bar.get_x() + bar.get_width() / 2, value + err + 1,
                 f"{value:.1f}%", ha="center", va="bottom", fontweight="bold")
    ax2.set_title(f"Complete loss of {short_name.lower()}")
    ax2.set_ylabel("Percentage (%)")
    ax2.set_ylim(0, max(100, loss_m + err_m + 10, loss_l + err_l + 10))
    sns.despine(top=True, right=True)

    plt.tight_layout()
    fig_path = f"{out_dir}/MacroComparison_{short_name}.{save_format}"
    plt.savefig(fig_path, dpi=300, bbox_inches="tight", format=save_format)
    plt.close()

    pd.DataFrame({label_m: main_scores.reset_index(drop=True),
                  label_l: l108g_scores.reset_index(drop=True)}).to_csv(
        f"{out_dir}/SourceData_MacroViolin_{short_name}.csv", index=False)

    print(f"   {short_name}: Mann-Whitney p = {p_val:.2e}; complete loss "
          f"{loss_m:.1f}% vs {loss_l:.1f}% ({fig_path})")


if __name__ == "__main__":
    import os

    parser = argparse.ArgumentParser(
        description="Macroscopic comparison of the functional and L108G libraries"
    )
    parser.add_argument("--main", required=True,
                        help="Main functional-variant library table (CSV/XLSX)")
    parser.add_argument("--l108g", required=True,
                        help="L108G-fixed control library table (CSV/XLSX)")
    parser.add_argument("--outdir", default="output_l108g")
    parser.add_argument("--format", default="svg", choices=["svg", "png"])
    args = parser.parse_args()

    def load(path):
        return (pd.read_excel(path) if path.lower().endswith(("xlsx", "xls"))
                else pd.read_csv(path))

    os.makedirs(args.outdir, exist_ok=True)
    main_df = match_id(load(args.main))
    l108g_df = match_id(load(args.l108g))
    print(f"Main library: {len(main_df)} rows; L108G control: {len(l108g_df)} rows.")

    for keyword, short in [("occupancy", "Affinity"), ("expression", "Expression")]:
        compare(main_df, l108g_df, keyword, short, args.outdir,
                save_format=args.format)
    print("Done.")
