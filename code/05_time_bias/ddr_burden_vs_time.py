"""
Within the metastatic arm, test whether DDR gnsRV burden correlates with
time-to-metastasis. A negative rho would suggest
that high-burden patients metastasize earlier, which would be evidence of
timing/selection bias. The null result rules this out.

Reported manuscript result:
    Spearman rho = -0.152, p = 0.458, n = 26.

Output:
    figures/ddr_burden_vs_time.pdf
    figures/ddr_burden_vs_time.png
"""

from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy import stats

REPO = Path(__file__).resolve().parents[2]
DATA_FILE = REPO / "data" / "processed" / "epc_cohort.tsv"
OUT_DIR = REPO / "figures"
OUT_DIR.mkdir(exist_ok=True)


def main() -> None:
    df = pd.read_csv(DATA_FILE, sep="\t")
    met = df[df["metastatic"] == 1].dropna(subset=["ddr_burden", "time_to_metastasis"])
    x, y = met["time_to_metastasis"].values, met["ddr_burden"].values

    rho, pval = stats.spearmanr(y, x)
    p_display = "< 0.001" if pval < 0.001 else f"= {pval:.3f}"
    print(f"Spearman rho = {rho:.3f}, p {p_display}, n = {len(met)}")

    fig, ax = plt.subplots(figsize=(5, 4))
    rng = np.random.default_rng(42)
    y_jit = y + rng.uniform(-0.08, 0.08, size=len(y))
    ax.scatter(x, y_jit, color="steelblue", edgecolors="white",
               linewidths=0.6, s=60, alpha=0.85, zorder=3)
    m, b = np.polyfit(x, y, 1)
    line_x = np.linspace(x.min(), x.max(), 100)
    ax.plot(line_x, m * line_x + b, color="gray", linewidth=1, linestyle="--", alpha=0.6)
    ax.text(0.97, 0.97,
            f"Spearman $\\rho$ = {rho:.3f}\n$p$ {p_display}\n$n$ = {len(met)}",
            transform=ax.transAxes, ha="right", va="top", fontsize=9,
            bbox=dict(boxstyle="round,pad=0.3", facecolor="white",
                      edgecolor="lightgray", alpha=0.9))
    ax.set_xlabel("Time to metastasis (months)", fontsize=11)
    ax.set_ylabel("DDR gnsRV burden (count)", fontsize=11)
    ax.set_yticks(sorted(met["ddr_burden"].unique()))
    ax.spines[["top", "right"]].set_visible(False)
    ax.grid(axis="y", linestyle=":", linewidth=0.5, alpha=0.5)
    plt.tight_layout()

    pdf, png = OUT_DIR / "ddr_burden_vs_time.pdf", OUT_DIR / "ddr_burden_vs_time.png"
    plt.savefig(pdf, dpi=300)
    plt.savefig(png, dpi=300)
    print(f"Saved: {pdf}\nSaved: {png}")


if __name__ == "__main__":
    main()
