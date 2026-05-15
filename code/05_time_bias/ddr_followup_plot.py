"""
Compare time-to-metastasis (metastatic arm) with follow-up duration
(non-metastatic arm) — a check for time-dependent selection bias.

Reported manuscript result:
    Mann-Whitney p < 0.001 — non-metastatic patients have substantially
    longer follow-up than the median time-to-metastasis in the metastatic
    arm, arguing against time-dependent selection bias.

Output:
    figures/followup_comparison.pdf
    figures/followup_comparison.png
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
    met = df[df["metastatic"] == 1]["time_to_metastasis"].dropna().values
    nonmet = df[df["metastatic"] == 0]["followup_duration_months"].dropna().values

    _, pval = stats.mannwhitneyu(met, nonmet, alternative="two-sided")
    p_display = "< 0.001" if pval < 0.001 else f"= {pval:.3f}"
    print(f"Mann-Whitney U: p {p_display}  (n_met={len(met)}, n_nm={len(nonmet)})")

    fig, ax = plt.subplots(figsize=(5, 5))
    rng = np.random.default_rng(42)
    for i, (values, color) in enumerate([(met, "#d62728"), (nonmet, "#1f77b4")]):
        x = i + rng.uniform(-0.12, 0.12, size=len(values))
        ax.scatter(x, values, color=color, edgecolors="white",
                   linewidths=0.5, s=55, alpha=0.85, zorder=3)
        med = float(np.median(values))
        ax.plot([i - 0.25, i + 0.25], [med, med], color=color, linewidth=2.5)
        ax.text(i + 0.28, med, f"median\n{med:.0f}",
                va="center", fontsize=8, color=color)

    y_top = max(met.max(), nonmet.max()) * 1.08
    ax.plot([0, 0, 1, 1], [y_top * 0.97, y_top, y_top, y_top * 0.97], color="black", linewidth=1)
    ax.text(0.5, y_top * 1.01, f"p {p_display}", ha="center", va="bottom", fontsize=9)
    ax.set_xticks([0, 1])
    ax.set_xticklabels(
        ["Metastatic\n(time to metastasis)", "Non-metastatic\n(follow-up)"], fontsize=10
    )
    ax.set_ylabel("Time (months)", fontsize=11)
    ax.set_xlim(-0.6, 1.6)
    ax.spines[["top", "right"]].set_visible(False)
    ax.grid(axis="y", linestyle=":", linewidth=0.5, alpha=0.5)
    plt.tight_layout()

    pdf, png = OUT_DIR / "followup_comparison.pdf", OUT_DIR / "followup_comparison.png"
    plt.savefig(pdf, dpi=300)
    plt.savefig(png, dpi=300)
    print(f"Saved: {pdf}\nSaved: {png}")


if __name__ == "__main__":
    main()
