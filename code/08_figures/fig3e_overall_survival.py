"""
fig3e_overall_survival.py — Fig. 3E Kaplan-Meier overall survival (Altered vs Unaltered)
----------------------------------------------------------------------------------------
Source: data/raw/cbioportal/fig3e_overall_survival.tsv  (cBioPortal export)
Output: figures/fig3e_overall_survival.{pdf,png}

Plots pre-computed KM curves from cBioPortal survival rates.
Log-rank p-value computed from raw event data.
"""

from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
import numpy as np
import pandas as pd

REPO = Path(__file__).resolve().parents[2]

# ── Toggle ────────────────────────────────────────────────────────────────────
SHOW_AT_RISK = True  # set False to hide number-at-risk table
AT_RISK_TIMES = [0, 24, 48, 72, 96, 120, 144, 168]  # months


# ── Parse two-section file ────────────────────────────────────────────────────
def parse_overall(path):
    groups = {}
    current = None
    rows = []
    with open(path) as f:
        for line in f:
            line = line.rstrip("\n")
            if line.startswith("Altered group"):
                current = "Altered"
                rows = []
            elif line.startswith("Unaltered group"):
                if current and rows:
                    groups[current] = pd.DataFrame(
                        rows,
                        columns=[
                            "case_id",
                            "study",
                            "n_risk",
                            "status",
                            "surv",
                            "time",
                        ],
                    )
                current = "Unaltered"
                rows = []
            elif (
                line.startswith("Case ID")
                or line.strip() == ""
                or line.startswith("Overall")
            ):
                continue
            else:
                parts = line.split("\t")
                if len(parts) == 6:
                    rows.append(parts)
    if current and rows:
        groups[current] = pd.DataFrame(
            rows, columns=["case_id", "study", "n_risk", "status", "surv", "time"]
        )
    for g in groups:
        groups[g]["time"] = groups[g]["time"].astype(float)
        groups[g]["surv"] = groups[g]["surv"].astype(float)
        groups[g]["event"] = (groups[g]["status"] == "deceased").astype(int)
        groups[g]["n_risk"] = groups[g]["n_risk"].astype(int)
        groups[g] = groups[g].sort_values("time").reset_index(drop=True)
    return groups


# ── Log-rank test (manual) ────────────────────────────────────────────────────
def logrank_test(df1, df2):
    from scipy.stats import chi2

    all_times = np.union1d(
        df1.loc[df1["event"] == 1, "time"].values,
        df2.loc[df2["event"] == 1, "time"].values,
    )
    O1 = E1 = var = 0.0
    for t in all_times:
        n1 = (df1["time"] >= t).sum()
        n2 = (df2["time"] >= t).sum()
        d1 = ((df1["time"] == t) & (df1["event"] == 1)).sum()
        d2 = ((df2["time"] == t) & (df2["event"] == 1)).sum()
        n = n1 + n2
        d = d1 + d2
        if n < 2:
            continue
        O1 += d1
        E1 += n1 * d / n
        var += n1 * n2 * d * (n - d) / (n**2 * (n - 1)) if n > 1 else 0
    if var == 0:
        return 1.0
    stat = (O1 - E1) ** 2 / var
    return chi2.sf(stat, df=1)


# ── KM step function with censoring marks ─────────────────────────────────────
def plot_km(ax, df, color, label):
    # Plot pre-computed survival rate as step function
    times = [0] + list(df["time"])
    surv = [1.0] + list(df["surv"])
    ax.step(times, surv, where="post", color=color, linewidth=2, label=label)

    # Censoring tick marks
    censored = df[df["event"] == 0]
    ax.plot(
        censored["time"],
        censored["surv"],
        "|",
        color=color,
        markersize=6,
        markeredgewidth=1.5,
    )


# ── Main ──────────────────────────────────────────────────────────────────────
groups = parse_overall(REPO / "data" / "raw" / "cbioportal" / "fig3e_overall_survival.tsv")
df_alt = groups["Altered"]
df_unt = groups["Unaltered"]

p_value = logrank_test(df_alt, df_unt)

n_alt = len(df_alt)
n_unt = len(df_unt)


# ── P-value formatting ────────────────────────────────────────────────────────
def fmt_pvalue(p):
    if p >= 0.001:
        return f"p = {p:.3f}"
    exp = int(np.floor(np.log10(p)))
    coef = p / 10**exp
    return f"p = {coef:.2f} × 10$^{{{exp}}}$"  # e.g. p = 1.88 × 10⁻⁹


# ── Number-at-risk helper ─────────────────────────────────────────────────────
def n_at_risk(df, t):
    return int((df["time"] >= t).sum())


# ── Plot ──────────────────────────────────────────────────────────────────────
COLOR_ALT = "#C0392B"  # red — altered (carrier)
COLOR_UNT = "#2471A3"  # blue — unaltered

fig_height = 5.8 if SHOW_AT_RISK else 4.8
fig, ax = plt.subplots(figsize=(5.5, fig_height))

plot_km(ax, df_alt, COLOR_ALT, f"Altered (n={n_alt})")
plot_km(ax, df_unt, COLOR_UNT, f"Unaltered (n={n_unt})")

# ── P-value annotation ────────────────────────────────────────────────────────
p_text = fmt_pvalue(p_value)
ax.text(
    0.97,
    0.55,
    p_text,
    transform=ax.transAxes,
    ha="right",
    va="top",
    fontsize=10,
    fontweight="bold",
    bbox=dict(
        boxstyle="round,pad=0.3", facecolor="white", edgecolor="#cccccc", linewidth=0.8
    ),
    zorder=10,
)

# ── Axes ──────────────────────────────────────────────────────────────────────
ax.set_xlabel("Time (months)", fontsize=10)
ax.set_ylabel("Overall survival probability", fontsize=10)
ax.set_xlim(left=0)
ax.set_ylim(0, 1.05)
ax.yaxis.set_major_formatter(mticker.PercentFormatter(xmax=1, decimals=0))
ax.yaxis.grid(True, linestyle="--", linewidth=0.5, alpha=0.5, zorder=0)
ax.set_axisbelow(True)
ax.spines["top"].set_visible(False)
ax.spines["right"].set_visible(False)

ax.set_title(
    "Overall survival: gnsRV carriers vs non-carriers\n(cBioPortal, TCGA + SU2C)",
    fontsize=10,
    fontweight="bold",
    pad=8,
)
ax.legend(fontsize=9, frameon=False, loc="lower left")

# ── Number-at-risk table ──────────────────────────────────────────────────────
if SHOW_AT_RISK:
    nars_alt = [n_at_risk(df_alt, t) for t in AT_RISK_TIMES]
    nars_unt = [n_at_risk(df_unt, t) for t in AT_RISK_TIMES]

    x_max = ax.get_xlim()[1]
    row_y = [-0.18, -0.26]  # in axes coordinates below x-axis
    colors = [COLOR_ALT, COLOR_UNT]
    labels = ["Altered", "Unaltered"]
    nar_data = [nars_alt, nars_unt]

    ax.text(
        -0.22,
        -0.14,
        "No. at risk",
        transform=ax.transAxes,
        ha="left",
        va="top",
        fontsize=8,
        fontweight="bold",
        color="#333",
    )

    for row, (ry, color, label, nars) in enumerate(
        zip(row_y, colors, labels, nar_data)
    ):
        ax.text(
            -0.22,
            ry,
            label,
            transform=ax.transAxes,
            ha="left",
            va="top",
            fontsize=8,
            color=color,
            fontweight="bold",
        )
        for t, n in zip(AT_RISK_TIMES, nars):
            x_pos = t / x_max  # normalise to axes coordinates
            ax.text(
                x_pos,
                ry,
                str(n),
                transform=ax.transAxes,
                ha="center",
                va="top",
                fontsize=8,
                color=color,
            )

    plt.subplots_adjust(bottom=0.22)

plt.tight_layout(rect=[0, 0.0 if not SHOW_AT_RISK else 0.05, 1, 1])
plt.savefig(REPO / "figures" / "fig3e_overall_survival.pdf", dpi=300, bbox_inches="tight")
plt.savefig(REPO / "figures" / "fig3e_overall_survival.png", dpi=300, bbox_inches="tight")
print(
    f"Saved: figures/fig3e_overall_survival.pdf / .png  |  log-rank {p_text}  |  p_raw = {p_value:.3e}"
)
