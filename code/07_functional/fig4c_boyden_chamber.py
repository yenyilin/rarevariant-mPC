"""
fig4c_boyden_chamber.py — Boyden chamber invasion and migration, rigorous analysis.
LNCaP parental vs KDM6B K973Q knock-in clones C3 (homozygous) and C6 (heterozygous).

Source: data/raw/functional_assays/fig4c_boyden_chamber.tsv
Output: figures/fig4c_boyden_chamber.{pdf,png}

Design: 3 biological replicates × 2 technical replicates per clone.

Intermediate step files saved to data/processed/functional/:
  Step 1 → boyden_step1_br_means.tsv       average technical replicates per BR
  Step 2 → boyden_step2_log_br_means.tsv   log-transform BR means
  Step 3 → boyden_step3_log_differences.tsv paired log differences within each BR
  boyden_stats_summary.txt                 full stats summary

Statistics:
  Parental vs C3, Parental vs C6 — one-sided paired t-test on log BR means.
    Directional hypothesis: K973Q increases invasion/migration.
    Bonferroni × 2 applied to the two primary comparisons.
  C3 vs C6 — two-sided paired t-test on log BR means (secondary, uncorrected).
             Pairing within each BR removes between-replicate noise; do NOT use
             unpaired Welch on raw log means (swamped by between-BR variation).

  Standard output shows raw p AND Bonferroni-corrected p for all three comparisons.
  No omnibus test: 2 pre-specified primary comparisons with Bonferroni already
  controls family-wise error; an omnibus p adds nothing and risks confusion.

  Effect size: geometric mean fold change = exp(mean log difference),
               95% CI back-transformed from log scale.
"""

from pathlib import Path

# ── USER TOGGLES ──────────────────────────────────────────────────────────────
SHOW_C3_C6 = False  # True  → show C3 vs C6 bracket on plots
P_DISPLAY = "bonf"  # "bonf" | "raw"  — p-value shown inside plot brackets

# ─────────────────────────────────────────────────────────────────────────────
import matplotlib

matplotlib.use("Agg")
import csv

import matplotlib.pyplot as plt
import numpy as np
from scipy import stats

REPO = Path(__file__).resolve().parents[2]
DATA_FILE = REPO / "data" / "raw" / "functional_assays" / "fig4c_boyden_chamber.tsv"
PROCESSED_DIR = REPO / "data" / "processed" / "functional"
PROCESSED_DIR.mkdir(parents=True, exist_ok=True)
FIGURES_DIR = REPO / "figures"
FIGURES_DIR.mkdir(exist_ok=True)

ASSAYS = ["Invasion", "Migration"]
CLONES = ["Parental", "C3", "C6"]

# ══════════════════════════════════════════════════════════════════════════════
# STEP 1 — load raw counts, average technical replicates → 3 BR means per clone
# ══════════════════════════════════════════════════════════════════════════════

with open(DATA_FILE) as f:
    rows = list(csv.reader(f, delimiter="\t"))


def parse_block(rows, data_start, n=3):
    out = {}
    for i in range(n):
        row = rows[data_start + i]
        name = row[0].strip()
        vals = [float(x) for x in row[1:] if x.strip()]
        out[name] = vals
    return out


inv_raw = parse_block(rows, data_start=2)
mig_raw = parse_block(rows, data_start=7)


def br_means(raw6):
    return [(raw6[0] + raw6[1]) / 2, (raw6[2] + raw6[3]) / 2, (raw6[4] + raw6[5]) / 2]


br_data = {}
for assay, raw in [("Invasion", inv_raw), ("Migration", mig_raw)]:
    br_data[assay] = {c: br_means(raw[c]) for c in CLONES}

with open(PROCESSED_DIR / "boyden_step1_br_means.tsv", "w") as f:
    f.write("Assay\tClone\tBR1_mean\tBR2_mean\tBR3_mean\n")
    for assay in ASSAYS:
        for clone in CLONES:
            v = br_data[assay][clone]
            f.write(f"{assay}\t{clone}\t{v[0]:.2f}\t{v[1]:.2f}\t{v[2]:.2f}\n")
print(f"Saved: {PROCESSED_DIR / 'boyden_step1_br_means.tsv'}")

# ══════════════════════════════════════════════════════════════════════════════
# STEP 2 — log-transform BR means
# ══════════════════════════════════════════════════════════════════════════════

log_data = {}
for assay in ASSAYS:
    log_data[assay] = {c: np.log(br_data[assay][c]) for c in CLONES}

with open(PROCESSED_DIR / "boyden_step2_log_br_means.tsv", "w") as f:
    f.write("Assay\tClone\tlog_BR1\tlog_BR2\tlog_BR3\n")
    for assay in ASSAYS:
        for clone in CLONES:
            v = log_data[assay][clone]
            f.write(f"{assay}\t{clone}\t{v[0]:.4f}\t{v[1]:.4f}\t{v[2]:.4f}\n")
print(f"Saved: {PROCESSED_DIR / 'boyden_step2_log_br_means.tsv'}")

# ══════════════════════════════════════════════════════════════════════════════
# STEP 3 — paired log differences (log fold change per BR)
# ══════════════════════════════════════════════════════════════════════════════

log_diff = {}
for assay in ASSAYS:
    par = log_data[assay]["Parental"]
    log_diff[assay] = {
        clone: [log_data[assay][clone][i] - par[i] for i in range(3)]
        for clone in ["C3", "C6"]
    }

with open(PROCESSED_DIR / "boyden_step3_log_differences.tsv", "w") as f:
    f.write(
        "Assay\tClone\tBR1_logFC\tBR2_logFC\tBR3_logFC\tgeom_mean_FC\tlog_mean\tlog_SE\n"
    )
    for assay in ASSAYS:
        for clone in ["C3", "C6"]:
            d = log_diff[assay][clone]
            f.write(
                f"{assay}\t{clone}\t{d[0]:.4f}\t{d[1]:.4f}\t{d[2]:.4f}\t"
                f"{np.exp(np.mean(d)):.4f}\t{np.mean(d):.4f}\t{stats.sem(d):.4f}\n"
            )
print(f"Saved: {PROCESSED_DIR / 'boyden_step3_log_differences.tsv'}")

# ══════════════════════════════════════════════════════════════════════════════
# STATISTICS
# ══════════════════════════════════════════════════════════════════════════════

t_crit = stats.t.ppf(0.975, df=2)  # for 95% CI, df = n-1 = 2

results = {}

for assay in ASSAYS:
    results[assay] = {}
    par_log = np.array(log_data[assay]["Parental"])

    # ── Parental vs C3 and Parental vs C6: one-sided paired t-test ────────────
    for clone in ["C3", "C6"]:
        cl_log = np.array(log_data[assay][clone])
        diffs = log_diff[assay][clone]
        mean_d = np.mean(diffs)
        se_d = stats.sem(diffs)
        t_stat, _ = stats.ttest_rel(cl_log, par_log)
        p_raw = float(stats.t.sf(t_stat, df=2))  # one-sided upper-tail
        p_bonf = min(p_raw * 2, 1.0)  # Bonferroni × 2 primary
        results[assay][("Parental", clone)] = dict(
            t=t_stat,
            df=2,
            p_raw=p_raw,
            p_bonf=p_bonf,
            geom_fc=np.exp(mean_d),
            ci_lo=np.exp(mean_d - t_crit * se_d),
            ci_hi=np.exp(mean_d + t_crit * se_d),
            test="one-sided paired t-test",
        )

    # ── C3 vs C6: two-sided paired t-test on log BR means (secondary) ───────────
    # Must NOT compare raw log BR means — between-BR noise (BR3 high for all
    # clones) swamps the C3 vs C6 difference (gives t=-0.72, p=0.51).
    # Pairing within each BR removes that shared variation. Equivalent to Welch
    # on log fold changes (R approach), but more powerful as it uses the design.
    c3_log = np.array(log_data[assay]["C3"])
    c6_log = np.array(log_data[assay]["C6"])
    t_stat, p_two = stats.ttest_rel(c6_log, c3_log)  # paired, two-sided
    diffs = c6_log - c3_log
    mean_d = float(np.mean(diffs))
    se_d = float(stats.sem(diffs))
    results[assay][("C3", "C6")] = dict(
        t=t_stat,
        df=2,
        p_raw=p_two,
        p_bonf=p_two,  # secondary: uncorrected
        geom_fc=np.exp(mean_d),
        ci_lo=np.exp(mean_d - t_crit * se_d),
        ci_hi=np.exp(mean_d + t_crit * se_d),
        test="two-sided paired t-test",
    )

# ══════════════════════════════════════════════════════════════════════════════
# STANDARD OUTPUT — both raw and Bonferroni for all comparisons
# ══════════════════════════════════════════════════════════════════════════════


def fmt_star(p):
    return "***" if p < 0.001 else "**" if p < 0.01 else "*" if p < 0.05 else "ns"


def fmt_p_plain(p):
    if p >= 0.001:
        return f"{p:.3f}"
    e = int(np.floor(np.log10(p)))
    return f"{p / 10**e:.2f}e{e}"


def fmt_p(p):
    if p >= 0.001:
        return f"p = {p:.3f}"
    e = int(np.floor(np.log10(p)))
    return f"p = {p / 10**e:.2f}×10$^{{{e}}}$"


lines = [
    "Boyden Chamber Assay — Statistical Summary",
    "=" * 62,
    "Design  : 3 biological replicates × 2 technical replicates",
    "Scale   : log-transformed BR mean counts",
    "Primary : one-sided paired t-test, Bonferroni × 2",
    "          (directional: K973Q increases invasion/migration)",
    "Secondary: two-sided paired t-test, uncorrected (C3 vs C6)",
    f"Settings: SHOW_C3_C6={SHOW_C3_C6}  P_DISPLAY={P_DISPLAY!r}",
    "",
    f"{'Comparison':<22} {'Geom FC':>7} {'95% CI':>14}  "
    f"{'t':>6} {'df':>4}  {'raw p':>7}  {'Bonf p':>7}  sig",
    "-" * 82,
]

PAIR_ORDER = [("Parental", "C3"), ("Parental", "C6"), ("C3", "C6")]
LABELS_ASSAY = {
    ("Parental", "C3"): "Parental vs C3 *",
    ("Parental", "C6"): "Parental vs C6 *",
    ("C3", "C6"): "C3 vs C6",
}

for assay in ASSAYS:
    lines.append(f"\n  {assay.upper()}")
    for pair in PAIR_ORDER:
        r = results[assay][pair]
        df = f"{r['df']:.0f}" if r["df"] is not None else "—"
        lab = LABELS_ASSAY[pair]
        bonf_note = "" if pair != ("C3", "C6") else "(uncorr)"
        lines.append(
            f"  {lab:<22} {r['geom_fc']:>7.3f} "
            f"[{r['ci_lo']:.3f}–{r['ci_hi']:.3f}]  "
            f"{r['t']:>6.3f} {df:>4}  "
            f"{fmt_p_plain(r['p_raw']):>7}  "
            f"{fmt_p_plain(r['p_bonf']):>7} {bonf_note}  "
            f"{fmt_star(r['p_bonf'])}"
        )

lines += [
    "-" * 82,
    "* primary comparison; Bonferroni × 2",
    "",
    "RECOMMENDED REPORT TEXT",
    "=" * 62,
    "",
]

for assay in ASSAYS:
    r3 = results[assay][("Parental", "C3")]
    r6 = results[assay][("Parental", "C6")]
    lines.append(f"{assay}:")
    lines.append(
        f"  C3: {r3['geom_fc']:.2f}-fold (95% CI {r3['ci_lo']:.2f}–{r3['ci_hi']:.2f}), "
        f"paired t(2)={r3['t']:.2f}, raw p={fmt_p_plain(r3['p_raw'])}, "
        f"Bonf p={fmt_p_plain(r3['p_bonf'])} {fmt_star(r3['p_bonf'])}"
    )
    lines.append(
        f"  C6: {r6['geom_fc']:.2f}-fold (95% CI {r6['ci_lo']:.2f}–{r6['ci_hi']:.2f}), "
        f"paired t(2)={r6['t']:.2f}, raw p={fmt_p_plain(r6['p_raw'])}, "
        f"Bonf p={fmt_p_plain(r6['p_bonf'])} {fmt_star(r6['p_bonf'])}"
    )
    lines.append("")

summary = "\n".join(lines)
print("\n" + summary)
with open(PROCESSED_DIR / "boyden_stats_summary.txt", "w") as f:
    f.write(summary + "\n")
print(f"Saved: {PROCESSED_DIR / 'boyden_stats_summary.txt'}")

# ══════════════════════════════════════════════════════════════════════════════
# PLOT HELPERS
# ══════════════════════════════════════════════════════════════════════════════

# Arithmetic fold changes per BR (bars = arithmetic mean, err = SEM)
fc_data = {}
for assay in ASSAYS:
    par_br = br_data[assay]["Parental"]
    fc_data[assay] = {
        "Parental": [1.0, 1.0, 1.0],
        "C3": [br_data[assay]["C3"][i] / par_br[i] for i in range(3)],
        "C6": [br_data[assay]["C6"][i] / par_br[i] for i in range(3)],
    }

SUBLABEL = ["Parental\n(LNCaP)", "Homozygous\n(C3)", "Heterozygous\n(C6)"]


def annotation_label():
    p_name = "Bonferroni p" if P_DISPLAY == "bonf" else "raw p"
    return f"Paired t-test, log scale\n{p_name} shown, \nBonf × 2 (primary)"


def add_brackets(ax, assay, labels, x, means, sems, style="default"):
    bracket_col = "#333333" if style == "default" else "#888888"
    text_col = "#222222" if style == "default" else "#666666"
    y_max = max(m + s for m, s in zip(means, sems))
    h = y_max * 0.030
    idx = {l: i for i, l in enumerate(labels)}

    primary = [("Parental", "C3", y_max * 1.12), ("Parental", "C6", y_max * 1.26)]
    secondary = [("C3", "C6", y_max * 1.40)] if SHOW_C3_C6 else []

    for a, b, yb in primary + secondary:
        r = results[assay][(a, b)]
        p = r["p_bonf"] if P_DISPLAY == "bonf" else r["p_raw"]
        tag = "Bonf " if P_DISPLAY == "bonf" else "raw "
        xa, xb = x[idx[a]], x[idx[b]]
        ax.plot(
            [xa, xa, xb, xb], [yb, yb + h, yb + h, yb], color=bracket_col, linewidth=0.9
        )
        ax.text(
            (xa + xb) / 2,
            yb + h * 1.5,
            f"{fmt_star(p)}  ({tag}{fmt_p(p)})",
            ha="center",
            va="bottom",
            fontsize=10,  # 9.5,
            fontweight="bold",
            color=text_col,
        )
    return y_max


# ══════════════════════════════════════════════════════════════════════════════
# PLOT A — default style (matching plot_4c.py)
# ══════════════════════════════════════════════════════════════════════════════

AMBER = "#E8A838"
GREEN = "#27AE60"
BLUE = "#2980B9"

fig, axes = plt.subplots(1, 2, figsize=(9.0, 5.2), sharey=False)
for ax, assay in zip(axes, ASSAYS):
    labels = CLONES
    x = np.arange(len(labels))
    means = [np.mean(fc_data[assay][l]) for l in labels]
    sems = [stats.sem(fc_data[assay][l]) for l in labels]

    ax.bar(
        x,
        means,
        width=0.55,
        color=[AMBER, GREEN, BLUE],
        edgecolor="#333333",
        linewidth=0.8,
        alpha=0.85,
        zorder=3,
    )
    ax.errorbar(
        x,
        means,
        yerr=sems,
        fmt="none",
        color="#222222",
        capsize=5,
        capthick=1.2,
        linewidth=1.2,
        zorder=6,
    )
    for i, lab in enumerate(labels):
        if lab == "Parental":
            continue
        jitter = np.linspace(-0.10, 0.10, 3)
        ax.scatter(
            [x[i] + j for j in jitter],
            fc_data[assay][lab],
            color="#333333",
            s=38,
            zorder=5,
            edgecolors="white",
            linewidths=0.5,
        )

    ax.axhline(1.0, color="#888888", linewidth=0.8, linestyle=":", zorder=2)
    y_max = add_brackets(ax, assay, labels, x, means, sems, style="default")
    ax.set_xticks(x)
    ax.set_xticklabels(SUBLABEL, fontsize=12)
    ax.set_ylabel(f"{assay} (fold change vs. Parental)", fontsize=12)
    ax.set_ylim(0, y_max * 1.72)
    ax.yaxis.grid(True, linestyle="--", linewidth=0.5, alpha=0.5, zorder=0)
    ax.set_axisbelow(True)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.tick_params(axis="y", labelsize=11)
    ax.set_title(f"KDM6B K973Q: {assay.lower()}", fontsize=12, fontweight="bold", pad=8)
    ax.text(
        0.98,
        0.97,
        annotation_label(),
        transform=ax.transAxes,
        ha="right",
        va="top",
        fontsize=8.5,
        color="#666666",
        style="italic",
    )

fig.suptitle(
    "Cell Invasion and Migration\nAssessed by Boyden Chamber Assay",
    fontsize=12,
    fontweight="bold",
    y=1.01,
)
plt.tight_layout()
plt.savefig(FIGURES_DIR / "fig4c_boyden_chamber_alt.pdf", dpi=300, bbox_inches="tight")
plt.savefig(FIGURES_DIR / "fig4c_boyden_chamber_alt.png", dpi=300, bbox_inches="tight")
print(f"\nSaved: {FIGURES_DIR / 'fig4c_boyden_chamber_alt.pdf'}")
print(f"Saved: {FIGURES_DIR / 'fig4c_boyden_chamber_alt.png'}")

# ══════════════════════════════════════════════════════════════════════════════
# PLOT B — Style 2 editorial (terracotta · teal · indigo)
# ══════════════════════════════════════════════════════════════════════════════

TERRA = "#C0513A"
TEAL = "#2E7D6A"
INDIGO = "#4B5EA8"
S2_BARS = [TERRA, TEAL, INDIGO]
S2_SCATTER = ["#E8937E", "#5DB09A", "#8A98D0"]
rng = np.random.default_rng(42)

fig2, axes2 = plt.subplots(1, 2, figsize=(9.0 * 0.8, 5.2), sharey=False)
for ax, assay in zip(axes2, ASSAYS):
    labels = CLONES
    x = np.arange(len(labels))
    means = [np.mean(fc_data[assay][l]) for l in labels]
    sems = [stats.sem(fc_data[assay][l]) for l in labels]

    ax.bar(x, means, width=0.60, color=S2_BARS, edgecolor="none", alpha=0.88, zorder=3)
    for i, col in enumerate(S2_BARS):
        ax.errorbar(
            x[i],
            means[i],
            yerr=sems[i],
            fmt="none",
            color="#222222",
            capsize=5,
            capthick=1.4,
            linewidth=1.4,
            zorder=6,
        )
    for i, (lab, sc, bc) in enumerate(zip(labels, S2_SCATTER, S2_BARS)):
        if lab == "Parental":
            continue
        jitter = rng.uniform(-0.10, 0.10, 3)
        ax.scatter(
            [x[i] + j for j in jitter],
            fc_data[assay][lab],
            color=sc,
            edgecolors=bc,
            linewidths=0.8,
            s=42,
            zorder=5,
            alpha=0.95,
        )

    ax.axhline(1.0, color="#AAAAAA", linewidth=0.8, linestyle=":", zorder=2)
    y_max = add_brackets(ax, assay, labels, x, means, sems, style="editorial")
    ax.set_xticks(x)
    ax.set_xticklabels(SUBLABEL, fontsize=10)
    ax.set_ylabel(f"{assay} (fold change vs. Parental)", fontsize=12, labelpad=8)
    ax.set_ylim(0, y_max * 1.72)
    ax.yaxis.grid(True, linestyle="-", linewidth=0.35, color="#E5E5E5", zorder=0)
    ax.set_axisbelow(True)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.spines["left"].set_linewidth(0.6)
    ax.spines["bottom"].set_linewidth(0.6)
    ax.tick_params(axis="x", length=0, labelsize=9)
    ax.tick_params(axis="y", length=3, labelsize=11)
    ax.set_title(
        f"KDM6B K973Q: {assay.lower()}",
        fontsize=12,
        fontweight="bold",
        pad=8,
        loc="left",
    )
    ax.text(
        0.98,
        0.97,
        annotation_label(),
        transform=ax.transAxes,
        ha="right",
        va="top",
        fontsize=10,  # 8.5,
        fontweight="bold",
        color="#AAAAAA",
        style="italic",
    )

fig2.suptitle(
    "Cell Invasion and Migration\nAssessed by Boyden Chamber Assay",
    fontsize=12,
    fontweight="bold",
    y=1.01,
)
plt.tight_layout(pad=1.1)
plt.savefig(FIGURES_DIR / "fig4c_boyden_chamber.pdf", dpi=300, bbox_inches="tight")
plt.savefig(FIGURES_DIR / "fig4c_boyden_chamber.png", dpi=300, bbox_inches="tight")
print(f"Saved: {FIGURES_DIR / 'fig4c_boyden_chamber.pdf'}")
print(f"Saved: {FIGURES_DIR / 'fig4c_boyden_chamber.png'}")
