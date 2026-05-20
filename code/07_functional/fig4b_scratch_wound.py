"""
fig4b_scratch_wound.py — Scratch wound 48 h confluency, rigorous analysis.
LNCaP parental vs KDM6B K973Q knock-in clones C3 (homozygous) and C6 (heterozygous).

Source: data/raw/functional_assays/fig4b_scratch_wound.tsv
Output: figures/fig4b_scratch_wound.{pdf,png}

Layout:
  Row 2  = Ctrl (Parental)
  Row 3  = C3 (homozygous KDM6B K973Q)
  Row 4  = C6 (heterozygous KDM6B K973Q)
  BR1: cols 2–4  (3 technical replicates)
  BR2: cols 5–8  (4 technical replicates)
  BR3: cols 9–11 (3 technical replicates)

Intermediate step files saved to data/processed/functional/:
  scratch_step1_br_means.tsv        average TRs within each BR → 3 BR means/clone
  scratch_step2_log_br_means.tsv    log-transform BR means
  scratch_step3_log_differences.tsv paired log fold changes per BR

Statistics (same logic as fig4c_boyden_chamber.py):
  Ctrl vs C3, Ctrl vs C6 — one-sided paired t-test, Bonferroni × 2
    (directional: K973Q increases wound closure)
  C3  vs C6              — two-sided paired t-test, uncorrected
    (paired within BR; do NOT use unpaired Welch on raw log means)

  Standard output shows raw p AND Bonferroni p for all three comparisons.

Toggles:
  SHOW_C3_C6 = True / False   — C3 vs C6 bracket on plot
  P_DISPLAY  = "bonf" / "raw" — p shown inside plot brackets
"""

from pathlib import Path

# ── USER TOGGLES ──────────────────────────────────────────────────────────────
SHOW_C3_C6 = True  # show C3 vs C6 bracket
P_DISPLAY = "bonf"  # "bonf" | "raw"

# ─────────────────────────────────────────────────────────────────────────────
import matplotlib

matplotlib.use("Agg")
import csv

import matplotlib.pyplot as plt
import numpy as np
from scipy import stats

REPO = Path(__file__).resolve().parents[2]
DATA_FILE = REPO / "data" / "raw" / "functional_assays" / "fig4b_scratch_wound.tsv"
PROCESSED_DIR = REPO / "data" / "processed" / "functional"
PROCESSED_DIR.mkdir(parents=True, exist_ok=True)
FIGURES_DIR = REPO / "figures"
FIGURES_DIR.mkdir(exist_ok=True)

CLONES = ["Ctrl", "C3", "C6"]

# ══════════════════════════════════════════════════════════════════════════════
# STEP 1 — load raw values, average TRs within each BR
# ══════════════════════════════════════════════════════════════════════════════

with open(DATA_FILE) as f:
    rows = list(csv.reader(f, delimiter="\t"))

# header row tells us which columns belong to which BR
header = rows[0]  # ['', 'B1T1','B1T2','B1T3','B2T1',...]


def br_col_groups(header):
    """Parse header → dict {BR_index: [col_indices]} (1-based BR label)."""
    groups = {}
    for ci, h in enumerate(header):
        h = h.strip()
        if not h or not h.startswith("B"):
            continue
        br = int(h[1])  # B1T1 → 1
        groups.setdefault(br, []).append(ci)
    return groups


col_groups = br_col_groups(header)  # {1: [1,2,3], 2: [4,5,6,7], 3: [8,9,10]}
br_ids = sorted(col_groups.keys())  # [1, 2, 3]


def parse_clone_row(row, col_groups, br_ids):
    """Return dict {br_id: mean_of_TRs} for one clone row."""
    out = {}
    for br in br_ids:
        vals = [float(row[ci]) for ci in col_groups[br]]
        out[br] = np.mean(vals)
    return out


# data[clone] = {br_id: mean}
data_raw = {}
for i, clone in enumerate(CLONES):
    data_raw[clone] = parse_clone_row(rows[i + 1], col_groups, br_ids)

# br_data[clone] = list of 3 BR means (ordered by br_id)
br_data = {clone: [data_raw[clone][br] for br in br_ids] for clone in CLONES}

with open(PROCESSED_DIR / "scratch_step1_br_means.tsv", "w") as f:
    f.write("Clone\t" + "\t".join(f"BR{b}_mean" for b in br_ids) + "\n")
    for clone in CLONES:
        v = br_data[clone]
        f.write(clone + "\t" + "\t".join(f"{x:.4f}" for x in v) + "\n")
print(f"Saved: {PROCESSED_DIR / 'scratch_step1_br_means.tsv'}")

# ══════════════════════════════════════════════════════════════════════════════
# STEP 2 — log-transform
# ══════════════════════════════════════════════════════════════════════════════

log_data = {clone: np.log(br_data[clone]) for clone in CLONES}

with open(PROCESSED_DIR / "scratch_step2_log_br_means.tsv", "w") as f:
    f.write("Clone\t" + "\t".join(f"log_BR{b}" for b in br_ids) + "\n")
    for clone in CLONES:
        v = log_data[clone]
        f.write(clone + "\t" + "\t".join(f"{x:.4f}" for x in v) + "\n")
print(f"Saved: {PROCESSED_DIR / 'scratch_step2_log_br_means.tsv'}")

# ══════════════════════════════════════════════════════════════════════════════
# STEP 3 — paired log differences (log fold change per BR)
# ══════════════════════════════════════════════════════════════════════════════

par_log = log_data["Ctrl"]
log_diff = {
    clone: [log_data[clone][i] - par_log[i] for i in range(len(br_ids))]
    for clone in ["C3", "C6"]
}

with open(PROCESSED_DIR / "scratch_step3_log_differences.tsv", "w") as f:
    f.write(
        "Clone\t"
        + "\t".join(f"BR{b}_logFC" for b in br_ids)
        + "\tgeom_mean_FC\tlog_mean\tlog_SE\n"
    )
    for clone in ["C3", "C6"]:
        d = log_diff[clone]
        f.write(
            clone
            + "\t"
            + "\t".join(f"{x:.4f}" for x in d)
            + f"\t{np.exp(np.mean(d)):.4f}\t{np.mean(d):.4f}\t{stats.sem(d):.4f}\n"
        )
print(f"Saved: {PROCESSED_DIR / 'scratch_step3_log_differences.tsv'}")

# ══════════════════════════════════════════════════════════════════════════════
# STATISTICS
# ══════════════════════════════════════════════════════════════════════════════

n = len(br_ids)  # 3
t_crit = stats.t.ppf(0.975, df=n - 1)
results = {}

par_log_arr = np.array(log_data["Ctrl"])

# ── Ctrl vs C3, Ctrl vs C6: paired t-test, store both one- and two-sided ─────
for clone in ["C3", "C6"]:
    cl_log = np.array(log_data[clone])
    diffs = log_diff[clone]
    mean_d = np.mean(diffs)
    se_d = stats.sem(diffs)
    t_stat, p_two = stats.ttest_rel(cl_log, par_log_arr)
    p_one = float(stats.t.sf(t_stat, df=n - 1))  # one-sided upper-tail
    p_two = float(p_two)  # two-sided (= p_one × 2)
    results[("Ctrl", clone)] = dict(
        t=t_stat,
        df=n - 1,
        p_one=p_one,
        p_bonf_one=min(p_one * 2, 1.0),  # Bonferroni × 2
        p_two=p_two,
        p_bonf_two=min(p_two * 2, 1.0),
        p_raw=p_one,
        p_bonf=min(p_one * 2, 1.0),  # used by brackets
        geom_fc=np.exp(mean_d),
        ci_lo=np.exp(mean_d - t_crit * se_d),
        ci_hi=np.exp(mean_d + t_crit * se_d),
        test="paired t-test",
    )

# ── C3 vs C6: two-sided paired t-test on log BR means, uncorrected ───────────
c3_log = np.array(log_data["C3"])
c6_log = np.array(log_data["C6"])
t_stat, p_two = stats.ttest_rel(c3_log, c6_log)
diffs = c3_log - c6_log
mean_d = float(np.mean(diffs))
se_d = float(stats.sem(diffs))
p_two = float(p_two)
results[("C3", "C6")] = dict(
    t=t_stat,
    df=n - 1,
    p_one=p_two / 2,
    p_bonf_one=None,  # not meaningful for secondary
    p_two=p_two,
    p_bonf_two=p_two,  # uncorrected
    p_raw=p_two,
    p_bonf=p_two,
    geom_fc=np.exp(mean_d),
    ci_lo=np.exp(mean_d - t_crit * se_d),
    ci_hi=np.exp(mean_d + t_crit * se_d),
    test="two-sided paired t-test",
)

# ══════════════════════════════════════════════════════════════════════════════
# STANDARD OUTPUT
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
    "Scratch Wound 48 h — Statistical Summary",
    "=" * 62,
    f"Design  : 3 biological replicates (BR1: n=3 TR, BR2: n=4 TR, BR3: n=3 TR)",
    "Scale   : log-transformed BR mean confluency (%)",
    "Primary : one-sided paired t-test, Bonferroni × 2",
    "          (directional: K973Q increases wound closure)",
    "Secondary: two-sided paired t-test, uncorrected (C3 vs C6)",
    f"Settings: SHOW_C3_C6={SHOW_C3_C6}  P_DISPLAY={P_DISPLAY!r}",
    "",
    "BR means (% wound closure at 48 h):",
]
for clone in CLONES:
    v = br_data[clone]
    lines.append(
        f"  {clone:4s}: BR1={v[0]:.2f}  BR2={v[1]:.2f}  BR3={v[2]:.2f}"
        f"  mean={np.mean(v):.2f}  SEM={stats.sem(v):.2f}"
    )

lines += [
    "",
    "Log fold changes per BR (vs Ctrl):",
]
for clone in ["C3", "C6"]:
    d = log_diff[clone]
    lines.append(
        f"  {clone}: BR1={d[0]:.3f}  BR2={d[1]:.3f}  BR3={d[2]:.3f}"
        f"  → geom mean FC={np.exp(np.mean(d)):.3f}"
    )

COL = f"{'Comparison':<18} {'Geom FC':>7} {'95% CI':>14}  {'t':>6} {'df':>2}"
SEP = "-" * 78
PCOL = f"  {'':18}   {'one-sided':>11}  {'Bonf(×2)':>10}    {'two-sided':>11}  {'Bonf(×2)':>10}  sig"

lines += ["", COL, PCOL, SEP]

for pair in [("Ctrl", "C3"), ("Ctrl", "C6"), ("C3", "C6")]:
    r = results[pair]
    is_primary = pair != ("C3", "C6")
    label = {
        ("Ctrl", "C3"): "Ctrl vs C3",
        ("Ctrl", "C6"): "Ctrl vs C6",
        ("C3", "C6"): "C3 vs C6 (2°)",
    }[pair]

    lines.append(
        f"  {label:<18} {r['geom_fc']:>7.3f} "
        f"[{r['ci_lo']:.3f}–{r['ci_hi']:.3f}]  "
        f"{r['t']:>6.3f} {r['df']:>2}"
    )
    if is_primary:
        sig_one = fmt_star(r["p_bonf_one"])
        sig_two = fmt_star(r["p_bonf_two"])
        lines.append(
            f"  {'':18}   "
            f"raw {fmt_p_plain(r['p_one']):>7}  "
            f"Bonf {fmt_p_plain(r['p_bonf_one']):>7}    "
            f"raw {fmt_p_plain(r['p_two']):>7}  "
            f"Bonf {fmt_p_plain(r['p_bonf_two']):>7}  "
            f"1-sided:{sig_one}  2-sided:{sig_two}"
        )
    else:
        lines.append(
            f"  {'':18}   "
            f"{'(not applicable)':>12}  {'':>10}    "
            f"raw {fmt_p_plain(r['p_two']):>7}  {'(uncorrected)':>14}  "
            f"{fmt_star(r['p_two'])}"
        )
lines.append(SEP)

lines += [
    "-" * 74,
    "* primary; Bonferroni × 2",
    "",
    "RECOMMENDED REPORT TEXT",
    "=" * 62,
]
r3 = results[("Ctrl", "C3")]
r6 = results[("Ctrl", "C6")]
lines.append(
    f"KDM6B K973Q homozygous (C3) and heterozygous (C6) LNCaP clones showed "
    f"{r3['geom_fc']:.2f}-fold (95% CI {r3['ci_lo']:.2f}–{r3['ci_hi']:.2f}) "
    f"and {r6['geom_fc']:.2f}-fold (95% CI {r6['ci_lo']:.2f}–{r6['ci_hi']:.2f}) "
    f"greater wound closure at 48 h vs parental, respectively "
    f"(one-sided paired t-test on log-transformed confluency, n=3 biological replicates; "
    f"C3: raw p={fmt_p_plain(r3['p_raw'])}, Bonf p={fmt_p_plain(r3['p_bonf'])} {fmt_star(r3['p_bonf'])}; "
    f"C6: raw p={fmt_p_plain(r6['p_raw'])}, Bonf p={fmt_p_plain(r6['p_bonf'])} {fmt_star(r6['p_bonf'])})."
)

summary = "\n".join(lines)
print("\n" + summary)
with open(PROCESSED_DIR / "scratch_stats_summary.txt", "w") as f:
    f.write(summary + "\n")
print(f"\nSaved: {PROCESSED_DIR / 'scratch_stats_summary.txt'}")

# ══════════════════════════════════════════════════════════════════════════════
# PLOT — Style 2 editorial (terracotta · teal · indigo)
# ══════════════════════════════════════════════════════════════════════════════

# Arithmetic fold changes per BR for plotting
par_br = br_data["Ctrl"]
fc_data = {
    "Ctrl": [1.0, 1.0, 1.0],
    "C3": [br_data["C3"][i] / par_br[i] for i in range(n)],
    "C6": [br_data["C6"][i] / par_br[i] for i in range(n)],
}

TERRA = "#C0513A"
TEAL = "#2E7D6A"
INDIGO = "#4B5EA8"
BARS = [TERRA, TEAL, INDIGO]
SCATTER = ["#E8937E", "#5DB09A", "#8A98D0"]
SUBLABEL = ["Parental\n(LNCaP)", "Homozygous\n(C3)", "Heterozygous\n(C6)"]

labels = CLONES
x = np.arange(len(labels))
means = [np.mean(fc_data[l]) for l in labels]
sems = [stats.sem(fc_data[l]) for l in labels]

rng = np.random.default_rng(42)
fig, ax = plt.subplots(figsize=(5.0, 5.4))

ax.bar(x, means, width=0.60, color=BARS, edgecolor="none", alpha=0.88, zorder=3)
for i, col in enumerate(BARS):
    ax.errorbar(
        x[i],
        means[i],
        yerr=sems[i],
        fmt="none",
        color="#222222",  # color=col,
        capsize=5,
        capthick=1.4,
        linewidth=1.4,
        zorder=6,
    )
for i, (lab, sc, bc) in enumerate(zip(labels, SCATTER, BARS)):
    if lab == "Ctrl":
        continue
    jitter = rng.uniform(-0.10, 0.10, n)
    ax.scatter(
        [x[i] + j for j in jitter],
        fc_data[lab],
        color=sc,
        edgecolors=bc,
        linewidths=0.8,
        s=44,
        zorder=5,
        alpha=0.95,
    )

ax.axhline(1.0, color="#AAAAAA", linewidth=0.8, linestyle=":", zorder=2)

# ── Significance brackets ─────────────────────────────────────────────────────
y_max = max(m + s for m, s in zip(means, sems))
h = y_max * 0.030
idx = {l: i for i, l in enumerate(labels)}

primary = [("Ctrl", "C3", y_max * 1.12), ("Ctrl", "C6", y_max * 1.26)]
secondary = [("C3", "C6", y_max * 1.40)] if SHOW_C3_C6 else []

for a, b, yb in primary + secondary:
    r = results[(a, b)]
    p = r["p_bonf"] if P_DISPLAY == "bonf" else r["p_raw"]
    tag = "Bonf " if P_DISPLAY == "bonf" else "raw "
    xa, xb = x[idx[a]], x[idx[b]]
    ax.plot([xa, xa, xb, xb], [yb, yb + h, yb + h, yb], color="#888888", linewidth=0.9)
    ax.text(
        (xa + xb) / 2,
        yb + h * 1.5,
        f"{fmt_star(p)}  ({tag}{fmt_p(p)})",
        ha="center",
        va="bottom",
        fontsize=12,  # fontsize=9.5,
        fontweight="bold",
        color="#666666",
    )

# ── Axes ──────────────────────────────────────────────────────────────────────
ax.set_xticks(x)
ax.set_xticklabels(SUBLABEL, fontsize=11)
ax.set_ylabel("Wound closure (fold change vs. Parental)", fontsize=12, labelpad=8)
ax.set_ylim(0, y_max * 1.72)
ax.yaxis.grid(True, linestyle="-", linewidth=0.35, color="#E5E5E5", zorder=0)
ax.set_axisbelow(True)
ax.spines["top"].set_visible(False)
ax.spines["right"].set_visible(False)
ax.spines["left"].set_linewidth(0.6)
ax.spines["bottom"].set_linewidth(0.6)
ax.tick_params(axis="x", length=0, labelsize=11)
ax.tick_params(axis="y", length=3, labelsize=10.5)

p_name = "Bonferroni p" if P_DISPLAY == "bonf" else "raw p"
ax.text(
    0.98,
    0.97,
    f"Paired t-test, log scale\n{p_name} shown, Bonf × 2 (primary)",
    transform=ax.transAxes,
    ha="right",
    va="top",
    fontsize=12,  # 8.5,
    fontweight="bold",
    color="#AAAAAA",
    style="italic",
)
ax.set_title(
    "KDM6B K973Q Promotes Wound Closure at 48 h",  # Wound Healing (48H) Assessed by Scratch Assay",
    fontsize=12,
    fontweight="bold",
    pad=10,
    loc="left",
)

plt.tight_layout(pad=1.1)
plt.savefig(FIGURES_DIR / "fig4b_scratch_wound.pdf", dpi=300, bbox_inches="tight")
plt.savefig(FIGURES_DIR / "fig4b_scratch_wound.png", dpi=300, bbox_inches="tight")
print(f"\nSaved: {FIGURES_DIR / 'fig4b_scratch_wound.pdf'}")
print(f"Saved: {FIGURES_DIR / 'fig4b_scratch_wound.png'}")
