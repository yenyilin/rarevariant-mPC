#!/usr/bin/env python3
"""
generate_table1.py
──────────────────
Reads a TSV file (one patient per row) and writes:
  1. table1_output.xlsx   — Table 1 (main MS) + Supplementary
  2. swimmer_plot.pdf     — patient follow-up timelines
  3. love_plot.pdf        — covariate balance (standardized mean differences)

Input format
  Column 1 : patient identifier (any header)
  Column 2 : arm  1 = metastatic  0 = non-metastatic
  Column 3+: clinical variables (headers matched case-insensitively)

Usage
  python generate_table1.py cohort.tsv
  python generate_table1.py cohort.tsv -o results/

Requirements
  pip install pandas scipy openpyxl matplotlib
"""

import argparse
import os
import sys

import matplotlib
import numpy as np
import pandas as pd
from scipy import stats

matplotlib.use("Agg")
import matplotlib.patches as mpatches
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

# ─────────────────────────────────────────────────────────────────────────────
# 1a.  Header rename map — applied at load time
#
# The clinical TSV (e.g. data/synthetic/clinical_synthetic.tsv) uses sentence-case headers
# and ships pre-derived "months" columns (no dates). All downstream code
# below uses the original ALL-CAPS constants; this dict translates the new
# headers into the constants the rest of the script expects.
#
# To extend: add `"<new header>": "<OLD CONSTANT>"` entries.
# To rename a constant: edit COLUMN_RENAME *and* the matching constant
# in the registry (Section 1).
# ─────────────────────────────────────────────────────────────────────────────

COLUMN_RENAME = {
    "Patient ID": "PATIENT_ID",
    "Metastasis (1=yes, 0=no)": "METASTASIS",
    "Ancestry (inferred)": "ETHNICITY (INFERRED)",
    "Age at diagnosis (years)": "AGE AT DX",
    "PSA at diagnosis (ng/mL)": "PSA AT DX",
    "Biopsy Gleason score 1": "BIOPSY GLEASON SCORE 1",
    "Biopsy Gleason score 2": "BIOPSY GLEASON SCORE 2",
    "Biopsy total Gleason score": "BIOPSY TOTAL GLEASON SCORE",
    "Biopsy positive cores": "BX POSITIVE CORES",
    "Biopsy total cores": "BX TOTAL CORES",
    "Clinical stage": "CLINICAL STAGE",
    "Main treatment": "MAIN TREATMENT",
    "Adjuvant therapy": "ADJUVANT THERAPY TYPE",
    "Pathological Gleason score 1": "PATH GLEASON SCORE 1",
    "Pathological Gleason score 2": "PATH GLEASON SCORE 2",
    "Pathological total Gleason score": "PATH GLEASON TOTAL SCORE",
    "Surgical margins": "SURGICAL MARGINS",
    "Lymph node status": "LYMPH NODE STATUS",
    "Pathological stage": "PATHOLOGICAL STAGE",
    "Follow-up duration (months)": "LAST FOLLOW UP",
    "Status at follow-up": "STATUS AT FOLLOW UP",
    "Time to PSA recurrence (months)": "TIME TO PSA RECURRENCE",
    "Time to metastasis (months)": "TIME TO METASTASIS",
    "Metastasis location": "METASTASIS LOCATION",
    "Metastasis treatment": "METASTASIS TREATMENT",
    "Time to castration resistance (months)": "TIME TO CRPC",
    "Castration resistance treatment": "CASTRATION RESISTANCE TREATMENT",
}

# ─────────────────────────────────────────────────────────────────────────────
# 1.  Variable registry
# ─────────────────────────────────────────────────────────────────────────────

SKIP_VARS = {
    "DIAGNOSIS DATE",
    "DATE OF SURGERY/RADIATION",
    "DATE LAST SEEN",
    "LAST SEEN 4-D",
    "LAST FOLLOW UP DATE",  # date source → LAST FOLLOW UP (months) derived below
    "DATE OF PSA RECURRENCE",
    "DATE OF METASTASIS DX",
    "CASTRATION RESISTANCE DATE",
    "DATE OF DEATH",
    "ROBOT CASE NUMBER",
    "REASON NO FURTHER FOLLOW UP",
}

CONTINUOUS_VARS = {
    "AGE AT DX",
    "PSA AT DX",
    # "TESTOSTERONE AT DX",
    "BX POSITIVE CORES",
    "BX TOTAL CORES",
    "PCT POSITIVE CORES",
    "BIOPSY GLEASON SCORE 1",
    "BIOPSY GLEASON SCORE 2",
    "BIOPSY TOTAL GLEASON SCORE",
    "PATH GLEASON SCORE 1",
    "PATH GLEASON SCORE 2",
    "PATH GLEASON TOTAL SCORE",
    "MONTHS OF NEO-ADJUVANT THERAPY",
    "IAS NUMBER OF CYCLES",
    "PSA AT RECURRENCE",
    "TIME TO PSA RECURRENCE",
    "LAST FOLLOW UP",
    "TIME TO METASTASIS",  # derived from dates if available
}

CATEGORICAL_VARS = {
    "ETHNICITY (INFERRED)",
    "ANCESTRY GROUPED",  # derived: European / East Asian / Other
    "CLINICAL STAGE",
    "CLINICAL STAGE BINARY",  # derived: T1/T2 vs T3
    "MAIN TREATMENT",
    "SURGICAL APPROACH",  # derived: Open RP vs Robotic RP
    "NEO-ADJUVANT THERAPY DRUGS",
    "ADJUVANT THERAPY TYPE",
    # "IS THERE A TERTIARY GLEASON COMPONENT",
    # "IF TERTIARY GLEASON PATTERN PRESENT",
    "SURGICAL MARGINS",
    "LYMPH NODE STATUS",
    "PATHOLOGICAL STAGE",
    "PATH STAGE BINARY",  # derived: pT2 vs pT3
    "ALIVE AT FOLLOW UP",  # derived: Alive vs Deceased
    "STATUS AT FOLLOW UP",
    "PSA RECURRENCE TREATMENT",
    "METASTASIS LOCATION",
    "METASTASIS TREATMENT",
    "CASTRATION RESISTANCE TREATMENT",
}

MATCHED_VARS = {
    "PSA AT DX",
    "BIOPSY TOTAL GLEASON SCORE",
    "BX POSITIVE CORES",
    "BX TOTAL CORES",
    "PCT POSITIVE CORES",
    "SURGICAL MARGINS",
    "ETHNICITY (INFERRED)",
    "ANCESTRY GROUPED",  # derived from matched variable
}

METASTATIC_ONLY = {
    "PSA AT RECURRENCE",
    "TIME TO PSA RECURRENCE",
    "PSA RECURRENCE TREATMENT",
    "METASTASIS LOCATION",
    "METASTASIS TREATMENT",
    "CASTRATION RESISTANCE TREATMENT",
    "TIME TO METASTASIS",
}

# Table 1 variable order
TABLE1_VARS = [
    # Demographics
    "AGE AT DX",
    "ANCESTRY GROUPED",
    # Pre-treatment tumour characteristics
    "PSA AT DX",
    # "TESTOSTERONE AT DX",
    "BIOPSY TOTAL GLEASON SCORE",
    "PCT POSITIVE CORES",
    "CLINICAL STAGE BINARY",  # T1/T2 vs T3
    # Post-operative pathology (all patients had RP)
    "SURGICAL MARGINS",
    "LYMPH NODE STATUS",
    "PATH STAGE BINARY",  # pT2 vs pT3
    # Follow-up
    "LAST FOLLOW UP",
    "ALIVE AT FOLLOW UP",  # binary: Alive vs Deceased
]

SUPPLEMENTARY_EXTRA_VARS = [
    # Full status at follow-up breakdown (ungrouped)
    "STATUS AT FOLLOW UP",
    # Full ancestry breakdown (ungrouped)
    "ETHNICITY (INFERRED)",
    # Full T-stage subcategories (ungrouped)
    "CLINICAL STAGE",
    "PATHOLOGICAL STAGE",
    # Surgical approach
    "MAIN TREATMENT",
    "SURGICAL APPROACH",
    # Gleason subscores
    "BIOPSY GLEASON SCORE 1",
    "BIOPSY GLEASON SCORE 2",
    "BX POSITIVE CORES",
    "BX TOTAL CORES",
    # Pathological Gleason (expected to differ by design)
    "PATH GLEASON TOTAL SCORE",
    "PATH GLEASON SCORE 1",
    "PATH GLEASON SCORE 2",
    # Treatment detail
    "NEO-ADJUVANT THERAPY DRUGS",
    "MONTHS OF NEO-ADJUVANT THERAPY",
    "ADJUVANT THERAPY TYPE",
    "IAS NUMBER OF CYCLES",
    # Metastatic-arm outcomes
    "PSA AT RECURRENCE",
    "TIME TO PSA RECURRENCE",
    "PSA RECURRENCE TREATMENT",
    "TIME TO METASTASIS",
    "METASTASIS LOCATION",
    "METASTASIS TREATMENT",
    "CASTRATION RESISTANCE TREATMENT",
]

# ─────────────────────────────────────────────────────────────────────────────
# 1b. Publication-ready display labels (column name → table/figure label)
#     Edit these if your preferred wording differs.
# ─────────────────────────────────────────────────────────────────────────────

DISPLAY_NAMES = {
    "AGE AT DX": "Age at diagnosis, yr",
    "ETHNICITY (INFERRED)": "Ancestry",
    "PSA AT DX": "PSA at diagnosis, ng/mL",
    # "TESTOSTERONE AT DX": "Testosterone at diagnosis, nmol/L",
    "BIOPSY TOTAL GLEASON SCORE": "Biopsy Gleason score (total)",
    "BIOPSY GLEASON SCORE 1": "Biopsy Gleason primary pattern",
    "BIOPSY GLEASON SCORE 2": "Biopsy Gleason secondary pattern",
    "PCT POSITIVE CORES": "Positive biopsy cores, %",
    "BX POSITIVE CORES": "Positive biopsy cores, n",
    "BX TOTAL CORES": "Total biopsy cores, n",
    "ANCESTRY GROUPED": "Ancestry",
    "CLINICAL STAGE": "Clinical T-stage",
    "MAIN TREATMENT": "Primary treatment",
    # "IS THERE A TERTIARY GLEASON COMPONENT": "Tertiary Gleason component",
    # "IF TERTIARY GLEASON PATTERN PRESENT": "Tertiary Gleason pattern",
    "PATH GLEASON TOTAL SCORE": "Pathological Gleason score (total)",
    "PATH GLEASON SCORE 1": "Pathological Gleason primary pattern",
    "PATH GLEASON SCORE 2": "Pathological Gleason secondary pattern",
    "SURGICAL MARGINS": "Surgical margins",
    "LYMPH NODE STATUS": "Lymph node involvement",
    "PATHOLOGICAL STAGE": "Pathological T-stage (full)",
    "CLINICAL STAGE BINARY": "Clinical T-stage",
    "PATH STAGE BINARY": "Pathological T-stage",
    "SURGICAL APPROACH": "Surgical approach",
    "LAST FOLLOW UP": "Follow-up duration, months",
    "STATUS AT FOLLOW UP": "Status at last follow-up",
    "ALIVE AT FOLLOW UP": "Alive at last follow-up",
    "MONTHS OF NEO-ADJUVANT THERAPY": "Neoadjuvant therapy duration, months",
    "NEO-ADJUVANT THERAPY DRUGS": "Neoadjuvant therapy agent",
    "ADJUVANT THERAPY TYPE": "Adjuvant therapy type",
    "IAS NUMBER OF CYCLES": "IAS cycles, n",
    "PSA AT RECURRENCE": "PSA at biochemical recurrence, ng/mL",
    "TIME TO PSA RECURRENCE": "Time to PSA recurrence, months",
    "PSA RECURRENCE TREATMENT": "PSA recurrence treatment",
    "TIME TO METASTASIS": "Time to metastasis, months",
    "METASTASIS LOCATION": "Metastasis site",
    "METASTASIS TREATMENT": "Metastasis treatment",
    "CASTRATION RESISTANCE TREATMENT": "Castration-resistant PCa treatment",
}


def display(var):
    """Return publication-ready label; fall back to title-case column name."""
    return DISPLAY_NAMES.get(var, var.title())


# Variables shown in love plot (must be computable as continuous or binary SMD)
LOVE_VARS = [
    "AGE AT DX",
    "PSA AT DX",
    # "TESTOSTERONE AT DX",
    "BIOPSY TOTAL GLEASON SCORE",
    "PCT POSITIVE CORES",
    "LAST FOLLOW UP",
    "ANCESTRY GROUPED",  # 3-cat: European / East Asian / Other
    "CLINICAL STAGE BINARY",  # binary: T1/T2 vs T3
    "SURGICAL MARGINS",  # binary
    "LYMPH NODE STATUS",  # binary
    "PATH STAGE BINARY",  # binary: pT2 vs pT3
    "IS THERE A TERTIARY GLEASON COMPONENT",  # binary
]


# ─────────────────────────────────────────────────────────────────────────────
# 2.  Type classifier
# ─────────────────────────────────────────────────────────────────────────────


def classify(var, df):
    if var in SKIP_VARS:
        return "skip"
    if var in METASTATIC_ONLY:
        return "metonly"
    if var in CONTINUOUS_VARS:
        return "continuous"
    if var in CATEGORICAL_VARS:
        return "categorical"
    if var not in df.columns:
        return "skip"
    col = pd.to_numeric(df[var], errors="coerce")
    if col.notna().sum() > 0 and col.nunique() > 6:
        return "continuous"
    return "categorical"


# ─────────────────────────────────────────────────────────────────────────────
# 3.  Formatting helpers
# ─────────────────────────────────────────────────────────────────────────────


def fmt_median_iqr(series):
    s = pd.to_numeric(series, errors="coerce").dropna()
    if len(s) == 0:
        return "—"
    med, q1, q3 = np.median(s), np.percentile(s, 25), np.percentile(s, 75)
    note = f" [n={len(s)}]" if series.isna().any() else ""
    return f"{med:.1f} ({q1:.1f}–{q3:.1f}){note}"


def fmt_n_pct(count, total):
    return "—" if total == 0 else f"{int(count)} ({100 * count / total:.1f}%)"


def fmt_p(p):
    if p is None or (isinstance(p, float) and np.isnan(p)):
        return "—"
    if p < 0.001:
        return "<0.001"
    return f"{p:.3f}"


# ─────────────────────────────────────────────────────────────────────────────
# 4.  Statistical tests
# ─────────────────────────────────────────────────────────────────────────────


def mw_test(a, b):
    a = pd.to_numeric(a, errors="coerce").dropna()
    b = pd.to_numeric(b, errors="coerce").dropna()
    if len(a) < 2 or len(b) < 2:
        return np.nan
    return stats.mannwhitneyu(a, b, alternative="two-sided").pvalue


def cat_test(a, b):
    a = a.astype(str).replace("nan", np.nan)
    b = b.astype(str).replace("nan", np.nan)
    cats = sorted(set(a.dropna()) | set(b.dropna()), key=str)
    if not cats:
        return np.nan, "—"
    ct = (
        pd.DataFrame(
            {"M": a.value_counts(dropna=True), "N": b.value_counts(dropna=True)},
            index=cats,
        )
        .fillna(0)
        .astype(int)
    )
    if ct.shape == (2, 2):
        return stats.fisher_exact(ct.values).pvalue, "Fisher's exact"
    _, p, _, exp = stats.chi2_contingency(ct.values)
    warn = " *" if (exp < 5).any() else ""
    return p, f"Chi-square{warn}"


# ─────────────────────────────────────────────────────────────────────────────
# 5.  Table row builders
# ─────────────────────────────────────────────────────────────────────────────

COLS = [
    "Variable",
    "Metastatic (n=26)",
    "Non-metastatic (n=26)",
    "p-value",
    "Statistical test",
]


def rows_continuous(var, met, nmet):
    flag = "†" if var in MATCHED_VARS else ""
    m = met[var] if var in met.columns else pd.Series(dtype=float)
    n = nmet[var] if var in nmet.columns else pd.Series(dtype=float)
    p = mw_test(m, n)
    return [
        {
            "Variable": f"{display(var)}{flag},  median (IQR)",
            "Metastatic (n=26)": fmt_median_iqr(m),
            "Non-metastatic (n=26)": fmt_median_iqr(n),
            "p-value": fmt_p(p),
            "Statistical test": "Mann-Whitney U",
        }
    ]


def rows_categorical(var, met, nmet):
    flag = "†" if var in MATCHED_VARS else ""
    m = (
        met[var].astype(str).replace("nan", np.nan)
        if var in met.columns
        else pd.Series(dtype=object)
    )
    n = (
        nmet[var].astype(str).replace("nan", np.nan)
        if var in nmet.columns
        else pd.Series(dtype=object)
    )
    cats = sorted(set(m.dropna()) | set(n.dropna()), key=str)
    print(var)
    print(cats)
    if not cats:
        return []
    p, test = cat_test(m, n)
    mt, nt = m.notna().sum(), n.notna().sum()
    rows = [
        {
            "Variable": f"{display(var)}{flag},  n (%)",
            "Metastatic (n=26)": f"[n={mt}]" if mt < 26 else "",
            "Non-metastatic (n=26)": f"[n={nt}]" if nt < 26 else "",
            "p-value": fmt_p(p),
            "Statistical test": test,
        }
    ]
    for cat in cats:
        rows.append(
            {
                "Variable": f"    {cat}",
                "Metastatic (n=26)": fmt_n_pct((m == cat).sum(), mt),
                "Non-metastatic (n=26)": fmt_n_pct((n == cat).sum(), nt),
                "p-value": "",
                "Statistical test": "",
            }
        )
    return rows


def rows_metonly(var, met):
    if var not in met.columns:
        return []
    s = met[var]
    if s.notna().sum() == 0:
        return []
    if var in CONTINUOUS_VARS:
        return [
            {
                "Variable": f"{display(var)}  (metastatic arm),  median (IQR)",
                "Metastatic (n=26)": fmt_median_iqr(pd.to_numeric(s, errors="coerce")),
                "Non-metastatic (n=26)": "n/a",
                "p-value": "—",
                "Statistical test": "Descriptive only",
            }
        ]
    s2 = s.astype(str).replace("nan", np.nan)
    cats = sorted(s2.dropna().unique(), key=str)
    total = s2.notna().sum()
    rows = [
        {
            "Variable": f"{display(var)}  (metastatic arm),  n (%)",
            "Metastatic (n=26)": f"[n={total}]" if total < 26 else "",
            "Non-metastatic (n=26)": "n/a",
            "p-value": "—",
            "Statistical test": "Descriptive only",
        }
    ]
    for cat in cats:
        rows.append(
            {
                "Variable": f"    {cat}",
                "Metastatic (n=26)": fmt_n_pct((s2 == cat).sum(), total),
                "Non-metastatic (n=26)": "n/a",
                "p-value": "",
                "Statistical test": "",
            }
        )
    return rows


def build_table(var_list, df, met, nmet, include_metonly=False):
    rows = []
    for var in var_list:
        vtype = classify(var, df)
        if vtype == "skip":
            continue
        if vtype == "metonly":
            if include_metonly:
                rows.extend(rows_metonly(var, met))
            continue
        if vtype == "continuous":
            rows.extend(rows_continuous(var, met, nmet))
        else:
            rows.extend(rows_categorical(var, met, nmet))
        print("Finish " + var)
    return pd.DataFrame(rows, columns=COLS)


# ─────────────────────────────────────────────────────────────────────────────
# 6.  Love plot
# ─────────────────────────────────────────────────────────────────────────────


def smd_continuous(a, b):
    """Standardized mean difference: (mean_a − mean_b) / pooled SD."""
    a = pd.to_numeric(a, errors="coerce").dropna()
    b = pd.to_numeric(b, errors="coerce").dropna()
    if len(a) < 2 or len(b) < 2:
        return np.nan
    pool_var = (a.var(ddof=1) + b.var(ddof=1)) / 2
    return 0.0 if pool_var == 0 else (a.mean() - b.mean()) / np.sqrt(pool_var)


def smd_binary(a, b, cat):
    """SMD for a single binary category: proportion difference / pooled proportion SD."""
    a = a.astype(str).replace("nan", np.nan)
    b = b.astype(str).replace("nan", np.nan)
    pm = (a.dropna() == cat).mean()
    pn = (b.dropna() == cat).mean()
    pp = (pm + pn) / 2
    return np.nan if pp in (0, 1) else (pm - pn) / np.sqrt(pp * (1 - pp))


def plot_love(met, nmet, df, output_path):
    MET_COLOR = "#c0392b"
    NMET_COLOR = "#2980b9"
    records = []  # (label, smd, is_matched)

    for var in LOVE_VARS:
        if var not in df.columns:
            continue
        vtype = classify(var, df)
        matched = var in MATCHED_VARS

        if vtype == "continuous":
            s = smd_continuous(
                met[var] if var in met.columns else pd.Series(dtype=float),
                nmet[var] if var in nmet.columns else pd.Series(dtype=float),
            )
            if not np.isnan(s):
                records.append((display(var), s, matched))

        elif vtype == "categorical":
            ms = (
                met[var].astype(str).replace("nan", np.nan)
                if var in met.columns
                else pd.Series(dtype=object)
            )
            ns = (
                nmet[var].astype(str).replace("nan", np.nan)
                if var in nmet.columns
                else pd.Series(dtype=object)
            )
            cats = sorted(set(ms.dropna()) | set(ns.dropna()), key=str)
            if len(cats) == 2:
                s = smd_binary(ms, ns, cats[1])
                if not np.isnan(s):
                    records.append((display(var), s, matched))
            else:
                # Multi-category: one dot per category (all shown for balance assessment)
                for cat in cats:
                    s = smd_binary(ms, ns, cat)
                    if not np.isnan(s):
                        records.append((f"  {display(var)}: {cat}", s, matched))

    if not records:
        print("Warning: no variables available for love plot — skipping.")
        return

    # Sort by |SMD| descending (most imbalanced at top)
    records.sort(key=lambda x: abs(x[1]), reverse=True)
    labels = [r[0] for r in records]
    smds = [r[1] for r in records]
    colors = [MET_COLOR if r[2] else NMET_COLOR for r in records]

    fig_h = max(5, len(records) * 0.45)
    fig, ax = plt.subplots(figsize=(8, fig_h))
    y = list(range(len(records)))

    ax.scatter(smds, y, c=colors, s=70, zorder=4, edgecolors="white", linewidths=0.5)

    # Reference lines
    ax.axvline(0, color="black", linewidth=0.9, zorder=3)
    ax.axvline(0.10, color="gray", linewidth=0.6, linestyle="--", zorder=2)
    ax.axvline(-0.10, color="gray", linewidth=0.6, linestyle="--", zorder=2)
    ax.axvline(0.25, color="darkorange", linewidth=0.6, linestyle=":", zorder=2)
    ax.axvline(-0.25, color="darkorange", linewidth=0.6, linestyle=":", zorder=2)
    ax.axvspan(-0.10, 0.10, alpha=0.07, color="green", zorder=1)
    ax.axvspan(0.10, 0.25, alpha=0.05, color="darkorange", zorder=1)
    ax.axvspan(-0.25, -0.10, alpha=0.05, color="darkorange", zorder=1)

    ax.set_yticks(y)
    ax.set_yticklabels(labels, fontsize=8)
    ax.set_xlabel(
        "Standardized Mean Difference (SMD)\nPositive = higher in metastatic arm",
        fontsize=9,
    )
    ax.set_title(
        "Covariate Balance: Metastatic vs Non-metastatic",
        fontsize=10,
        fontweight="bold",
        loc="center",
    )

    legend_elements = [
        Line2D(
            [0],
            [0],
            marker="o",
            color="w",
            markerfacecolor=MET_COLOR,
            markersize=9,
            label="Matching criterion (†)",
        ),
        Line2D(
            [0],
            [0],
            marker="o",
            color="w",
            markerfacecolor=NMET_COLOR,
            markersize=9,
            label="Unmatched covariate",
        ),
        mpatches.Patch(facecolor="green", alpha=0.15, label="|SMD| < 0.10  (balanced)"),
        mpatches.Patch(
            facecolor="darkorange", alpha=0.15, label="|SMD| 0.10–0.25  (borderline)"
        ),
        mpatches.Patch(
            facecolor="darkorange", alpha=0.40, label="|SMD| > 0.25  (imbalanced)"
        ),
    ]
    ax.legend(handles=legend_elements, fontsize=10, loc="upper right", framealpha=0.9)
    ax.spines[["top", "right"]].set_visible(False)

    plt.tight_layout()
    plt.savefig(output_path, dpi=300, bbox_inches="tight")
    plt.close()
    print(f"Saved → {output_path}")


# ─────────────────────────────────────────────────────────────────────────────
# 7.  Swimmer plot
# ─────────────────────────────────────────────────────────────────────────────


def date_diff_months(df, date_col, ref_col="DIAGNOSIS DATE"):
    """Return months between ref_col and date_col. Returns None if columns absent."""
    if date_col not in df.columns or ref_col not in df.columns:
        return None
    try:
        ref = pd.to_datetime(df[ref_col].reset_index(drop=True), errors="coerce")
        evt = pd.to_datetime(df[date_col].reset_index(drop=True), errors="coerce")
        return ((evt - ref).dt.days / 30.44).where(evt.notna() & ref.notna())
    except Exception:
        return None


def arm_records(arm_df, arm_label, color):
    """Extract one dict per patient for swimmer plot."""
    arm_df = arm_df.reset_index(drop=True)
    followup = pd.to_numeric(arm_df.get("LAST FOLLOW UP", pd.Series()), errors="coerce")
    # followup = date_diff_months(arm_df, "LAST FOLLOW UP", ref_col="DIAGNOSIS DATE")
    psa_t = (
        pd.to_numeric(
            arm_df.get("TIME TO PSA RECURRENCE", pd.Series()), errors="coerce"
        )
        if "TIME TO PSA RECURRENCE" in arm_df.columns
        else None
    )
    # Time-to-metastasis is now provided pre-derived as TIME TO METASTASIS
    # (months from diagnosis); previously computed from DIAGNOSIS DATE +
    # DATE OF METASTASIS DX, which are no longer in the de-identified TSV.
    met_t = (
        pd.to_numeric(arm_df.get("TIME TO METASTASIS", pd.Series()), errors="coerce")
        if arm_label == "Metastatic" and "TIME TO METASTASIS" in arm_df.columns
        else None
    )

    dead = None
    if "STATUS AT FOLLOW UP" in arm_df.columns:
        dead = (
            arm_df["STATUS AT FOLLOW UP"]
            .astype(str)
            .str.lower()
            .str.contains("dead|died|death|deceased", na=False)
        )

    records = []
    for i in range(len(arm_df)):
        fu = followup.iloc[i] if i < len(followup) else np.nan
        if pd.isna(fu):
            continue
        records.append(
            {
                "arm": arm_label,
                "color": color,
                "followup": fu,
                "psa": psa_t.iloc[i]
                if psa_t is not None and i < len(psa_t)
                else np.nan,
                "metastasis": met_t.iloc[i]
                if met_t is not None and i < len(met_t)
                else np.nan,
                "dead": bool(dead.iloc[i])
                if dead is not None and i < len(dead)
                else False,
            }
        )
    return records


def plot_swimmer(met, nmet, output_path):
    # ── TOGGLE ────────────────────────────────────────────────────────────────
    # Set to True to draw a dashed vertical line at the median PSA recurrence
    # time for the non-metastatic arm; False to hide it.
    SHOW_NMET_PSA_MEDIAN = True  # <── change here

    MET_COLOR = "#c0392b"
    NMET_COLOR = "#2980b9"

    if "LAST FOLLOW UP" not in pd.concat([met, nmet]).columns:
        print("Warning: LAST FOLLOW UP column missing — skipping swimmer plot.")
        return

    # Build records: metastatic first (sorted ascending by FU), then non-met
    met_rec = sorted(
        arm_records(met, "Metastatic", MET_COLOR), key=lambda r: r["followup"]
    )
    nmet_rec = sorted(
        arm_records(nmet, "Non-metastatic", NMET_COLOR), key=lambda r: r["followup"]
    )
    all_rec = met_rec + nmet_rec
    n = len(all_rec)

    if n == 0:
        print("Warning: no valid follow-up data for swimmer plot.")
        return

    # ── Patient ID extraction ─────────────────────────────────────────────────
    # Replicates the NaN-followup filtering in arm_records to recover which
    # rows (and therefore which patient IDs) were included, in the same order.
    pid_col = met.columns[0]  # first column = patient ID (UNNAMED: 0)

    def _pid_fu_pairs(arm_df):
        """Return (patient_id, followup_months) pairs, NaN-followup rows skipped,
        in the same iteration order as arm_records.

        LAST FOLLOW UP is now provided directly as months in the new TSV;
        previously it was derived from LAST FOLLOW UP DATE - DIAGNOSIS DATE.
        """
        fu_series = pd.to_numeric(
            arm_df.get("LAST FOLLOW UP", pd.Series()), errors="coerce"
        )
        pairs = []
        for i in range(len(arm_df)):
            fu = fu_series.iloc[i] if i < len(fu_series) else np.nan
            if not pd.isna(fu):
                pairs.append((arm_df[pid_col].iloc[i], fu))
        return pairs

    met_pid_pairs = sorted(_pid_fu_pairs(met), key=lambda x: x[1])
    nmet_pid_pairs = sorted(_pid_fu_pairs(nmet), key=lambda x: x[1])
    all_pids = [p[0] for p in met_pid_pairs] + [p[0] for p in nmet_pid_pairs]

    fig_h = max(6, n * 0.32)
    fig, ax = plt.subplots(figsize=(11, fig_h))

    for y_pos, r in enumerate(all_rec):
        ax.barh(
            y_pos,
            r["followup"],
            height=0.65,
            color=r["color"],
            alpha=0.70,
            edgecolor="none",
        )

        # End-of-bar marker
        marker = "x" if r["dead"] else "|"
        mcolor = "black" if r["dead"] else r["color"]
        ax.scatter(
            r["followup"],
            y_pos,
            marker=marker,
            color=mcolor,
            s=40,
            zorder=5,
            linewidths=1.5,
        )

        # PSA recurrence
        if not pd.isna(r["psa"]) and 0 < r["psa"] <= r["followup"]:
            ax.scatter(r["psa"], y_pos, marker="^", color="darkorange", s=55, zorder=6)

        # Metastasis event
        if not pd.isna(r["metastasis"]) and 0 < r["metastasis"] <= r["followup"]:
            ax.scatter(
                r["metastasis"], y_pos, marker="*", color="darkred", s=90, zorder=6
            )

    # Arm separator and labels
    n_met = len(met_rec)
    if 0 < n_met < n:
        ax.axhline(n_met - 0.5, color="gray", linewidth=0.8, linestyle=":")

    max_fu = max(r["followup"] for r in all_rec)
    label_x = -max_fu * 0.03

    if n_met > 0:
        ax.text(
            label_x,
            n_met / 2 - 0.5,
            f"Metastatic\n(n={n_met})",
            ha="right",
            va="center",
            fontsize=10,
            color=MET_COLOR,
            fontweight="bold",
        )
    n_nmet = n - n_met
    if n_nmet > 0:
        ax.text(
            label_x,
            n_met + n_nmet / 2 - 0.5,
            f"Non-metastatic\n(n={n_nmet})",
            ha="right",
            va="center",
            fontsize=10,
            color=NMET_COLOR,
            fontweight="bold",
        )

    # ── Median annotation lines (arm-specific vertical spans) ────────────────
    met_ymin = -0.45
    met_ymax = n_met - 0.55
    nmet_ymin = n_met - 0.45
    nmet_ymax = n - 0.55

    # Non-met: median follow-up
    nmet_fu_vals = [r["followup"] for r in nmet_rec if not pd.isna(r["followup"])]
    if nmet_fu_vals:
        med_fu = np.median(nmet_fu_vals)
        ax.vlines(
            med_fu,
            nmet_ymin,
            nmet_ymax,
            color=NMET_COLOR,
            linewidth=1.4,
            linestyle="--",
            zorder=7,
        )
        ax.text(
            med_fu,
            nmet_ymax + 0.25,
            f"Median FU\n{med_fu:.0f} mo",
            ha="center",
            va="bottom",
            fontsize=10,
            color=NMET_COLOR,
            fontweight="bold",
        )

    # Non-met: median PSA recurrence (optional toggle)
    nmet_psa_vals = [r["psa"] for r in nmet_rec if not pd.isna(r["psa"])]
    if SHOW_NMET_PSA_MEDIAN and nmet_psa_vals:
        med_nmet_psa = np.median(nmet_psa_vals)
        ax.vlines(
            med_nmet_psa,
            nmet_ymin,
            nmet_ymax,
            color="darkorange",
            linewidth=1.4,
            linestyle="-.",
            zorder=7,
        )
        ax.text(
            med_nmet_psa,
            nmet_ymax + 0.25,  # nmet_ymin - 0.35,
            f"Median PSA (non-met)\n{med_nmet_psa:.0f} mo",
            ha="center",
            va="bottom",  # va="top",
            fontsize=10,
            color="darkorange",
            fontweight="bold",
        )

    # Met: median PSA recurrence time
    met_psa_vals = [r["psa"] for r in met_rec if not pd.isna(r["psa"])]
    if met_psa_vals:
        med_psa = np.median(met_psa_vals)
        ax.vlines(
            med_psa,
            met_ymin,
            met_ymax,
            color="darkorange",
            linewidth=1.4,
            linestyle="--",
            zorder=7,
        )
        ax.text(
            med_psa - 20,
            met_ymin - 0.25,  # met_ymax + 0.25,
            f"Median PSA\n{med_psa:.0f} mo",
            ha="left",
            va="top",  # va="bottom",
            fontsize=10,
            color="darkorange",
            fontweight="bold",
        )

    # Met: median time to metastasis
    met_mets_vals = [r["metastasis"] for r in met_rec if not pd.isna(r["metastasis"])]

    # ── Time-dependent bias verification ─────────────────────────────────────
    # Mann-Whitney: non-met follow-up > met time-to-metastasis (one-sided)
    if nmet_fu_vals and met_mets_vals:
        stat, p_mw = stats.mannwhitneyu(
            nmet_fu_vals, met_mets_vals, alternative="greater"
        )
        p_str = "<0.001" if p_mw < 0.001 else f"{p_mw:.3f}"
        print()
        print("── Time-dependent bias check (swimmer plot) ─────────────────────────")
        print(
            f"  Non-met  median follow-up:         {np.median(nmet_fu_vals):.1f} months  (n={len(nmet_fu_vals)})"
        )
        print(
            f"  Met      median time-to-metastasis: {np.median(met_mets_vals):.1f} months  (n={len(met_mets_vals)})"
        )
        print(f"  Mann-Whitney U = {stat:.0f},  p (one-sided, non-met > met) = {p_str}")
        print("─────────────────────────────────────────────────────────────────────")

        # Stats annotation box on the plot
        box_text = (
            f"Non-met median FU:   {np.median(nmet_fu_vals):.0f} mo\n"
            f"Met median mets:      {np.median(met_mets_vals):.0f} mo\n"
            f"Mann-Whitney p {p_str}"
        )
        ax.text(
            0.98,
            0.8,  # 0.98,
            box_text,
            transform=ax.transAxes,
            fontsize=16,  # 7.5,
            va="top",
            ha="right",
            linespacing=1.6,
            bbox=dict(
                boxstyle="round,pad=0.4", facecolor="white", edgecolor="gray", alpha=0.9
            ),
        )

    if met_mets_vals:
        med_mets = np.median(met_mets_vals)
        ax.vlines(
            med_mets,
            met_ymin,
            met_ymax,
            color="darkred",
            linewidth=1.4,
            linestyle="--",
            zorder=7,
        )
        ax.text(
            med_mets + 5,
            met_ymin - 0.25,  # met_ymax + 0.25,
            f"Median mets\n{med_mets:.0f} mo",
            ha="center",
            va="top",  # va="bottom",
            fontsize=10,
            color="darkred",
            fontweight="bold",
        )

    ax.set_xlabel("Months from diagnosis", fontsize=10)

    # ── Patient ID y-axis labels ──────────────────────────────────────────────
    ax.set_yticks(range(n))
    ax.set_yticklabels(all_pids, fontsize=7)
    ax.set_ylabel("Patient ID", fontsize=9)

    ax.set_title("Patient Follow-up Timeline", fontsize=11, fontweight="bold")
    ax.set_xlim(left=0)

    legend_elements = [
        mpatches.Patch(facecolor=MET_COLOR, alpha=0.7, label="Metastatic arm"),
        mpatches.Patch(facecolor=NMET_COLOR, alpha=0.7, label="Non-metastatic arm"),
        Line2D(
            [0],
            [0],
            marker="^",
            color="w",
            markerfacecolor="darkorange",
            markersize=9,
            label="PSA recurrence",
        ),
        Line2D(
            [0],
            [0],
            marker="*",
            color="w",
            markerfacecolor="darkred",
            markersize=11,
            label="Metastasis (from date)",
        ),
        Line2D(
            [0],
            [0],
            marker="x",
            color="black",
            markersize=7,
            linewidth=0,
            label="Death",
        ),
        Line2D(
            [0],
            [0],
            marker="|",
            color="gray",
            markersize=7,
            linewidth=0,
            label="Alive at last follow-up",
        ),
        Line2D(
            [0],
            [0],
            color=NMET_COLOR,
            linewidth=1.4,
            linestyle="--",
            label="Median follow-up (non-met)",
        ),
        Line2D(
            [0],
            [0],
            color="darkorange",
            linewidth=1.4,
            linestyle="--",
            label="Median PSA recurrence (met)",
        ),
        Line2D(
            [0],
            [0],
            color="darkred",
            linewidth=1.4,
            linestyle="--",
            label="Median time to metastasis (met)",
        ),
    ]
    if SHOW_NMET_PSA_MEDIAN and nmet_psa_vals:
        legend_elements.append(
            Line2D(
                [0],
                [0],
                color="darkorange",
                linewidth=1.4,
                linestyle="-.",
                label="Median PSA recurrence (non-met)",
            )
        )
    ax.legend(handles=legend_elements, fontsize=16, loc="lower right", framealpha=0.9)
    ax.spines[["top", "right", "left"]].set_visible(False)

    plt.tight_layout()
    plt.savefig(output_path, dpi=300, bbox_inches="tight")
    plt.close()
    print(f"Saved → {output_path}")


# ─────────────────────────────────────────────────────────────────────────────
# 8.  Excel writer
# ─────────────────────────────────────────────────────────────────────────────


def plot_table1_pdf(t1_df, output_path):
    """
    Render Table 1 as a publication-ready PDF using matplotlib.
    Columns: Variable | Metastatic (n=26) | Non-metastatic (n=26) | p-value
    Row types detected automatically:
      - continuous    : Variable contains 'median (IQR)'
      - categorical header : Variable contains 'n (%)'  → bold, p-value shown
      - subcategory   : Variable starts with spaces     → indented
    """
    SHOW_COLS = ["Variable", "Metastatic (n=26)", "Non-metastatic (n=26)", "p-value"]
    df = t1_df[SHOW_COLS].fillna("").copy()
    n_rows = len(df)

    # ── Figure dimensions ─────────────────────────────────────────────────────
    ROW_H = 0.27  # inches per data row
    HDR_H = 0.55  # header block height
    FOOT_H = 0.65  # footnote block height
    FIG_W = 9.0
    FIG_H = HDR_H + n_rows * ROW_H + FOOT_H + 0.20

    fig, ax = plt.subplots(figsize=(FIG_W, FIG_H))
    ax.set_xlim(0, FIG_W)
    ax.set_ylim(0, FIG_H)
    ax.axis("off")

    # ── Column x positions (inches from left) and alignments ──────────────────
    COL_X = [0.15, 5.10, 6.80, 8.30]
    COL_HA = ["left", "center", "center", "center"]
    COL_HDR = [
        "Variable",
        "Metastatic\n(n=26)",
        "Non-metastatic\n(n=26)",
        "p-value",
    ]

    # ── Title ─────────────────────────────────────────────────────────────────
    y_title = FIG_H - 0.12
    ax.text(
        COL_X[0],
        y_title,
        "Table 1.  Clinical and pathological characteristics of the EPC cohort",
        ha="left",
        va="top",
        fontsize=10,
        fontweight="bold",
    )

    # ── Borders and header ────────────────────────────────────────────────────
    y_top = FIG_H - 0.28
    ax.axhline(y_top, color="black", linewidth=1.4, xmin=0.015, xmax=0.985)

    y_hdr_mid = y_top - HDR_H / 2
    for hdr, x, ha in zip(COL_HDR, COL_X, COL_HA):
        ax.text(
            x,
            y_hdr_mid,
            hdr,
            ha=ha,
            va="center",
            fontsize=9,
            fontweight="bold",
            multialignment="center",
        )

    y_after_hdr = y_top - HDR_H
    ax.axhline(y_after_hdr, color="black", linewidth=0.9, xmin=0.015, xmax=0.985)

    # ── Data rows ─────────────────────────────────────────────────────────────
    for i, (_, row) in enumerate(df.iterrows()):
        y = y_after_hdr - (i + 0.5) * ROW_H

        var_raw = str(row["Variable"])
        is_sub = var_raw.startswith("    ")  # subcategory value
        is_cat = (not is_sub) and ("n (%)" in var_raw)  # categorical header
        is_cont = (not is_sub) and ("median (IQR)" in var_raw)

        label = var_raw.strip()
        x_label = COL_X[0] + (0.30 if is_sub else 0)
        fs = 8.5
        fw = "bold" if is_cat else "normal"
        fc = "#333333"

        ax.text(
            x_label,
            y,
            label,
            ha="left",
            va="center",
            fontsize=fs,
            fontweight=fw,
            color=fc,
        )

        for col, x, ha in zip(
            ["Metastatic (n=26)", "Non-metastatic (n=26)", "p-value"],
            COL_X[1:],
            COL_HA[1:],
        ):
            val = str(row[col]).strip()
            if val and val.lower() not in ("nan", ""):
                ax.text(x, y, val, ha=ha, va="center", fontsize=fs, color=fc)

    # ── Bottom border ─────────────────────────────────────────────────────────
    y_bottom = y_after_hdr - n_rows * ROW_H
    ax.axhline(y_bottom, color="black", linewidth=1.4, xmin=0.015, xmax=0.985)

    # ── Footnotes ─────────────────────────────────────────────────────────────
    footnotes = [
        "† Matching criterion; p-values shown for completeness.",
        "[n=X] after median (IQR): X patients had available data (when <26).",
        "Fisher's exact: binary/2-category variables.  "
        "Chi-square: >2 categories (* expected cell count <5 in ≥1 cell).",
        "Abbreviations: EPC, extreme phenotype cohort; IQR, interquartile range.",
    ]
    for k, line in enumerate(footnotes):
        ax.text(
            COL_X[0],
            y_bottom - 0.10 - k * 0.16,  # 0.12,
            line,
            ha="left",
            va="top",
            fontsize=10,  # 12,
            style="italic",
            color="#555555",
        )

    plt.savefig(output_path, dpi=300, bbox_inches="tight")
    plt.close()
    print(f"Saved → {output_path}")


def write_excel(tables, path):
    with pd.ExcelWriter(path, engine="openpyxl") as writer:
        for sheet, df in tables:
            df.to_excel(writer, sheet_name=sheet, index=False)
            ws = writer.sheets[sheet]
            for col_cells in ws.columns:
                width = max(
                    (len(str(c.value)) for c in col_cells if c.value), default=10
                )
                ws.column_dimensions[col_cells[0].column_letter].width = min(
                    width + 4, 70
                )
    print(f"Saved → {path}")


# ─────────────────────────────────────────────────────────────────────────────
# 9.  Main
# ─────────────────────────────────────────────────────────────────────────────


def main():
    parser = argparse.ArgumentParser(
        description="Generate Table 1, Supplementary, Swimmer plot, and Love plot."
    )
    parser.add_argument("input_tsv", help="Input TSV file path")
    parser.add_argument(
        "-o",
        "--outdir",
        default=".",
        help="Output directory (default: current directory)",
    )
    args = parser.parse_args()

    os.makedirs(args.outdir, exist_ok=True)
    stem = os.path.join(args.outdir, "table1")

    # ── Load ──────────────────────────────────────────────────────────────────
    try:
        df = pd.read_csv(args.input_tsv, sep="\t", dtype=str)
    except FileNotFoundError:
        sys.exit(f"ERROR: file not found — {args.input_tsv}")

    # Strip whitespace from header cells, then translate sentence-case headers
    # into the ALL-CAPS constants used throughout this script. Any column not
    # in COLUMN_RENAME is uppercased as a fallback (preserves the original
    # case-insensitive matching behaviour).
    df.columns = [c.strip() for c in df.columns]
    df = df.rename(columns=COLUMN_RENAME)
    df.columns = [
        c.upper() if c not in COLUMN_RENAME.values() else c for c in df.columns
    ]

    arm_col = df.columns[1]
    df[arm_col] = pd.to_numeric(df[arm_col], errors="coerce")
    met = df[df[arm_col] == 1].reset_index(drop=True)
    nmet = df[df[arm_col] == 0].reset_index(drop=True)
    print(
        f"Loaded {len(df)} patients: {len(met)} metastatic, {len(nmet)} non-metastatic."
    )

    # ── Derived variables ─────────────────────────────────────────────────────
    # % positive cores
    if "BX POSITIVE CORES" in df.columns and "BX TOTAL CORES" in df.columns:
        for frame in [df, met, nmet]:
            pos = pd.to_numeric(frame["BX POSITIVE CORES"], errors="coerce")
            tot = pd.to_numeric(frame["BX TOTAL CORES"], errors="coerce")
            frame["PCT POSITIVE CORES"] = (pos / tot * 100).where(tot > 0)
        print("Derived: PCT POSITIVE CORES")

    # Alive at last follow-up: binary (Alive vs Deceased)
    if "STATUS AT FOLLOW UP" in df.columns:
        for frame in [df, met, nmet]:
            frame["ALIVE AT FOLLOW UP"] = (
                frame["STATUS AT FOLLOW UP"]
                .astype(str)
                .str.upper()
                .str.strip()
                .map(
                    lambda x: (
                        "Alive"
                        if any(w in x for w in ["ALIVE", "LIVING"])
                        else "Deceased"
                        if any(w in x for w in ["DEAD", "DIED", "DEATH", "DECEASED"])
                        else np.nan
                    )
                )
            )
        print("Derived: ALIVE AT FOLLOW UP (Alive vs Deceased)")

    # Follow-up duration: compute from LAST FOLLOW UP DATE - DIAGNOSIS DATE
    if "LAST FOLLOW UP DATE" in df.columns and "DIAGNOSIS DATE" in df.columns:
        for frame in [df, met, nmet]:
            t = date_diff_months(frame, "LAST FOLLOW UP DATE")
            if t is not None:
                frame["LAST FOLLOW UP"] = t
        print(
            "Derived: LAST FOLLOW UP (months from DIAGNOSIS DATE to LAST FOLLOW UP DATE)"
        )

    # Ancestry: European / East Asian / Other
    if "ETHNICITY (INFERRED)" in df.columns:
        for frame in [df, met, nmet]:
            frame["ANCESTRY GROUPED"] = (
                frame["ETHNICITY (INFERRED)"]
                .astype(str)
                .str.upper()
                .str.strip()
                .map(
                    lambda x: (
                        "European"
                        if "EUROPEAN" in x
                        else "East Asian"
                        if "EAST ASIAN" in x
                        else "Other"
                        if x not in ("NAN", "")
                        else np.nan
                    )
                )
            )
        print("Derived: ANCESTRY GROUPED (European / East Asian / Other)")

    # Clinical T-stage: collapse to T1/T2 vs T3
    if "CLINICAL STAGE" in df.columns:
        for frame in [df, met, nmet]:
            frame["CLINICAL STAGE BINARY"] = (
                frame["CLINICAL STAGE"]
                .astype(str)
                .str.upper()
                .str.strip()
                .map(
                    lambda x: (
                        "T3"
                        if "T3" in x
                        else "T1/T2"
                        if any(t in x for t in ["T1", "T2"])
                        else np.nan
                    )
                )
            )
        print("Derived: CLINICAL STAGE BINARY (T1/T2 vs T3)")

    # Pathological T-stage: collapse to pT2 vs pT3
    if "PATHOLOGICAL STAGE" in df.columns:
        for frame in [df, met, nmet]:
            frame["PATH STAGE BINARY"] = (
                frame["PATHOLOGICAL STAGE"]
                .astype(str)
                .str.upper()
                .str.strip()
                .map(lambda x: "pT3" if "T3" in x else "pT2" if "T2" in x else np.nan)
            )
        print("Derived: PATH STAGE BINARY (pT2 vs pT3)")

    # Surgical approach: Open RP vs Robotic RP
    if "MAIN TREATMENT" in df.columns:
        for frame in [df, met, nmet]:
            frame["SURGICAL APPROACH"] = (
                frame["MAIN TREATMENT"]
                .astype(str)
                .str.upper()
                .str.strip()
                .map(
                    lambda x: (
                        "Robotic RP"
                        if "ROBOTIC" in x
                        else "Open RP"
                        if "RADICAL" in x
                        else np.nan
                    )
                )
            )
        print("Derived: SURGICAL APPROACH (Open RP vs Robotic RP)")

    # Time to metastasis (months from diagnosis)
    if "DATE OF METASTASIS DX" in df.columns and "DIAGNOSIS DATE" in df.columns:
        for frame in [df, met, nmet]:
            t = date_diff_months(frame, "DATE OF METASTASIS DX")
            if t is not None:
                frame["TIME TO METASTASIS"] = t
        print("Derived: TIME TO METASTASIS")

    # Warn about missing expected columns
    all_expected = set(TABLE1_VARS + SUPPLEMENTARY_EXTRA_VARS + LOVE_VARS)
    missing = [
        v
        for v in all_expected
        if v not in df.columns
        and v not in SKIP_VARS
        and v
        not in {
            "PCT POSITIVE CORES",
            "TIME TO METASTASIS",
            "CLINICAL STAGE BINARY",
            "PATH STAGE BINARY",
            "SURGICAL APPROACH",
            "ANCESTRY GROUPED",
            "LAST FOLLOW UP",
            "ALIVE AT FOLLOW UP",
        }
    ]
    if missing:
        print(f"\nNOTE — {len(missing)} expected column(s) not found in TSV:")
        for c in sorted(missing):
            print(f"  · {c}")

    # ── Tables ────────────────────────────────────────────────────────────────
    t1_vars = [v for v in TABLE1_VARS if v in df.columns or v == "PCT POSITIVE CORES"]
    supp_vars = list(dict.fromkeys(TABLE1_VARS + SUPPLEMENTARY_EXTRA_VARS))

    t1 = build_table(t1_vars, df, met, nmet, include_metonly=False)
    supp = build_table(supp_vars, df, met, nmet, include_metonly=True)
    write_excel([("Table1", t1), ("Supplementary", supp)], stem + "_output.xlsx")

    # ── Figures ───────────────────────────────────────────────────────────────
    plot_table1_pdf(t1, stem + "_table1.pdf")
    plot_swimmer(met, nmet, stem + "_swimmer.pdf")
    plot_love(met, nmet, df, stem + "_love.pdf")

    # ── Legend ────────────────────────────────────────────────────────────────
    print()
    print("── Table footnotes ──────────────────────────────────────────────────")
    print("  †    Matching criterion; p-values shown for completeness.")
    print("  [n=X] follow median (IQR) or in header row: X of 26 patients had data.")
    print("  Fisher's exact: binary/2-category; Chi-square: >2 categories.")
    print("  Chi-square *: expected cell count <5 in ≥1 cell — interpret cautiously.")
    print("── Love plot notes ──────────────────────────────────────────────────")
    print("  Shaded band: |SMD| < 0.10 (conventional balance threshold).")
    print("  Red dots: matching criteria (expected near 0).")
    print("  Blue dots: unmatched covariates.")
    print("── Swimmer plot notes ───────────────────────────────────────────────")
    print("  Metastasis markers read from the TIME TO METASTASIS column (months).")
    print("  Bar length = LAST FOLLOW UP (assumed to be in months).")
    print("────────────────────────────────────────────────────────────────────")


if __name__ == "__main__":
    main()
