"""
Firth penalized logistic regression of DDR gnsRV predictors of metastatic
phenotype (Methods §2.3.2; Supplementary Methods §4.7). Reported as a
multivariable sensitivity analysis.

Fits four model combinations matching the supplement's table:

    Ancestry-adjusted (age + ancestry)
      DDR carrier (binary):  OR = 26.80, 95% CI 3.07-3528.42, p = 0.0006
      DDR burden  (count):   OR =  4.44, 95% CI 2.13-13.05,   p < 0.001

    Fully adjusted (+ biopsy Gleason + surgical margins)
      DDR carrier (binary):  OR = 24.93, 95% CI 2.91-3261.56, p = 0.0008
      DDR burden  (count):   OR =  4.26, 95% CI 2.06-12.31,   p < 0.001

Firth (1993) penalized likelihood handles the quasi-complete separation
that arises when 26/26 metastatic patients carry >=1 DDR gnsRV.

The count predictor (DDR burden as a continuous variable) is fit in
parallel to the binary carrier predictor and reported as a sensitivity
analysis: the count model captures dose-dependence of the DDR-metastasis
association beyond the carrier/non-carrier dichotomy.

Usage
-----
    # All four combinations (default)
    python code/04_firth_regression/ddr_logistic_regression.py

    # Restrict to specific predictor or adjustment level
    python code/04_firth_regression/ddr_logistic_regression.py --predictor binary
    python code/04_firth_regression/ddr_logistic_regression.py --predictor count
    python code/04_firth_regression/ddr_logistic_regression.py --adjustment primary
    python code/04_firth_regression/ddr_logistic_regression.py --adjustment extended
"""

import argparse
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats
from scipy.optimize import brentq, minimize

REPO = Path(__file__).resolve().parents[2]


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    p.add_argument(
        "--data", type=Path,
        default=REPO / "data" / "processed" / "epc_cohort.tsv",
    )
    p.add_argument(
        "--predictor", choices=("binary", "count", "both"), default="both",
        help="DDR predictor type: binary (carrier ≥1 vs 0), count (DDR burden), "
             "or both (default).",
    )
    p.add_argument(
        "--adjustment", choices=("primary", "extended", "both"), default="both",
        help="Adjustment level: primary (age+ancestry), extended "
             "(+biopsy Gleason+surgical margins), or both (default).",
    )
    p.add_argument(
        "--out", type=Path,
        default=REPO / "data" / "processed" / "firth_results.csv",
    )
    return p.parse_args()


def firth_negloglik(beta, X, y):
    mu = 1 / (1 + np.exp(-X @ beta))
    mu = np.clip(mu, 1e-15, 1 - 1e-15)
    W = mu * (1 - mu)
    XtWX = X.T @ (X * W[:, None])
    sign, logdet = np.linalg.slogdet(XtWX)
    penalty = 0.5 * logdet if sign > 0 else 0.0
    loglik = np.sum(y * np.log(mu) + (1 - y) * np.log(1 - mu))
    return -(loglik + penalty)


def firth_fit(X, y):
    res = minimize(firth_negloglik, np.zeros(X.shape[1]), args=(X, y),
                   method="BFGS", options={"maxiter": 500, "gtol": 1e-8})
    return res.x, -res.fun


def firth_lrt_pval(X, y, ll_full, j):
    idx = [i for i in range(X.shape[1]) if i != j]
    _, ll_rest = firth_fit(X[:, idx], y)
    lrt = 2 * (ll_full - ll_rest)
    return float(stats.chi2.sf(max(lrt, 0), df=1))


def firth_profile_ci(X, y, beta, ll_full, j, alpha=0.05):
    chi2_crit = stats.chi2.ppf(1 - alpha, 1) / 2

    def profile_at(b):
        idx_rest = [i for i in range(X.shape[1]) if i != j]

        def neg(beta_rest):
            full = np.insert(beta_rest, j, b)
            return firth_negloglik(full, X, y)

        r = minimize(neg, beta[idx_rest], method="BFGS",
                     options={"maxiter": 200, "gtol": 1e-6})
        return -r.fun

    def boundary(b):
        return ll_full - profile_at(b) - chi2_crit

    mu = 1 / (1 + np.exp(-X @ beta))
    W = mu * (1 - mu)
    try:
        se = float(np.sqrt(np.linalg.inv(X.T @ (X * W[:, None]))[j, j]))
    except np.linalg.LinAlgError:
        se = 1.0
    try:
        lo = brentq(boundary, beta[j] - 10 * se, beta[j], xtol=1e-4)
        hi = brentq(boundary, beta[j], beta[j] + 10 * se, xtol=1e-4)
    except ValueError:
        lo, hi = beta[j] - 1.96 * se, beta[j] + 1.96 * se
    return lo, hi


def build_design(df: pd.DataFrame, predictor_type: str, extended: bool):
    """Build design matrix.

    predictor_type:
        - 'binary' uses ddr_carrier = (ddr_burden >= 1) — primary manuscript model.
        - 'count'  uses ddr_burden directly as a continuous predictor.
    """
    df = df.copy()
    if predictor_type == "binary":
        df["ddr_predictor"] = (df["ddr_burden"] >= 1).astype(int)
        pred_name = "ddr_carrier"
    elif predictor_type == "count":
        df["ddr_predictor"] = df["ddr_burden"].astype(float)
        pred_name = "ddr_burden_count"
    else:
        raise ValueError(f"Unknown predictor_type: {predictor_type}")

    df["ancestry"] = pd.Categorical(
        df["ancestry"],
        categories=["European"] + sorted(set(df["ancestry"]) - {"European"}),
    )
    ancestry_dummies = pd.get_dummies(df["ancestry"], drop_first=True)

    cols = ["ddr_predictor", "age_at_diagnosis"]
    extra = pd.DataFrame(index=df.index)
    if extended:
        cols += ["biopsy_gleason"]
        extra["surgical_margins_binary"] = df["surgical_margins_binary"]

    X_df = pd.concat([df[cols], ancestry_dummies, extra], axis=1)
    X_df = X_df.rename(columns={"ddr_predictor": pred_name})
    complete = X_df.notna().all(axis=1)
    X_df = X_df.loc[complete]
    y = df.loc[complete, "metastatic"].values
    feature_names = ["intercept"] + X_df.columns.tolist()
    X = np.column_stack([np.ones(len(y)), X_df.values]).astype(np.float64)
    return X, y, feature_names, pred_name


def fit_one(df: pd.DataFrame, predictor_type: str, extended: bool) -> dict:
    X, y, names, pred_name = build_design(df, predictor_type, extended)

    beta, ll_full = firth_fit(X, y)
    j = names.index(pred_name)
    lo, hi = firth_profile_ci(X, y, beta, ll_full, j)
    p = firth_lrt_pval(X, y, ll_full, j)

    OR, ci_lo, ci_hi = float(np.exp(beta[j])), float(np.exp(lo)), float(np.exp(hi))
    epv = int(y.sum()) / (len(names) - 1)

    return {
        "adjustment": "extended" if extended else "ancestry-adjusted",
        "predictor_type": predictor_type,
        "predictor_name": pred_name,
        "n": len(y),
        "n_metastatic": int(y.sum()),
        "n_predictors": len(names) - 1,
        "EPV": round(epv, 2),
        "OR": OR,
        "CI_lo": ci_lo,
        "CI_hi": ci_hi,
        "p_LRT": p,
    }


def format_p(p: float) -> str:
    """Match the supplement's table conventions: explicit value if >= 0.0001,
    '<0.001' otherwise."""
    if p < 0.001:
        return "<0.001"
    return f"{p:.4f}"


def main() -> None:
    args = parse_args()
    df = pd.read_csv(args.data, sep="\t")

    adjustments = (
        ["primary", "extended"]
        if args.adjustment == "both"
        else [args.adjustment]
    )
    predictors = (
        ["binary", "count"]
        if args.predictor == "both"
        else [args.predictor]
    )

    results = []
    for adj in adjustments:
        ext = adj == "extended"
        for pred in predictors:
            adj_label = "extended (age+ancestry+Gleason+margins)" if ext else "ancestry-adjusted (age+ancestry)"
            print(f"\n=== Firth penalized logistic regression ===")
            print(f"  Adjustment: {adj_label}")
            print(f"  Predictor:  {pred} ({'ddr_carrier ≥1' if pred == 'binary' else 'ddr_burden count'})")
            r = fit_one(df, pred, ext)
            print(f"  n           = {r['n']} ({r['n_metastatic']} metastatic)")
            print(f"  predictors  = {r['n_predictors']}  (EPV ~ {r['EPV']:.1f})")
            print(f"  OR          = {r['OR']:.2f}")
            print(f"  95% prof CI = {r['CI_lo']:.2f} – {r['CI_hi']:.2f}")
            print(f"  LRT p       = {format_p(r['p_LRT'])}")
            results.append(r)

    # ── Summary table (matches Supplementary Methods §4.7 / table block) ────
    print(f"\n{'=' * 78}")
    print("SUMMARY — matches the table reported in Supplementary Methods §4.7")
    print(f"{'=' * 78}")
    print(f"  {'Adjustment':<20} {'Predictor':<14} {'OR':>7}  {'95% CI':>18}  {'p':>10}")
    print(f"  {'-' * 76}")
    for r in results:
        ci = f"{r['CI_lo']:.2f}–{r['CI_hi']:.2f}"
        pred_lbl = "DDR carrier (bin)" if r["predictor_type"] == "binary" else "DDR burden (cnt)"
        print(f"  {r['adjustment']:<20} {pred_lbl:<14} {r['OR']:>7.2f}  {ci:>18}  {format_p(r['p_LRT']):>10}")
    print(f"{'=' * 78}")

    # ── Save CSV ───────────────────────────────────────────────────────────
    args.out.parent.mkdir(parents=True, exist_ok=True)
    pd.DataFrame(results).to_csv(args.out, index=False)
    print(f"\nWrote {args.out}")


if __name__ == "__main__":
    main()
