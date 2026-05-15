"""
Direct between-arm comparison of DDR gnsRV burden.

Compares per-patient DDR burden in the metastatic arm versus the
non-metastatic arm. Reports descriptive statistics for each arm, a
Mann-Whitney U test (non-parametric two-sample), and a label-permutation
test that makes no distributional assumptions.

Inputs default to the burden vectors below. Override with --met / --nonmet
to pass your own (comma-separated values, or a path to a one-value-per-line
text file).

Usage
-----
    python code/03_ddr_enrichment/ddr_arm_compare.py
    python code/03_ddr_enrichment/ddr_arm_compare.py --n-permutations 100000
    python code/03_ddr_enrichment/ddr_arm_compare.py \\
        --met    "2,4,3,5,3,1,5,1,7,1,2,1,3,1,5,1,4,4,3,7,5,3,6,2,1,1" \\
        --nonmet "1,1,0,0,0,0,0,1,2,0,1,0,3,1,2,0,0,2,1,2,0,0,1,1,2,0"
"""

import argparse
from pathlib import Path

import numpy as np
from scipy import stats

# Default burden vectors — DDR (25-gene panel) per-patient counts
DEFAULT_MET = (
    2,
    4,
    3,
    5,
    3,
    1,
    5,
    1,
    7,
    1,
    2,
    1,
    3,
    1,
    5,
    1,
    4,
    4,
    3,
    7,
    5,
    3,
    6,
    2,
    1,
    1,
)
DEFAULT_NONMET = (
    1,
    1,
    0,
    0,
    0,
    0,
    0,
    1,
    2,
    0,
    1,
    0,
    3,
    1,
    2,
    0,
    0,
    2,
    1,
    2,
    0,
    0,
    1,
    1,
    2,
    0,
)


def parse_vector(arg: str | None, default: tuple[int, ...]) -> np.ndarray:
    if arg is None:
        return np.asarray(default, dtype=int)
    p = Path(arg)
    if p.exists():
        # one value per line, blank lines ignored
        values = [
            int(line.strip()) for line in p.read_text().splitlines() if line.strip()
        ]
    else:
        values = [int(v.strip()) for v in arg.split(",") if v.strip()]
    return np.asarray(values, dtype=int)


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    p.add_argument(
        "--met",
        default=None,
        help="Met-arm burden values: comma-separated or path to a text file.",
    )
    p.add_argument(
        "--nonmet",
        default=None,
        help="Non-met-arm burden values: comma-separated or path to a text file.",
    )
    p.add_argument(
        "--n-permutations",
        type=int,
        default=1_000_000,
        help="Iterations for the label-permutation test (default 1,000,000 to reach p < 10^-6 precision).",
    )
    p.add_argument("--seed", type=int, default=42)
    return p.parse_args()


def describe(label: str, x: np.ndarray) -> None:
    carriers = int((x > 0).sum())
    print(
        f"  {label:>10}:  n={len(x)}  sum={int(x.sum())}  "
        f"mean={x.mean():.2f}  median={int(np.median(x))}  "
        f"IQR=({int(np.percentile(x, 25))}, {int(np.percentile(x, 75))})  "
        f"range={int(x.min())}-{int(x.max())}  "
        f"carriers={carriers}/{len(x)} ({100 * carriers / len(x):.1f}%)"
    )


def permutation_pvalues(
    met: np.ndarray,
    nm: np.ndarray,
    n_permutations: int,
    seed: int,
) -> tuple[float, float, np.ndarray, float]:
    """Label-shuffle the 52 combined values, recompute mean difference,
    return one-sided (greater) and two-sided empirical p-values from the
    same null distribution."""
    rng = np.random.default_rng(seed)
    combined = np.concatenate([met, nm])
    n_met = len(met)
    obs = met.mean() - nm.mean()

    null = np.empty(n_permutations)
    for i in range(n_permutations):
        rng.shuffle(combined)
        null[i] = combined[:n_met].mean() - combined[n_met:].mean()

    p_greater = (int((null >= obs).sum()) + 1) / (n_permutations + 1)
    p_two = (int((np.abs(null) >= abs(obs)).sum()) + 1) / (n_permutations + 1)
    return p_greater, p_two, null, obs


def main() -> None:
    args = parse_args()
    met = parse_vector(args.met, DEFAULT_MET)
    nm = parse_vector(args.nonmet, DEFAULT_NONMET)

    print("DDR per-patient burden — between-arm comparison\n")
    print("Per-arm summary:")
    describe("met", met)
    describe("nonmet", nm)

    obs_diff = float(met.mean() - nm.mean())
    fold = (met.sum() / nm.sum()) if nm.sum() > 0 else float("inf")
    print(f"\n  observed mean difference  (met - nm) = {obs_diff:+.3f}")
    print(
        f"  observed sum ratio        (met / nm) = {fold:.2f}×  "
        f"({int(met.sum())} vs {int(nm.sum())})"
    )

    # ── Mann-Whitney U (non-parametric two-sample) ──────────────────────────
    u_stat, p_mwu_g = stats.mannwhitneyu(met, nm, alternative="greater")
    _, p_mwu_2 = stats.mannwhitneyu(met, nm, alternative="two-sided")
    print(f"\nMann-Whitney U  (U = {u_stat:.1f}):")
    print(f"  greater    (met > nm):  p = {p_mwu_g:.4e}")
    print(f"  two-sided  (met ≠ nm):  p = {p_mwu_2:.4e}")

    # ── Permutation test (label shuffle, both alternatives from same null) ──
    p_perm_g, p_perm_2, null_diffs, obs_diff_perm = permutation_pvalues(
        met,
        nm,
        args.n_permutations,
        args.seed,
    )
    print(f"\nLabel-permutation test  ({args.n_permutations:,} iterations):")
    print(
        f"  null mean diff distribution: "
        f"mean={null_diffs.mean():+.3f}, std={null_diffs.std():.3f}, "
        f"range=({null_diffs.min():+.2f}, {null_diffs.max():+.2f})"
    )
    print(f"  observed diff = {obs_diff_perm:+.3f}")
    print(f"  greater    (met > nm):  p = {p_perm_g:.4e}")
    print(f"  two-sided  (|met−nm|):  p = {p_perm_2:.4e}")

    # ── Welch t-test (parametric, for completeness — counts violate normality) ──
    t_stat, p_t_g = stats.ttest_ind(met, nm, equal_var=False, alternative="greater")
    _, p_t_2 = stats.ttest_ind(met, nm, equal_var=False, alternative="two-sided")
    print(f"\nWelch t-test  (t = {t_stat:.3f}, parametric — interpret cautiously):")
    print(f"  greater    (met > nm):  p = {p_t_g:.4e}")
    print(f"  two-sided  (met ≠ nm):  p = {p_t_2:.4e}")


if __name__ == "__main__":
    main()
