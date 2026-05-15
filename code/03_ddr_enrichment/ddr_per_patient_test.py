"""
Per-patient empirical p-values for DDR observed burden against an in-script
bootstrap null distribution, with Fisher's combined test.

For each of the 26 patients in each arm, ask: "is this patient's observed
DDR burden higher than what 25 random non-DDR genes produce in this same
patient?" The null for patient i is built by sampling 25 random genes from
the protein-coding gene pool (excluding DDR and optionally FLAGS), summing
their gnsRV counts for that patient, and repeating N times. Patient-level
empirical p-values are then combined into a single arm-level p using
Fisher's method.

The bootstrap null is computed on-the-fly from the gene × patient matrix
(the same input used by ddr_burden_bootstrap.py); with a fixed seed
(default 42) the null distribution is fully reproducible.

This complements the sum-based total-burden test: instead of asking "is
the DDR total larger than typical random totals?", it asks "is the
per-patient pattern of DDR burden uniformly enriched, or driven by
outliers?"

Inputs
------
  --matrix         gnsrv_per_gene_per_patient.tsv     wide gene × patient matrix
  --ddr-panel      ddr_panel_25genes.tsv              DDR genes (excluded from random pool)
  --exclude-genes  flags_shyr2014.tsv                 additional exclusions (FLAGS)
  Observed vectors are hardcoded as DEFAULT_MET / DEFAULT_NONMET.

Usage
-----
    python code/03_ddr_enrichment/ddr_per_patient_test.py
    python code/03_ddr_enrichment/ddr_per_patient_test.py --n-iterations 50000
    python code/03_ddr_enrichment/ddr_per_patient_test.py --no-exclude
    python code/03_ddr_enrichment/ddr_per_patient_test.py \\
        --save-null per_patient_null.tsv
"""

import argparse
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

REPO = Path(__file__).resolve().parents[2]
HERE = Path(__file__).resolve().parent

DEFAULT_MET = (
    2, 4, 3, 5, 3, 1, 5, 1, 7, 1, 2, 1, 3, 1,
    5, 1, 4, 4, 3, 7, 5, 3, 6, 2, 1, 1,
)
DEFAULT_NONMET = (
    1, 1, 0, 0, 0, 0, 0, 1, 2, 0, 1, 0, 3, 1,
    2, 0, 0, 2, 1, 2, 0, 0, 1, 1, 2, 0,
)

PANEL_SIZE = 25  # number of genes per random panel (matches DDR panel size)
N_MET = 26       # first N_MET cols of the matrix  = metastatic arm
N_NONMET = 26    # next  N_NONMET cols              = non-metastatic arm


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    p.add_argument(
        "--matrix", type=Path,
        default=REPO / "data" / "raw" / "gnsrv_per_gene_per_patient.tsv",
        help=f"Wide gene × patient TSV (col 1 = gene name → row index). "
             f"First {N_MET} patient columns = met arm; next {N_NONMET} = non-met.",
    )
    p.add_argument(
        "--ddr-panel", type=Path,
        default=REPO / "data" / "derived" / "ddr_panel_25genes.tsv",
        help="TSV with a 'gene' column listing the 25 DDR genes. "
             "These are excluded from the random pool.",
    )
    p.add_argument(
        "--include-ddr-in-pool", action="store_true",
        help="Allow DDR genes to be drawn into random panels (default: exclude them).",
    )
    p.add_argument(
        "--exclude-genes", type=Path,
        default=REPO / "data" / "derived" / "flags_shyr2014.tsv",
        help="TSV (header `gene`) of additional genes to exclude from the "
             "random pool — typically a FLAGS list (Shyr 2014). "
             "Lines starting with `#` are ignored. Use --no-exclude to disable.",
    )
    p.add_argument(
        "--no-exclude", action="store_true",
        help="Disable the FLAGS-based exclusion (use the full gene pool minus DDR).",
    )
    p.add_argument("--n-iterations", type=int, default=10_000)
    p.add_argument("--seed", type=int, default=42)
    p.add_argument(
        "--save-null", type=Path, default=None,
        help="Optional: write the (n_iterations × n_patients) null matrix to this TSV.",
    )
    return p.parse_args()


def load_exclude_list(path: Path) -> set[str]:
    """Load a one-column TSV of gene symbols to exclude. Tolerates `#` comments
    and either a header line `gene` or no header."""
    if not path.exists():
        print(f"WARNING: exclusion list {path} not found — proceeding with no exclusions.")
        return set()
    genes: set[str] = set()
    for line in path.read_text().splitlines():
        s = line.strip()
        if not s or s.startswith("#") or s.lower() == "gene":
            continue
        genes.add(s.split("\t")[0])
    return genes


def build_null_matrix(
    matrix: pd.DataFrame,
    pool: list[str],
    n_iterations: int,
    seed: int,
) -> np.ndarray:
    """Sample n_iterations random PANEL_SIZE-gene panels and compute per-patient
    burdens. Returns an (n_iterations, n_patients) integer array.

    This mirrors the bootstrap loop in ddr_burden_bootstrap.py: the same seed
    produces the same gene draws, so per-patient null distributions here are
    numerically identical to what that script would produce per-arm."""
    rng = np.random.default_rng(seed)
    n_patients = matrix.shape[1]
    null_matrix = np.empty((n_iterations, n_patients), dtype=np.int64)
    for i in range(n_iterations):
        sampled = rng.choice(pool, size=PANEL_SIZE, replace=False)
        null_matrix[i] = matrix.loc[sampled].to_numpy().sum(axis=0)
    return null_matrix


def per_patient_pvalues(observed: np.ndarray, null_matrix: np.ndarray) -> np.ndarray:
    """Empirical one-sided p_i = (#null_rows where null_i >= observed_i + 1) / (N+1)."""
    n_iter = null_matrix.shape[0]
    n_ge = (null_matrix >= observed[None, :]).sum(axis=0).astype(int)
    return (n_ge + 1) / (n_iter + 1)


def fisher_combined(p_values: np.ndarray) -> tuple[float, float, int]:
    """Fisher's combined probability test.
    Returns (chi2_stat, combined_p, df). df = 2 * k."""
    p_clipped = np.clip(p_values, 1e-300, 1.0)
    chi2_stat = float(-2.0 * np.sum(np.log(p_clipped)))
    df = 2 * len(p_values)
    combined_p = float(stats.chi2.sf(chi2_stat, df))
    return chi2_stat, combined_p, df


def analyse_arm(
    label: str,
    observed: tuple[int, ...],
    null_matrix: np.ndarray,
    patient_ids: list[str],
) -> dict:
    obs = np.asarray(observed, dtype=int)
    n_iter = null_matrix.shape[0]

    p_per = per_patient_pvalues(obs, null_matrix)
    chi2, p_comb, df_chi = fisher_combined(p_per)

    n_lt_05 = int((p_per < 0.05).sum())
    n_lt_001 = int((p_per < 0.001).sum())
    n_min = int((p_per <= 1.0 / (n_iter + 1) + 1e-12).sum())  # capped at most extreme

    print(f"\n=== {label} arm ===")
    print(f"Null:       {n_iter:,} random {PANEL_SIZE}-gene panels × {null_matrix.shape[1]} patients")
    print(f"Observed:   sum = {int(obs.sum())} across {len(obs)} patients\n")

    print(f"  {'pid':>5} {'obs':>5} {'null_mean':>9} {'null_max':>8} {'rank':>6} {'p_emp':>10}")
    null_mean = null_matrix.mean(axis=0)
    null_max = null_matrix.max(axis=0)
    rank = (null_matrix < obs[None, :]).sum(axis=0)
    for i in range(len(obs)):
        marker = "  *" if p_per[i] < 0.05 else ""
        print(f"  {patient_ids[i]:>5} {int(obs[i]):>5} {null_mean[i]:>9.2f} "
              f"{int(null_max[i]):>8d} {int(rank[i]):>6d}/{n_iter} "
              f"{p_per[i]:>10.4e}{marker}")

    print(f"\n  Distribution of per-patient p-values:")
    print(f"    median = {float(np.median(p_per)):.3e}")
    print(f"    geomean= {float(np.exp(np.log(np.clip(p_per, 1e-300, 1)).mean())):.3e}")
    print(f"    p < 0.05:    {n_lt_05}/{len(obs)} patients")
    print(f"    p < 0.001:   {n_lt_001}/{len(obs)} patients")
    print(f"    at floor (≤ 1/N+1):  {n_min}/{len(obs)} patients")

    smallest = int(np.argmin(p_per))
    largest = int(np.argmax(p_per))
    print(f"    smallest p: patient {patient_ids[smallest]}  "
          f"obs={int(obs[smallest])}  p={p_per[smallest]:.3e}")
    print(f"    largest  p: patient {patient_ids[largest]}  "
          f"obs={int(obs[largest])}  p={p_per[largest]:.3e}")

    print(f"\n  Fisher's combined test (k={len(obs)} patients):")
    print(f"    χ² = {chi2:.2f},  df = {df_chi},  combined p = {p_comb:.4e}")

    # Sanity binomial: under H0, expect ~5 % of patients to have p < 0.05.
    expected_05 = 0.05 * len(obs)
    binom_p = stats.binomtest(n_lt_05, len(obs), p=0.05, alternative="greater").pvalue
    print(f"  Sanity check (binomial): under H0 expect ~{expected_05:.1f} patients with "
          f"p < 0.05; observed {n_lt_05} → binom p = {binom_p:.4e}")

    return {
        "label": label,
        "n_iter": n_iter,
        "obs_sum": int(obs.sum()),
        "p_per": p_per,
        "n_lt_05": n_lt_05,
        "n_lt_001": n_lt_001,
        "n_at_floor": n_min,
        "fisher_chi2": chi2,
        "fisher_p": p_comb,
        "fisher_df": df_chi,
        "binom_p_lt_05": binom_p,
    }


def main() -> None:
    args = parse_args()
    print("DDR per-patient empirical test (with in-script bootstrap)")
    print(f"Bootstrap: {args.n_iterations:,} random {PANEL_SIZE}-gene panels, seed = {args.seed}\n")

    # ── Load matrix and build the random pool ──────────────────────────────
    matrix = pd.read_csv(args.matrix, sep="\t", index_col=0)
    matrix = matrix.groupby(level=0).sum()  # collapse duplicate gene rows
    if matrix.shape[1] < N_MET + N_NONMET:
        raise SystemExit(
            f"Matrix has {matrix.shape[1]} patient columns; "
            f"need at least {N_MET + N_NONMET}."
        )

    ddr_genes = set(pd.read_csv(args.ddr_panel, sep="\t")["gene"].tolist())
    flags = set() if args.no_exclude else load_exclude_list(args.exclude_genes)
    excluded = (set() if args.include_ddr_in_pool else ddr_genes) | flags
    pool = [g for g in matrix.index if g not in excluded]

    if len(pool) < PANEL_SIZE:
        raise SystemExit(
            f"Random pool has only {len(pool)} genes; need at least {PANEL_SIZE}."
        )

    print(f"Matrix:                  {matrix.shape[0]} genes × {matrix.shape[1]} patients "
          f"(cols 1-{N_MET} = met; cols {N_MET+1}-{N_MET+N_NONMET} = non-met)")
    print(f"DDR genes in matrix:     {sum(1 for g in ddr_genes if g in matrix.index)}/{len(ddr_genes)}  "
          f"({'kept' if args.include_ddr_in_pool else 'excluded'} from pool)")
    print(f"FLAGS in matrix:         {sum(1 for g in flags if g in matrix.index)}/{len(flags)}  "
          f"({'excluded' if flags else 'no exclusion list applied'})")
    print(f"Random pool:             {len(pool)} genes")

    # ── Single bootstrap covering both arms (52 patients) ───────────────────
    null_matrix_full = build_null_matrix(matrix, pool, args.n_iterations, args.seed)

    # ── Optional: save null matrix to TSV ───────────────────────────────────
    if args.save_null is not None:
        args.save_null.parent.mkdir(parents=True, exist_ok=True)
        save_df = pd.DataFrame(
            null_matrix_full,
            columns=list(matrix.columns),
        )
        save_df.insert(0, "iter", range(args.n_iterations))
        save_df.to_csv(args.save_null, sep="\t", index=False)
        print(f"\nSaved null matrix to: {args.save_null}")

    # ── Split into per-arm null matrices ────────────────────────────────────
    met_null = null_matrix_full[:, :N_MET]
    nm_null = null_matrix_full[:, N_MET : N_MET + N_NONMET]
    met_patient_ids = [str(c) for c in matrix.columns[:N_MET]]
    nm_patient_ids = [str(c) for c in matrix.columns[N_MET : N_MET + N_NONMET]]

    # ── Run per-arm analysis ────────────────────────────────────────────────
    met = analyse_arm("Metastatic", DEFAULT_MET, met_null, met_patient_ids)
    nm = analyse_arm("Non-metastatic", DEFAULT_NONMET, nm_null, nm_patient_ids)

    # ── Side-by-side summary ────────────────────────────────────────────────
    print("\n=== Side-by-side summary ===")
    print(f"  {'arm':<15} {'sum':>5} {'p<0.05':>9} {'p<0.001':>9} {'Fisher p':>13} {'binom p':>12}")
    for r in (met, nm):
        print(f"  {r['label']:<15} {r['obs_sum']:>5} "
              f"{r['n_lt_05']:>4}/{len(DEFAULT_MET):<2}    "
              f"{r['n_lt_001']:>4}/{len(DEFAULT_MET):<2}     "
              f"{r['fisher_p']:>10.3e}   {r['binom_p_lt_05']:>10.3e}")


if __name__ == "__main__":
    main()
