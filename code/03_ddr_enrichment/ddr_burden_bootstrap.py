"""
Within-arm competitive null test for DDR gnsRV burden in the metastatic arm.

Question
--------
Is the total gnsRV burden of the 25 prespecified DDR genes (summed across
the 26 metastatic patients) greater than what 25 random non-DDR genes
would produce in the same 26 patients?

Method
------
- Observed: DDR burden vector across 26 metastatic patients is provided
  inline (DDR_OBSERVED below). Total burden = sum of that vector = 80.
- Null:    for each iteration, draw 25 random genes from the gene pool
           (excluding the DDR panel), sum their burden across the same
           26 met patients, store the total.
- p-value: (#iterations where random_total >= DDR_total + 1) / (N + 1).

Caveat
------
This is a *within-arm* competitive null. It tests "DDR > random in met
patients", not "DDR enrichment is metastatic-specific". A positive result
here does not by itself support the metastatic-specific claim — DDR
genes are larger than typical genes (more mutational target) and may
also be enriched over random sets in the non-metastatic arm. To address
the metastatic-specific question, run this script on the non-met arm
too and compare; or use ddr_bootstrap.py which contrasts the two arms
directly.

Usage
-----
    python code/03_ddr_enrichment/ddr_burden_bootstrap.py
    python code/03_ddr_enrichment/ddr_burden_bootstrap.py --n-iterations 100000
    python code/03_ddr_enrichment/ddr_burden_bootstrap.py \\
        --matrix /path/to/gnsrv_per_gene_per_patient.tsv \\
        --n-iterations 50000 --seed 7
"""

import argparse
from pathlib import Path

import numpy as np
import pandas as pd

REPO = Path(__file__).resolve().parents[2]

# Observed DDR gnsRV burden across the 26 metastatic patients
# (one value per patient, in column order). Sum = 80.
DDR_OBSERVED = (
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
DDR_TOTAL_OBSERVED = sum(DDR_OBSERVED)  # 81

# Observed DDR gnsRV burden across the 26 NON-METASTATIC patients (sum = 21).
DDR_OBSERVED_NONMET = (
    1, 1, 0, 0, 0, 0, 0, 1, 2, 0, 1, 0, 3, 1,
    2, 0, 0, 2, 1, 2, 0, 0, 1, 1, 2, 0,
)
DDR_TOTAL_OBSERVED_NONMET = sum(DDR_OBSERVED_NONMET)  # 21

PANEL_SIZE = 25  # number of genes per random panel (matches DDR panel size)
N_MET = 26       # first N_MET cols of the matrix  = metastatic arm
N_NONMET = 26    # next  N_NONMET cols              = non-metastatic arm


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description=__doc__.split("\n\n")[1].strip())
    p.add_argument(
        "--matrix",
        type=Path,
        default=REPO / "data" / "raw" / "gnsrv_per_gene_per_patient.tsv",
        help="Wide gene × patient TSV (col 1 = gene name → row index). "
        f"First {N_MET} patient columns = met arm; next {N_NONMET} = non-met.",
    )
    p.add_argument(
        "--arm", choices=("met", "nonmet", "both"), default="met",
        help="Which arm to test: met (default), nonmet, or both.",
    )
    p.add_argument(
        "--ddr-panel",
        type=Path,
        default=REPO / "data" / "derived" / "ddr_panel_25genes.tsv",
        help="TSV with a 'gene' column listing the 25 DDR genes. "
        "These are excluded from the random pool.",
    )
    p.add_argument("--n-iterations", type=int, default=10_000)
    p.add_argument("--seed", type=int, default=42)
    p.add_argument(
        "--include-ddr-in-pool",
        action="store_true",
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
    p.add_argument(
        "--top-n", type=int, default=20,
        help="Show the top-N right-tail panels with their gene contents "
             "(diagnostic for FLAGS-style hyper-variable genes inflating the null).",
    )
    p.add_argument(
        "--save-panels", type=Path, default=None,
        help="Optional: write every iteration's sampled gene panel + total to this TSV.",
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


def run_arm(
    arm_label: str,
    observed_total: int,
    n_observed: int,
    arm_matrix: pd.DataFrame,
    pool: list[str],
    args: argparse.Namespace,
    save_panels_path: Path | None,
) -> None:
    """Run the bootstrap and report for one arm. Called once per --arm value."""
    rng = np.random.default_rng(args.seed)

    print(f"\n{'=' * 72}")
    print(f"=== {arm_label.upper()} ARM ===")
    print(f"{'=' * 72}")
    print(f"Patients in this arm:    {arm_matrix.shape[1]}")
    print(f"Observed DDR burden:     sum = {observed_total} across {n_observed} patients")
    print(f"Iterations:              {args.n_iterations:,}\n")

    null_totals = np.empty(args.n_iterations, dtype=np.int64)
    sampled_panels = np.empty((args.n_iterations, PANEL_SIZE), dtype=object)
    for i in range(args.n_iterations):
        sampled = rng.choice(pool, size=PANEL_SIZE, replace=False)
        null_totals[i] = int(arm_matrix.loc[sampled].to_numpy().sum())
        sampled_panels[i] = sampled

    n_ge = int((null_totals >= observed_total).sum())
    p_value = (n_ge + 1) / (args.n_iterations + 1)

    print(f"Random 25-gene burden distribution (across {arm_matrix.shape[1]} {arm_label} patients):")
    print(f"  median  = {int(np.median(null_totals))}")
    print(f"  mean    = {null_totals.mean():.1f}")
    print(f"  range   = {int(null_totals.min())} – {int(null_totals.max())}")
    print(f"  95 % CI = ({int(np.percentile(null_totals, 2.5))}, "
          f"{int(np.percentile(null_totals, 97.5))})")
    print()
    print(f"Iterations with random_total ≥ {observed_total}: "
          f"{n_ge:,}/{args.n_iterations:,}")
    print(f"Empirical p (one-sided, random ≥ DDR) = {p_value:.4e}")

    # ── Diagnostic: top-N right-tail panels and the genes driving them ───
    print(f"\n--- Top {args.top_n} right-tail panels ({arm_label} arm) ---")
    print("If the same handful of genes appear repeatedly here, the right tail")
    print("is driven by hyper-variable genes (FLAGS, e.g. TTN, MUC4, MUC16,")
    print("OR-family, HLA, PCLO) — or a non-gene metadata row mislabeled as a gene.")
    top_idx = np.argsort(null_totals)[-args.top_n:][::-1]
    print(f"\n  {'rank':>4}  {'iter':>6}  {'total':>6}   top-5 contributing genes (gene:count)")
    contributing_counter: dict[str, int] = {}
    for rank, i in enumerate(top_idx, start=1):
        panel_genes = sampled_panels[i]
        per_gene = arm_matrix.loc[panel_genes].sum(axis=1).sort_values(ascending=False)
        top5 = per_gene.head(5)
        contribs = "  ".join(f"{g}:{int(v)}" for g, v in top5.items())
        for g, v in per_gene.items():
            if int(v) > 0:
                contributing_counter[g] = contributing_counter.get(g, 0) + int(v)
        print(f"  {rank:>4}  {int(i):>6}  {int(null_totals[i]):>6}   {contribs}")

    print(f"\nGenes contributing most across the top {args.top_n} panels ({arm_label}):")
    print(f"  {'gene':<14} {'cumulative_burden':>18} {'panels_seen':>13}")
    panels_seen: dict[str, int] = {}
    for i in top_idx:
        for g in sampled_panels[i]:
            panels_seen[g] = panels_seen.get(g, 0) + 1
    sorted_contribs = sorted(contributing_counter.items(),
                              key=lambda kv: kv[1], reverse=True)
    for g, total in sorted_contribs[:15]:
        print(f"  {g:<14} {total:>18d} {panels_seen.get(g, 0):>13d}")

    if save_panels_path is not None:
        save_panels_path.parent.mkdir(parents=True, exist_ok=True)
        rows = [{
            "iter": i,
            "total": int(null_totals[i]),
            "genes": ",".join(sampled_panels[i]),
        } for i in range(args.n_iterations)]
        pd.DataFrame(rows).to_csv(save_panels_path, sep="\t", index=False)
        print(f"\nWrote per-iteration panel log ({arm_label}): {save_panels_path}")


def main() -> None:
    args = parse_args()

    # ── Load matrix once, then build the shared random pool ──────────────
    matrix = pd.read_csv(args.matrix, sep="\t", index_col=0)
    matrix = matrix.groupby(level=0).sum()
    needed_cols = N_MET + (N_NONMET if args.arm in ("nonmet", "both") else 0)
    if matrix.shape[1] < needed_cols:
        raise SystemExit(
            f"Matrix has {matrix.shape[1]} patient columns; "
            f"--arm {args.arm} needs at least {needed_cols}."
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

    # ── Per-arm dispatch ─────────────────────────────────────────────────
    arms = ("met", "nonmet") if args.arm == "both" else (args.arm,)
    for arm in arms:
        if arm == "met":
            arm_matrix = matrix.iloc[:, :N_MET]
            obs_total = DDR_TOTAL_OBSERVED
            n_obs = len(DDR_OBSERVED)
        else:  # nonmet
            arm_matrix = matrix.iloc[:, N_MET : N_MET + N_NONMET]
            obs_total = DDR_TOTAL_OBSERVED_NONMET
            n_obs = len(DDR_OBSERVED_NONMET)

        # If two arms requested and --save-panels given, suffix the path
        save_path = args.save_panels
        if save_path is not None and args.arm == "both":
            save_path = save_path.with_name(f"{save_path.stem}_{arm}{save_path.suffix}")

        run_arm(arm, obs_total, n_obs, arm_matrix, pool, args, save_path)


if __name__ == "__main__":
    main()
