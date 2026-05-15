"""
Australian EPC replication — one-sided binomial test of carrier frequency
against gnomAD-derived population expectation (Methods §2.4; Results;
Fig. 3A; Supplementary Methods §4.6; Table 2).

REPRODUCED MANUSCRIPT VALUES (default flags)

    p_carrier   : 0.311  ( = 1 - prod_i (1 - v4.1.1_AF_i) )

    | Subgroup                       |  n |  k | expected | p (one-sided) |
    |--------------------------------|----|----|----------|---------------|
    | Metastatic (bone/visceral)     | 17 |  9 |    5.30  |    0.051      |
    | Metastatic + BCR               | 29 | 14 |    9.03  |    0.040      |
    | Metastatic + BCR + nodes       | 33 | 16 |   10.28  |    0.028      |
    | Non-metastatic (control)       | 14 |  5 |    4.36  |    0.454      |

FILTER RULE

    Default       : v4.1.1_wes_af_frozen < 0.02
    --current     : v4.1.1_wes_af_current < 0.02 (post-update gnomAD values)

    "Frozen" reflects the gnomAD release state at the time the variant set
    was finalized for submission. One variant (SCARF1 rs61738531) was
    excluded under the frozen filter (WES >= 2%) but would be retained
    under the current filter (WES = 0.01948 after a subsequent gnomAD
    release).

INPUTS

    data/derived/gnsrv_56_screen.tsv       (56 variants with MAF + annotation)
    data/derived/australian_epc_counts.tsv (carrier counts per subgroup)

OUTPUT

    data/processed/australian_binomial_results.csv
"""

import argparse
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import binomtest

REPO = Path(__file__).resolve().parents[2]
LOOKUP = REPO / "data" / "derived" / "gnsrv_56_screen.tsv"
COUNTS = REPO / "data" / "derived" / "australian_epc_counts.tsv"
OUTPUT = REPO / "data" / "processed" / "australian_binomial_results.csv"

MAF_THRESHOLD = 0.02
FILTER_COL_DEFAULT = "v411_wes_af_frozen"
FILTER_COL_CURRENT = "v411_wes_af_current"
P_CARRIER_COL = "v411_af"


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    p.add_argument(
        "--current-snapshot",
        action="store_true",
        help="Filter on v4.1.1 WES values from the current gnomAD release "
        "instead of the frozen analysis-date snapshot. "
        "p_carrier=0.321 (sensitivity).",
    )
    p.add_argument("--lookup", type=Path, default=LOOKUP)
    p.add_argument("--counts", type=Path, default=COUNTS)
    p.add_argument("--output", type=Path, default=OUTPUT)
    return p.parse_args()


def load_lookup(path: Path) -> pd.DataFrame:
    return pd.read_csv(path, sep="\t", comment="#", na_values=["NA"])


def apply_filter(df: pd.DataFrame, filter_col: str) -> pd.DataFrame:
    mask = df[filter_col].notna() & (df[filter_col] < MAF_THRESHOLD)
    return df.loc[mask].copy()


def compute_p_carrier(retained: pd.DataFrame) -> float:
    f = retained[P_CARRIER_COL].astype(float).clip(0, 1).values
    return float(1 - np.prod(1 - f))


def main() -> None:
    args = parse_args()
    filter_col = FILTER_COL_CURRENT if args.current_snapshot else FILTER_COL_DEFAULT
    snapshot_label = (
        "current gnomAD" if args.current_snapshot else "frozen (analysis-date) gnomAD"
    )

    df = load_lookup(args.lookup)
    retained = apply_filter(df, filter_col)
    p_carrier = compute_p_carrier(retained)

    print(f"Snapshot     : {snapshot_label}")
    print(f"Filter rule  : {filter_col} < {MAF_THRESHOLD}")
    print(f"p_carrier    : {p_carrier:.4f}")
    print()

    counts = pd.read_csv(args.counts, sep="\t", comment="#")
    print(f"{'subgroup':<35} {'n':>3} {'k':>3} {'exp':>6} {'p (1-sided)':>12}")
    print("-" * 65)
    rows = []
    for _, r in counts.iterrows():
        n, k = int(r["n"]), int(r["observed_carriers"])
        expected = n * p_carrier
        p = binomtest(k=k, n=n, p=p_carrier, alternative="greater").pvalue
        print(f"{r['subgroup']:<35} {n:>3} {k:>3} {expected:>6.2f} {p:>12.4f}")
        rows.append(
            {
                **r.to_dict(),
                "snapshot": snapshot_label,
                "p_carrier": p_carrier,
                "expected_recomputed": expected,
                "p_recomputed": p,
            }
        )

    args.output.parent.mkdir(parents=True, exist_ok=True)
    pd.DataFrame(rows).to_csv(args.output, index=False)
    print(f"\nWrote {args.output.relative_to(REPO)}")


if __name__ == "__main__":
    main()
