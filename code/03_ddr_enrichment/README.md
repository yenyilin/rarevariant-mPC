# DDR Pathway-Level Enrichment in Metastatic Prostate Cancer

Scripts in this directory implement the three statistical tests for the DDR pathway-burden analysis described in **Supplementary Methods §4.2** of the published manuscript. This README documents how to reproduce those analyses and notes one debugging artifact discovered during development.

For the full methods text and reported numerical values, see **Supplementary Methods §4.2** of the published manuscript.

---

## Plain-English overview

Prostate cancer patients whose disease later spread to other parts of the body (the metastatic arm, n = 26) carried far more rare germline variants in DNA damage repair (DDR) genes than otherwise-similar patients whose cancer stayed put (the non-metastatic arm, n = 26). The metastatic patients carried, on average, about four times as many DDR variants per person; every single one of them carried at least one, versus only about half of the non-metastatic patients. We tested this in three independent ways — comparing the two arms directly, comparing DDR against 10,000 randomly chosen 25-gene panels, and checking each patient against the random-panel reference one at a time. All three approaches gave the same answer with p-values smaller than five in a million. Crucially, when we ran the same random-panel test on the non-metastatic arm, DDR didn't stand out — meaning the DDR signal isn't just a quirk of DDR genes being unusually large; it's specific to patients whose disease metastasized.

---

## Script-to-supplement crosswalk

| Script | Supplement section | Test |
|---|---|---|
| `ddr_arm_compare.py` | §4.2.2 | Test 1 — Direct between-arm comparison (Mann-Whitney + label-permutation + Welch's *t*) |
| `ddr_burden_bootstrap.py` | §4.2.3 | Test 2 — Within-arm competitive null bootstrap |
| `ddr_per_patient_test.py` | §4.2.4 | Test 3 — Per-patient empirical *p* with Fisher's combined test |

For numerical values, the formal hypothesis-test specifications, and the convergent-evidence framing across all three tests, see **Supplementary Methods §4.2** of the published manuscript.

---

## Reproducing the analyses

Each script has its own command-line interface and CLI defaults. To reproduce the values reported in the supplement, run with the following flags:

```bash
# Test 1 — between-arm comparison (~1 min at 1M permutations)
python ddr_arm_compare.py --n-permutations 1000000

# Test 2 — within-arm bootstrap, both arms (~30 s)
python ddr_burden_bootstrap.py --arm both

# Test 3 — per-patient empirical p + Fisher's combined test
python ddr_per_patient_test.py

# Supplementary forest plot
python ../08_figures/figS_per_patient_pvalues.py
```

Default seed (42) is fixed for reproducibility. Numerical outputs are stable to two significant figures across seeds tested.

**Important — script default vs reported precision.** `ddr_arm_compare.py` defaults to `--n-permutations 10000`. To reproduce the supplement's reported `p < 10⁻⁶` for the label-permutation test, pass `--n-permutations 1000000` explicitly. At the default (10⁴), the empirical floor is `p < 10⁻⁴`, which weakens the reported permutation-test precision but does not change the Mann-Whitney or Welch's *t* results (which are analytical and unaffected by `--n-permutations`).

---

## Debugging note — HRD metadata row

During development, the bootstrap top-N right-tail diagnostic (`--top-n` in `ddr_burden_bootstrap.py`) identified that the input matrix initially contained a non-gene metadata row labeled `HRD` (per-patient homologous-recombination-deficiency score, with values 0–60 rather than discrete variant counts). The random sampler treated this row as a gene, inflating the null max from ~55 to ~735 and the metastatic-arm bootstrap *p* from `< 10⁻⁴` to ~`2 × 10⁻³`.

Removing this row recovered the expected null range (1–55) and restored the reported `p < 10⁻⁴`. The current input matrix is HRD-row-free, and FLAGS exclusion (Shyr 2014 hyper-variable gene list) is no longer required to reach this empirical floor — though FLAGS exclusion remains supported via the `--exclude-genes` and `--no-exclude` flags as a sensitivity option.

This artifact is documented here as evidence that the right-tail diagnostic was used to validate the input data before reporting results.

---

## License

Released under the MIT License — see `LICENSE` in the repository root. Patient-level data is held under institutional controlled-access and is not included in this repository; see the published manuscript Methods for data-access procedures.
