# commensurability

Diagnostics for **cell-type label commensurability** in spatial omics: a toolkit for
testing whether two cell-type labels can legitimately be compared, and for measuring
how much a spatial conclusion depends on the definitions used to produce it.

Companion code for the manuscript:

> *Cell-type labels are not interchangeable measurements: definition alignment decides
> the sign of spatial neighbourhood effects* (manuscript in preparation)

---

## Why this exists

Cross-platform reproducibility is the usual test of a spatial-omics claim, but that
test assumes the two analyses measured **the same cell population**. The assumption is
almost never checked, and when it fails the resulting "platform difference" may be
nothing of the kind. In the datasets we examined, changing only the cell-type
definition — holding tissue, platform, partner definition and radius fixed — reversed
a significant spatial enrichment into a significant depletion.

Two further findings motivate the toolkit:

* The statistic usually offered as proof that a label is sound — the fraction of cells
  **inside** the label expressing its namesake marker — has no discriminative power.
  In our data the label named after a marker had a *lower* within-label marker rate
  than a different label. The useful quantity is the between-label ratio.
* Overlap with a published label (e.g. Jaccard) does **not** predict whether a
  re-derived definition reproduces another platform; the discriminative marker rate
  does. Marker fidelity is not construct validity.

---

## Installation

No package installation is required. The functions are a single sourced file.

```r
source("commensurability.R")
```

Dependencies: `data.table`, `RANN`, `metafor` (all on CRAN).

---

## Quick start

Run the self-contained demonstrations, which need no data download:

```bash
Rscript demo/make_demo_data.R        # writes demo/demo_cells.csv (synthetic, 16k cells)
Rscript demo/run_demo.R              # the three diagnostic layers
Rscript demo/sensitivity_demo.R      # null calibration and the sensitivity ladder
Rscript tests/run_tests.R            # test suite (18 checks)
```

## Constructing a valid null

`demo/sensitivity_demo.R` demonstrates the single most important practical point, on
synthetic data:

* A label drawn at random from **non-tag** cells carries no tag-specific information.
  Compared against an **annotated** reference population it nonetheless yields odds
  ratios up to 2.5 with **6 of 10 configurations significant** — because the reference
  is not exchangeable with the label.
* When the label and the reference are drawn **the same way**, the same estimator
  returns odds ratios of 0.84–1.14 with **0 of 10 significant**. That is the correct
  null.

The methodological consequence is that "there is no effect" cannot be inferred from a
diagnostic statistic such as a low Jaccard. It has to be established by construction:
if you want to know what your estimator reports when nothing is happening, you must
build a null in which the label and the reference are exchangeable by design.

## Random seeds

Every script sets a seed before drawing anything, and the seeds are listed here so that
results can be regenerated exactly.

| Script | Seed |
|---|---|
| `demo/make_demo_data.R` | 20260101 |
| `demo/sensitivity_demo.R` | 20260101 (data), 11, 12 (null constructions), 1000 + β (ladder) |
| `tests/run_tests.R` | inherits the demo data; no draws beyond it |

The analysis scripts that produced the manuscript fix their seeds in the same way; the
seed used by each is stated at the top of the corresponding script in the companion
analysis repository.

For real data, the three layers are called as follows.

```r
library(data.table)
source("commensurability.R")

# ---- Layer A: can the published label be reconstructed from marker expression? ----
ta <- label_diagnostics(
        dt             = cells,                 # one row per cell
        label          = "published_label",     # logical column
        rules          = c("rule_naive", "rule_aligned", "rule_strict"),
        reference_rule = NULL)                  # a rule that IS the label, if any
# Read R: only R > 1 indicates the rule points at the label.
# Ignore frac_rule_in_label unless p_out_label is also reported.

# ---- Layer B: how much does the conclusion depend on the definitions? ----
res <- effect_sweep(
        dt          = cells,                    # MUST contain all lineages
        sample_col  = "sample", x_col = "x", y_col = "y", ct_col = "final_CT",
        disease_col = "disease_status",
        exposures   = list(published = "published_label",
                           naive     = "rule_naive",
                           aligned   = "rule_aligned"),
        reference   = "reference_pop",
        partner_defs = list(activated = c("Activated Fibrotic FBs","Myofibroblasts"),
                            all       = c("Alveolar FBs","Activated Fibrotic FBs",
                                          "Inflammatory FBs","Myofibroblasts")),
        radii      = c(10, 20, 30, 50, 100))
conclusion_flips(res$meta)      # flips == TRUE means the sign depends on the definition

# ---- Layer C: do two platforms agree once definitions are aligned? ----
platform_compare(res_xenium$meta, res_codex$meta, "Xenium (RNA)", "CODEX (protein)")

# ---- Reporting audit (S1-S5 minimal reporting set) ----
check_report(my_report_table)
```

---

## The three diagnostic layers

| Layer | Function | Question it answers |
|---|---|---|
| **A — label** | `label_diagnostics()` | Can this label be migrated to another platform, and is it still the same population? |
| **B — effect** | `effect_sweep()` | How much does the reported effect depend on the exposure and partner definitions and on the radius? |
| **C — platform** | `platform_compare()` | Once definitions are aligned, do the platforms still disagree? |
| summary | `conclusion_flips()` | In how many configurations does the conclusion change sign? |
| audit | `check_report()` | Does a reported result meet the minimal reporting set? |

### Minimal reporting set (S1–S5)

* **S1** label prevalence and within-label marker-positive rate
* **S2** **outside-label** marker-positive rate and the discriminative ratio R
* **S3** the Jaccard ladder against at least three candidate marker rules
* **S4** the effect across ≥2 exposure definitions × ≥2 partner definitions × ≥3 radii
* **S5** for cross-platform claims: the alignment rule used, |Δθ|, and whether it passed

---

## Three implementation traps this code handles

These are the ways this analysis most often goes silently wrong. The code avoids all
three, and any reimplementation should too.

1. **`nn2()` takes coordinate matrices, never row indices.** Passing indices instead of
   coordinates produces plausible-looking but meaningless distances, with no error.
2. **The per-sample table must contain all lineages.** If the working table is filtered
   to epithelial cells before looking for fibroblasts, the partner type is absent and
   the sweep returns nothing rather than failing.
3. **Pooling must be per sample, not per cell.** At single-cell *n*, pooled p-values are
   not interpretable. This code estimates per sample and pools by random-effects
   meta-analysis (REML, Knapp–Hartung).

A fourth trap is guarded in the analysis scripts rather than here: a mistyped column
name in a definition ladder yields zero rows silently. The companion scripts include an
integrity gate that validates every ladder column against the data and fails loudly.

---

## Licence

MIT. See `LICENSE`.

## Data

The synthetic data in `demo/` are generated by `demo/make_demo_data.R` and contain no
real observations. Primary datasets used in the manuscript are public; see the
manuscript's data availability statement. This repository redistributes none of them.

## Citation

If you use this code, please cite the manuscript above (full citation to follow
publication) and the archived version of this repository, which carries a DOI.
