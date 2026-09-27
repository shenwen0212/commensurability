## run_demo.R -- end-to-end demonstration of the three-layer commensurability
## diagnostic on the synthetic data in demo/. No primary data download required.
##
##   Rscript demo/run_demo.R            (run from the repository root)
##
## The demo shows, in order:
##   Layer A  the published label is poorly reconstructed by marker rules, and the
##            within-label marker rate is not discriminative
##   Layer B  the effect size changes with the exposure definition, and the
##            conclusion can reverse sign
##   Layer C  reporting helpers (conclusion_flips, check_report)

suppressPackageStartupMessages(library(data.table))
source("commensurability.R")

D <- fread("demo/demo_cells.csv")
cat(sprintf("demo data: %d cells, %d samples, %d cell types\n\n",
            nrow(D), uniqueN(D$sample), uniqueN(D$final_CT)))

###############################################################################
cat("=========== Layer A: label diagnostics ===========\n")
ta <- label_diagnostics(D, label = "lab_published",
                        rules = c("lab_published", "rule_naive", "rule_aligned", "rule_strict"),
                        reference_rule = "lab_published")
print(ta[, .(rule, n_rule, n_label,
             jaccard = round(jaccard_vs_label, 3),
             p_in = round(p_in_label, 3), p_out = round(p_out_label, 3),
             R = round(R, 2), marker_in_label = round(frac_rule_in_label, 3))])
cat("\nInterpretation: for a rule to point at the label, R must exceed 1.")
cat("\nA high marker_in_label value alone is not evidence, because p_out may be equally high.\n")

###############################################################################
cat("\n=========== Layer B: effect sweep across the definition matrix ===========\n")
res <- effect_sweep(
  dt = D,
  sample_col = "sample", x_col = "x", y_col = "y", ct_col = "final_CT",
  disease_col = "disease_status",
  exposures = list(published = "lab_published",
                   naive     = "rule_naive",
                   aligned   = "rule_aligned"),
  reference = "ref_pop",
  partner_defs = list(activated = c("FibroblastAct"),
                      all       = c("FibroblastAct", "FibroblastAll")),
  radii = c(10, 20, 30, 50, 100),
  min_cells = 3)

cat(sprintf("per-sample estimates: %d rows | pooled estimates: %d rows\n",
            nrow(res$sweep), nrow(res$meta)))
if (nrow(res$meta)) {
  cat("\npooled odds ratios by exposure definition and partner:\n")
  print(dcast(res$meta, exposure ~ partner + radius, value.var = "OR")[
          , lapply(.SD, function(x) if (is.numeric(x)) round(x, 3) else x)], row.names = FALSE)
}

###############################################################################
cat("\n=========== Does the conclusion depend on the definition? ===========\n")
if (nrow(res$meta)) {
  fl <- conclusion_flips(res$meta)
  print(fl)
  cat(sprintf("\nconfigurations with BOTH significant enrichment and depletion: %d of %d\n",
              sum(fl$flips), nrow(fl)))
}

###############################################################################
cat("\n=========== Reporting audit ===========\n")
audit_tab <- data.table(
  label_in_marker_rate = ta[rule == "rule_naive", frac_rule_in_label],
  marker_out_rate      = ta[rule == "rule_naive", p_out_label],
  R                    = ta[rule == "rule_naive", R],
  mean_jaccard         = mean(ta[!is.na(jaccard_vs_label), jaccard_vs_label]),
  n_exposure_defs      = uniqueN(res$meta$exposure),
  n_partner_defs       = uniqueN(res$meta$partner),
  n_radii              = uniqueN(res$meta$radius),
  alignment_delta      = NA_real_,
  alignment_pass       = FALSE)
cat("\naudit of a report that omits the alignment check and reports a single definition:\n")
print(check_report(audit_tab))

cat("\nDONE\n")
