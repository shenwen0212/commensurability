## tests/run_tests.R -- self-contained test suite for commensurability.R
## Uses only the synthetic demo data, so it runs anywhere.
##   Rscript tests/run_tests.R      (run from the repository root)

suppressPackageStartupMessages(library(data.table))
source("commensurability.R")

fails <- 0
ok <- function(cond, label) {
  if (isTRUE(cond)) cat(sprintf("  PASS  %s\n", label))
  else { cat(sprintf("  FAIL  %s\n", label)); fails <<- fails + 1 }
}

D <- fread("demo/demo_cells.csv")
cat(sprintf("test data: %d cells, %d samples\n\n", nrow(D), uniqueN(D$sample)))

cat("--- label_diagnostics ---\n")
ta <- label_diagnostics(D, label = "lab_published",
                        rules = c("lab_published","rule_naive","rule_aligned","rule_strict"),
                        reference_rule = "lab_published")
ok(is.na(ta[rule == "lab_published", jaccard_vs_label]),
   "reference_rule yields NA Jaccard for the label itself")
ok(ta[rule == "rule_naive", R] > 1, "a discriminative rule has R > 1")
ok(all(c("p_in_label","p_out_label","R","lambda") %in% names(ta)),
   "all Layer-A statistics are returned")
ok(is.data.table(ta), "returns a data.table")

cat("\n--- effect_sweep ---\n")
res <- effect_sweep(D, sample_col = "sample", x_col = "x", y_col = "y",
                    ct_col = "final_CT", disease_col = "disease_status",
                    exposures = list(published = "lab_published",
                                     naive = "rule_naive",
                                     aligned = "rule_aligned"),
                    reference = "ref_pop",
                    partner_defs = list(activated = c("FibroblastAct"),
                                        all = c("FibroblastAct","FibroblastAll")),
                    radii = c(10, 30, 100), min_cells = 3)
ok(nrow(res$sweep) > 0, "produces per-sample estimates")
ok(nrow(res$meta) > 0, "produces pooled estimates")
ok(all(c("OR","lo","hi","p","I2","k") %in% names(res$meta)), "pooled table is complete")
ok(all(res$meta$OR > 0), "all pooled odds ratios are positive")

cat("\n--- effect_sweep with strata override ---\n")
D2 <- copy(D); D2[, `:=`(s1 = n_SFTPC > 0, a1 = final_CT == "Epithelial")]
res2 <- effect_sweep(D2, sample_col = "sample", x_col = "x", y_col = "y",
                     ct_col = "final_CT",
                     exposures = list(s1 = "s1"),
                     reference = "ref_pop",
                     strata = list(list(exposures = "s1", reference = "a1")),
                     partner_defs = list(activated = c("FibroblastAct")),
                     radii = c(10, 30), min_cells = 3)
ok(nrow(res2$sweep) > 0, "strata override runs and produces rows")

cat("\n--- conclusion_flips ---\n")
fl <- conclusion_flips(res$meta)
ok(nrow(fl) == uniqueN(paste(res$meta$partner, res$meta$radius)), "one row per partner x radius")
ok(all(c("ratio_max_min","n_sig_enrich","n_sig_deplete","flips") %in% names(fl)),
   "flip summary columns present")

cat("\n--- platform_compare ---\n")
a <- copy(res$meta); a[, OR := OR * 1.1]
pc <- platform_compare(res$meta, a, "A", "B")
ok(!is.null(pc$by_definition), "returns a by-definition comparison")
ok(pc$median_platform_diff > 0, "computes a positive median platform difference")

cat("\n--- check_report ---\n")
bad <- data.table(label_in_marker_rate = 0.5, marker_out_rate = NA, R = NA,
                  mean_jaccard = 0.1, n_exposure_defs = 1, n_partner_defs = 1,
                  n_radii = 1, alignment_delta = NA, alignment_pass = FALSE)
ok(nrow(check_report(bad)) == 4, "flags all four violations in an incomplete report")
good <- data.table(label_in_marker_rate = 0.6, marker_out_rate = 0.05, R = 10,
                   mean_jaccard = 0.5, n_exposure_defs = 3, n_partner_defs = 2,
                   n_radii = 5, alignment_delta = 0.05, alignment_pass = TRUE)
ok(nrow(check_report(good)) == 0, "passes a complete report")
ok(inherits(try(check_report(data.table(a = 1)), silent = TRUE), "try-error"),
   "errors on a table with missing required columns")

cat("\n--- guard: nn2 must receive coordinates, not indices ---\n")
# passing indices as coordinates must NOT reproduce the coordinate-based result
Di <- copy(D); Di[, `:=`(x = as.numeric(.I), y = as.numeric(.I))]
ri <- effect_sweep(Di, sample_col = "sample", x_col = "x", y_col = "y",
                   ct_col = "final_CT",
                   exposures = list(published = "lab_published"),
                   reference = "ref_pop",
                   partner_defs = list(activated = c("FibroblastAct")),
                   radii = c(10), min_cells = 3)
ok(nrow(ri$meta) == 0 || !isTRUE(all.equal(ri$meta$OR, res$meta[exposure=="published" &
      partner=="activated" & radius==10, OR])),
   "index-as-coordinate input does not reproduce the real result (demonstrates the trap)")

cat(sprintf("\n==== %s: %d test(s) failed ====\n", ifelse(fails == 0, "ALL TESTS PASSED", "FAILURES"), fails))
if (fails > 0) quit(status = 1)
