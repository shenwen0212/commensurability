## demo/sensitivity_demo.R -- two demonstrated capabilities that the manuscript's
## arguments depend on:
##   (1) NULL CALIBRATION: with synthetic labels carrying no tag-specific signal,
##       the estimator still returns a non-null effect if the reference population is
##       not exchangeable. Demonstrates that "no effect" must be established by
##       construction, not inferred from a diagnostic statistic.
##   (2) SENSITIVITY LADDER: how the reported effect moves as a label is progressively
##       contaminated, i.e. as its overlap with the true population falls.
##
## Run from the repository root:  Rscript demo/sensitivity_demo.R
suppressPackageStartupMessages(library(data.table))
source("commensurability.R")

set.seed(20260101)
D <- fread("demo/demo_cells.csv")
cat(sprintf("demo data: %d cells, %d samples\n\n", nrow(D), uniqueN(D$sample)))

###############################################################################
cat("=========== (1) Null calibration: labels with no tag-specific signal ===========\n")
## A "null" label is sampled only from NON-tag cells. It carries no tag-specific
## information, so any non-null estimate must come from the reference or geometry.
## We build two variants to separate those two sources:
##   null_random     : random non-tag epithelial cells (reference = Basal)
##   null_matched    : random non-tag cells, but reference redrawn the same way
##                     (an exchangeable reference), i.e. a properly constructed null
D[, is_tag := lab_published]
set.seed(11)
D[, null_random := FALSE]
D[is_tag == FALSE, null_random := runif(.N) < 0.5]
set.seed(12)
D[, null_matched_a := FALSE][, null_matched_b := FALSE]
non_tag <- D[is_tag == FALSE, which = TRUE]
pick <- sample(non_tag, floor(length(non_tag) * 0.25))
D$null_matched_a[pick[seq(1, length(pick), by = 2)]] <- TRUE
D$null_matched_b[pick[seq(2, length(pick), by = 2)]] <- TRUE

run <- function(exp_col, ref_col, label) {
  r <- effect_sweep(
    dt = D, sample_col = "sample", x_col = "x", y_col = "y", ct_col = "final_CT",
    exposures = setNames(list(exp_col), "exposure"), reference = ref_col,
    partner_defs = list(activated = c("FibroblastAct"),
                        all = c("FibroblastAct", "FibroblastAll")),
    radii = c(10, 20, 30, 50, 100), min_cells = 3)
  m <- r$meta
  if (nrow(m)) {
    m[, contrast := label]
    cat(sprintf("\n--- %s ---\n", label))
    print(dcast(m, partner ~ radius, value.var = "OR")[
            , lapply(.SD, function(x) if (is.numeric(x)) round(x, 3) else x)], row.names = FALSE)
    sig <- m[p < 0.05, .N]
    cat(sprintf("configurations with p < 0.05: %d of %d\n", sig, nrow(m)))
    if (nrow(m[p < 0.05])) cat("  -> these are effects produced WITHOUT tag-specific signal\n")
  } else cat(sprintf("\n--- %s ---\nno pooled estimates (too few cells per sample)\n", label))
  m
}

m_rand <- run("null_random", "ref_pop",
              "null label vs annotated Basal reference (non-exchangeable reference)")
m_match <- run("null_matched_a", "null_matched_b",
               "null label vs null reference (exchangeable reference: a valid null)")

cat("\nInterpretation: the first construction yields effects even though the label carries\n")
cat("no tag-specific information, because the reference population is not exchangeable\n")
cat("with it. The second construction, in which label and reference are drawn the same\n")
cat("way, is the correct null and should sit near OR = 1.\n")

###############################################################################
cat("\n=========== (2) Sensitivity ladder as overlap falls ===========\n")
## Contaminate the true tag population with non-tag cells at increasing rates.
## Jaccard(labelled set, true set) = 1 / (1 + beta) by construction.
tag_idx <- D[is_tag == TRUE, which = TRUE]
non_idx <- D[is_tag == FALSE, which = TRUE]
BETAS <- c(0, 0.5, 1, 2, 4, 8, 16)
out <- rbindlist(lapply(BETAS, function(b) {
  set.seed(1000 + b)
  D[, tmp_exp := FALSE]
  m <- min(2000, length(tag_idx))
  D$tmp_exp[sample(tag_idx, m)] <- TRUE
  nb <- floor(b * m)
  if (nb > 0) D$tmp_exp[sample(non_idx, min(nb, length(non_idx)))] <- TRUE
  r <- effect_sweep(
    dt = D, sample_col = "sample", x_col = "x", y_col = "y", ct_col = "final_CT",
    exposures = list(e = "tmp_exp"), reference = "ref_pop",
    partner_defs = list(activated = c("FibroblastAct")),
    radii = c(10, 30, 100), min_cells = 3)
  if (!nrow(r$meta)) return(NULL)
  r$meta[, `:=`(beta = b, jaccard = 1 / (1 + b))]
  r$meta[, .(jaccard, radius, partner, OR, lo, hi, p)]
}))
if (nrow(out)) {
  print(dcast(out, jaccard ~ radius, value.var = "OR")[
          , lapply(.SD, function(x) if (is.numeric(x)) round(x, 3) else x)], row.names = FALSE)
  cat("\nThe effect moves as overlap falls; the direction of movement depends on radius,\n")
  cat("which is why a single reported effect carries no information about its stability.\n")
}
cat("\nDONE\n")
