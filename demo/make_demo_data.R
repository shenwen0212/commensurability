## make_demo_data.R -- generate small synthetic datasets with the SAME schema as the
## real inputs, so the diagnostic pipeline can be run end-to-end without downloading
## any primary data. Two objects are written:
##   demo_cells.csv   : the per-sample working table (all lineages present)
##   demo_exposures.csv: logical exposure/reference columns aligned to demo_cells
suppressPackageStartupMessages(library(data.table))

set.seed(20260101)
outdir <- "demo"
if (!dir.exists(outdir)) dir.create(outdir, recursive = TRUE)

## Two synthetic samples laid out on separate 1 x 1 mm fields.
n_per <- 4000
mk <- function(sample, patient, status, ox) {
  n <- n_per
  # partner cells (fibroblasts) clustered in a band; two subtypes with different
  # spatial affinity for the "tag" population, which is what the diagnostic detects
  ct  <- sample(c("Basal", "TagLabel", "FibroblastAct", "FibroblastAll", "Epithelial"),
                n, replace = TRUE, prob = c(.18, .08, .06, .16, .52))
  x <- runif(n, 0, 1000) + ox
  y <- runif(n, 0, 1000)
  # activated fibroblasts sit slightly closer to TagLabel than to Basal,
  # but only when the tag is defined by the marker rule; the published label is
  # placed differently, which is the effect the diagnostic must recover
  fb_act <- which(ct == "FibroblastAct")
  tag    <- which(ct == "TagLabel")
  basa   <- which(ct == "Basal")
  if (length(fb_act)) {
    x[fb_act] <- pmin(1000, pmax(0, x[fb_act] - ox + 500 + rnorm(length(fb_act), 0, 220)))
    y[fb_act] <- pmin(1000, pmax(0, 500 + rnorm(length(fb_act), 0, 220)))
    x[fb_act] <- x[fb_act] + ox
  }
  data.table(cell = seq_len(n) + ox * n, sample = sample, patient = patient,
             disease_status = status, final_CT = ct, x = x, y = y)
}
D <- rbind(mk("S1", "P1", "Disease", 0), mk("S2", "P2", "Disease", 1))
D <- rbind(D, mk("S3", "P3", "Control", 2), mk("S4", "P4", "Control", 3))

## Marker counts: the published "TagLabel" is only loosely related to the marker
## that names it (mirroring the real data), while a gating rule is highly specific.
D[, `:=`(
  n_KRT17 = fifelse(final_CT == "TagLabel", rbinom(.N, 3, 0.55),
             fifelse(final_CT == "Basal",   rbinom(.N, 3, 0.65),
             rbinom(.N, 3, 0.02))),
  n_KRT5  = fifelse(final_CT %in% c("Basal", "TagLabel"), rbinom(.N, 3, 0.35),
             rbinom(.N, 3, 0.01)),
  n_TP63  = fifelse(final_CT %in% c("Basal", "TagLabel"), rbinom(.N, 2, 0.40),
             rbinom(.N, 2, 0.01)),
  n_SFTPC = fifelse(final_CT == "Epithelial", rbinom(.N, 5, 0.70), rbinom(.N, 5, 0.05)))]
D[, `:=`(
  lab_published = final_CT == "TagLabel",
  rule_naive    = n_KRT17 > 0,
  rule_aligned  = n_KRT17 > 0 & n_KRT5 == 0,
  rule_strict   = n_KRT17 >= 2 & n_KRT5 == 0,
  ref_pop       = final_CT == "Basal")]

fwrite(D, file.path(outdir, "demo_cells.csv"))
cat("wrote", file.path(outdir, "demo_cells.csv"), "with", nrow(D), "cells\n")
print(D[, .N, by = .(sample, disease_status, final_CT)][order(sample, final_CT)])
