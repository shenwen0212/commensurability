###############################################################################
# commensurability.R -- diagnostics for cell-type label commensurability
#
# Companion code for:
#   "Incommensurable cell-type annotations: an alignment diagnostic and its
#    effect on cross-platform spatial conclusions"
#
# Diagnostic layers:
#   label_diagnostics()  Layer A: overlap, discrimination, prevalence shift
#   effect_sweep()       Layer B: effect size across a definition matrix
#   platform_compare()   Layer C: does an aligned re-derivation rescue agreement?
# Reporting helpers:
#   conclusion_flips()   does the conclusion flip with the definition?
#   check_report()       audits an S1-S5 style reporting table
#
# Implementation contract (see Methods):
#  * nn2() takes COORDINATE matrices, never row indices. Passing indices yields
#    plausible-looking but meaningless distances without an error.
#  * the per-sample table must contain ALL cells (all lineages), otherwise the
#    partner cell type is absent and the sweep silently returns nothing.
#  * pooling is per sample followed by random-effects meta-analysis; pooling
#    cells directly makes p-values meaningless at single-cell n.
###############################################################################

suppressPackageStartupMessages({
  library(data.table); library(RANN); library(metafor)
})

#' Layer A -- label diagnostics
#'
#' @param dt data.table with one row per cell.
#' @param label character; name of the published annotation column (logical).
#' @param rules character vector; names of candidate gating-rule columns.
#' @param reference_rule optional character; a rule that IS the published label
#'   by construction (e.g. the annotation itself). Its Jaccard is undefined and
#'   returned as NA. Default NULL.
#' @return data.table with, per rule: n_rule, n_label, jaccard_vs_label,
#'   in/out positive rates, the discriminative ratio R, prevalence shift lambda,
#'   and frac_rule_in_label (the commonly reported but non-diagnostic statistic).
#' @export
label_diagnostics <- function(dt, label, rules, reference_rule = NULL) {
  stopifnot(is.data.table(dt), label %in% names(dt), all(rules %in% names(dt)))
  ann <- as.logical(dt[[label]])
  rbindlist(lapply(rules, function(r) {
    v <- as.logical(dt[[r]])
    tp <- sum(v & ann); fp <- sum(v & !ann); fn <- sum(!v & ann); tn <- sum(!v & !ann)
    p_in  <- tp / max(tp + fn, 1)
    p_out <- fp / max(fp + tn, 1)
    data.table(
      rule = r, n_rule = sum(v), n_label = sum(ann),
      jaccard_vs_label = if (!is.null(reference_rule) && identical(r, reference_rule))
        NA_real_ else tp / max(tp + fp + fn, 1),
      p_in_label = p_in, p_out_label = p_out,
      R = p_in / max(p_out, 1e-9),
      lambda = sum(v) / max(sum(ann), 1),
      frac_rule_in_label = p_in)
  }))
}

###############################################################################
# internal helpers
###############################################################################

.erode_hull <- function(pts, r) {
  if (nrow(pts) < 5) return(rep(TRUE, nrow(pts)))
  h <- chull(pts); poly <- pts[c(h, h[1]), , drop = FALSE]; n <- nrow(poly) - 1
  A <- 0; for (i in seq_len(n)) A <- A + (poly[i,1]*poly[i+1,2] - poly[i+1,1]*poly[i,2])
  sgn <- if (A > 0) 1 else -1
  ok <- rep(TRUE, nrow(pts))
  for (i in seq_len(n)) {
    x1 <- poly[i,1]; y1 <- poly[i,2]; x2 <- poly[i+1,1]; y2 <- poly[i+1,2]
    ex <- x2-x1; ey <- y2-y1; L <- sqrt(ex^2 + ey^2); if (L < 1e-9) next
    nx <- sgn*(-ey/L); ny <- sgn*(ex/L)
    ok <- ok & (nx*(pts[,1]-x1) + ny*(pts[,2]-y1) >= r)
  }
  ok
}

.or_se <- function(d_exp, d_ref, r) {
  a <- sum(d_exp <= r); b <- sum(d_ref <= r)
  nA <- length(d_exp); nB <- length(d_ref)
  if (nA < 1 || nB < 1) return(c(NA_real_, NA_real_))
  a2 <- a + .5; b2 <- b + .5; c2 <- (nA - a) + .5; d2 <- (nB - b) + .5
  c(log((a2*d2)/(b2*c2)), sqrt(1/a2 + 1/b2 + 1/c2 + 1/d2))
}

.as_flag <- function(z, dt) if (is.character(z)) as.logical(dt[[z]]) else as.logical(z)

#' Layer B -- effect sweep across a definition matrix
#'
#' @param dt data.table of ALL cells (every lineage present).
#' @param sample_col,x_col,y_col,ct_col column names.
#' @param exposures named list; each element is a logical vector or a column name.
#' @param reference logical vector or column name: the default reference.
#' @param partner_defs named list of character vectors of cell-type labels.
#' @param radii numeric vector of radii in the units of x/y.
#' @param disease_col optional column whose per-sample first value is recorded.
#' @param strata optional list of strata, each with elements `exposures`
#'   (character vector of exposure names) and `reference` (logical vector or
#'   column name) to override the default reference for those exposures.
#' @param erode logical; erode the sample hull by max(radii).
#' @param min_cells minimum cells required in exposure, reference and partner.
#' @return list(sweep = per-sample estimates, meta = pooled estimates)
#' @export
effect_sweep <- function(dt, sample_col, x_col, y_col, ct_col,
                         exposures, reference, partner_defs,
                         radii = c(10, 20, 30, 50, 100), disease_col = NULL,
                         strata = NULL, erode = TRUE, min_cells = 3) {
  stopifnot(is.data.table(dt))
  RMAX <- max(radii)
  ref_of <- function(nm) {
    if (is.null(strata)) return(reference)
    for (s in strata) if (nm %in% s$exposures) return(s$reference)
    reference
  }
  samples <- unique(dt[[sample_col]])
  sw <- list()
  for (s in samples) {
    A <- dt[get(sample_col) == s]
    if (nrow(A) < 50) next
    pts <- as.matrix(A[, c(x_col, y_col), with = FALSE])
    keep <- if (erode) .erode_hull(pts, RMAX) else rep(TRUE, nrow(A))
    dP <- list()
    for (pd in names(partner_defs)) {
      fi <- which(A[[ct_col]] %in% partner_defs[[pd]] & keep)
      dP[[pd]] <- if (length(fi) >= min_cells)
        nn2(pts[fi, , drop = FALSE], pts, k = 1, searchtype = "priority")$nn.dist[, 1] else NULL
    }
    dis <- if (!is.null(disease_col)) A[[disease_col]][1] else NA_character_
    for (en in names(exposures)) {
      ex_i <- which(.as_flag(exposures[[en]], A) & keep)
      ref_i <- which(.as_flag(ref_of(en), A) & keep)
      if (length(ex_i) < min_cells || length(ref_i) < min_cells) next
      for (pd in names(dP)) {
        if (is.null(dP[[pd]])) next
        for (r in radii) {
          o <- .or_se(dP[[pd]][ex_i], dP[[pd]][ref_i], r)
          if (is.na(o[1])) next
          sw[[length(sw) + 1]] <- data.table(
            sample = s, disease_status = dis, exposure = en, partner = pd, radius = r,
            n_exposure = length(ex_i), n_reference = length(ref_i),
            theta = o[1], se = o[2],
            median_exposure = median(dP[[pd]][ex_i]),
            median_reference = median(dP[[pd]][ref_i]))
        }
      }
    }
  }
  SW <- if (length(sw)) rbindlist(sw) else data.table()
  me <- list()
  if (nrow(SW)) {
    grid <- unique(SW[, .(exposure, partner, radius)])
    for (g in seq_len(nrow(grid))) {
      gr <- grid[g]
      d <- SW[exposure == gr$exposure & partner == gr$partner & radius == gr$radius]
      d <- d[!is.na(theta) & !is.na(se) & se > 0]
      if (nrow(d) < 3) next
      m <- tryCatch(rma(yi = theta, sei = se, data = d, method = "REML", test = "knha"),
                    error = function(e) NULL)
      if (is.null(m)) next
      me[[length(me) + 1]] <- data.table(
        exposure = gr$exposure, partner = gr$partner, radius = gr$radius, k = nrow(d),
        OR = exp(as.numeric(m$b)), lo = exp(m$ci.lb), hi = exp(m$ci.ub),
        p = m$pval, I2 = m$I2)
    }
  }
  list(sweep = SW, meta = if (length(me)) rbindlist(me) else data.table())
}

#' Layer B summary -- does the reported conclusion depend on the definition?
#' @export
conclusion_flips <- function(meta) {
  stopifnot(nrow(meta) > 0)
  meta[, {
    en <- sum(p < 0.05 & OR > 1); de <- sum(p < 0.05 & OR < 1)
    .(n_def = .N, OR_min = min(OR), OR_max = max(OR),
      ratio_max_min = max(OR) / min(OR),
      n_sig_enrich = en, n_sig_deplete = de,
      flips = (en > 0 && de > 0))
  }, by = .(partner, radius)]
}

#' Layer C -- compare platforms after aligning the exposure definition
#' @export
platform_compare <- function(meta_a, meta_b, label_a = "platform A", label_b = "platform B") {
  a <- copy(meta_a); b <- copy(meta_b)
  a[, platform := label_a]; b[, platform := label_b]
  common <- intersect(names(a), names(b))
  AB <- rbind(a[, ..common], b[, ..common])
  AB[, yi := log(OR)]
  cmp <- dcast(AB, partner + radius + exposure ~ platform,
               value.var = "yi", fun.aggregate = mean)
  pcols <- setdiff(names(cmp), c("partner", "radius", "exposure"))
  if (length(pcols) < 2) return(cmp)
  cmp[, platform_diff := abs(get(pcols[1]) - get(pcols[2]))]
  list(by_definition = cmp,
       median_platform_diff = median(cmp$platform_diff, na.rm = TRUE),
       max_platform_diff = max(cmp$platform_diff, na.rm = TRUE))
}

#' Audit a reporting table against the S1-S5 minimal reporting set
#' @export
check_report <- function(tab) {
  need <- c("label_in_marker_rate","marker_out_rate","R","mean_jaccard",
            "n_exposure_defs","n_partner_defs","n_radii","alignment_delta","alignment_pass")
  miss <- setdiff(need, names(tab))
  if (length(miss)) stop("missing columns: ", paste(miss, collapse = ", "))
  v <- list()
  add <- function(i, rule, msg) v[[length(v)+1]] <<- data.table(row = i, rule = rule, message = msg)
  for (i in seq_len(nrow(tab))) {
    r <- tab[i]
    if (!is.na(r$R) && r$R <= 1)
      add(i, "S1/S2", "R <= 1: the defining marker is not enriched in this label")
    if (is.na(r$marker_out_rate))
      add(i, "S2", "marker-out positive rate not reported: label-in rate alone is uninterpretable")
    if (is.na(r$mean_jaccard) || r$mean_jaccard < 0.3)
      add(i, "S3", sprintf("low overlap with marker rules (Jaccard %.3f): report the full ladder", r$mean_jaccard))
    if (r$n_exposure_defs < 2 || r$n_partner_defs < 2 || r$n_radii < 3)
      add(i, "S4", "definition matrix too small: need >=2 exposures x >=2 partners x >=3 radii")
    if (is.na(r$alignment_delta) || !isTRUE(r$alignment_pass))
      add(i, "S5", "cross-platform claim without a passing alignment check")
  }
  if (length(v)) rbindlist(v) else data.table(row = integer(), rule = character(), message = character())
}
