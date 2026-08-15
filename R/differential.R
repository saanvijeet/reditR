#' Run differential RNA editing analysis
#'
#' Identifies differentially RNA-edited (DRE) sites between two conditions
#' using up to three independent statistical tests: a binomial GLMM with a
#' per-sample random effect, Fisher's exact test on pooled per-sample counts,
#' and a Wilcoxon rank-sum test on per-sample editing ratios. All three tests
#' are run on ALL sites, independently of one another -- there is no gating
#' of one test on another's significance, and no combined verdict column.
#' Each requested test gets its own Benjamini-Hochberg FDR correction.
#' \code{test} selects which of the three are actually computed, since which
#' test is most appropriate depends on the experimental design (see
#' Details) -- this is deliberately a menu, not a pipeline that runs
#' everything and reduces it to one answer.
#'
#' @param data_path Path to the filtered and clustered editing-site count
#'   file. Required columns: \code{site}, \code{sample}, \code{edited},
#'   \code{total}. \code{edit_ratio} is computed if absent.
#' @param meta_path Path to the sample metadata file. Required columns:
#'   \code{sample}, \code{condition}. Ignored if \code{condition} is already
#'   present in the data file.
#' @param out_path File path for the tab-delimited results table. Pass
#'   \code{NULL} (default) to skip writing.
#' @param summary_path File path for a tab-delimited summary of key run
#'   metrics (sites tested, significant count per requested test). Pass
#'   \code{NULL} (default) to skip writing.
#' @param test Character vector selecting which test(s) to run. Any subset
#'   of \code{c("glmm", "fisher", "wilcoxon")}. Default runs all three. Only
#'   the requested test(s) are computed -- e.g. \code{test = "wilcoxon"}
#'   skips GLMM/Fisher entirely, useful when the design is known in advance
#'   to make GLMM unreliable (see Details).
#' @param reference_level Condition label for the reference/control group.
#'   Default \code{"control"}.
#' @param case_level Condition label for the case/disease group. Default
#'   \code{"diabetic"}.
#' @param random_effects Character string giving the random-effects side of
#'   the GLMM formula (everything after \code{condition +}), passed through
#'   to \code{lme4::glmer()}. Default \code{"(1 | sample)"} -- one random
#'   intercept per row, correct when each \code{sample} is one true
#'   biological replicate (bulk data). For single-cell pseudobulk data,
#'   where multiple pseudobulk units (e.g. cluster x library) share a
#'   library and/or a cluster identity, a single \code{sample}-level term
#'   pseudoreplicates -- pass crossed or nested terms instead, e.g.
#'   \code{"(1 | library) + (1 | cluster_id)"}, using whatever grouping
#'   columns are present in \code{data_path}. Ignored if \code{"glmm"} is
#'   not in \code{test}.
#' @param min_obs Minimum number of observations required to test a site
#'   with the GLMM. Default 4.
#' @param fdr_threshold FDR threshold applied to BH-adjusted p-values for
#'   every requested test. Default 0.05.
#' @param n_cores Number of cores for parallel per-site fitting (via
#'   \code{parallel::mclapply}). Default 1 (serial). Always serial on
#'   Windows, where \code{mclapply} cannot fork.
#' @param verbose Print progress and summary messages. Default \code{TRUE}.
#'
#' @return A \code{data.table} with column \code{site} plus, for each
#'   requested test, \verb{<test>_pvalue}, \verb{<Test>_FDR}, and
#'   \verb{<Test>_sig} (\code{TRUE} where FDR < \code{fdr_threshold}). No
#'   combined verdict column is produced -- the requested tests are
#'   independent options, not a consensus vote. Compare the columns you
#'   asked for directly, or intersect them yourself if you want an AND.
#'
#' @details
#' Which test to choose depends on the experimental design -- null-
#' calibration simulations (permutation nulls and parametric variance
#' sweeps) across three independent datasets found a consistent pattern:
#' \itemize{
#'   \item \strong{GLMM} (\code{cbind(edited, unedited) ~ condition +
#'     (1|sample)}) is the primary test (Srivastava et al. 2017), but can be
#'     severely anti-conservative when between-replicate variance is high
#'     relative to the number of independent control replicates -- in a
#'     confounded or very-small-n design this inflation reached ~9-10x the
#'     nominal false-positive rate. In a balanced, adequately-replicated
#'     design the same test was only mildly anti-conservative.
#'   \item \strong{Fisher's exact test} pools reads across samples per
#'     condition; independent of GLMM, but was anti-conservative to a
#'     similar degree whenever real between-replicate variance was present,
#'     regardless of dataset -- and unlike GLMM, adding replicates did not
#'     fix it.
#'   \item \strong{Wilcoxon rank-sum} on per-sample editing ratios was the
#'     most consistently well-calibrated test across every design checked,
#'     at the cost of lower power and a coarse p-value floor at small
#'     sample sizes (e.g. minimum ~0.05 at n=3 vs 3).
#' }
#' A design with few, unbalanced, or confounded replicates should lean on
#' Wilcoxon; a well-powered balanced design can reasonably use GLMM. Fisher
#' should be treated cautiously whenever real between-replicate variance is
#' plausible. Requesting more than one test lets you compare them directly
#' rather than trusting either in isolation.
#'
#' @examples
#' \dontrun{
#' res <- differential_editing(
#'   data_path = "filtered_sites_clustered.txt",
#'   meta_path = "sample_metadata.txt",
#'   test      = c("glmm", "wilcoxon"),
#'   out_path  = "DRE_results.txt"
#' )
#' res[GLMM_sig == TRUE]
#' res[GLMM_sig == TRUE & Wilcox_sig == TRUE]   # your own AND, if you want one
#' }
#'
#' @export
differential_editing <- function(data_path,
                                  meta_path       = NULL,
                                  out_path        = NULL,
                                  summary_path    = NULL,
                                  test            = c("glmm", "fisher", "wilcoxon"),
                                  reference_level = "control",
                                  case_level      = "diabetic",
                                  random_effects  = "(1 | sample)",
                                  min_obs         = 4L,
                                  fdr_threshold   = 0.05,
                                  n_cores         = 1L,
                                  verbose         = TRUE) {

  test <- match.arg(test, choices = c("glmm", "fisher", "wilcoxon"), several.ok = TRUE)

  data <- fread(data_path)
  if (!"condition" %in% names(data)) {
    if (is.null(meta_path))
      stop("'condition' column not found in data file and no meta_path supplied.")
    meta <- fread(meta_path)
    data <- merge(data, meta, by = "sample")
  }
  data <- data[condition %in% c(reference_level, case_level)]

  absent <- setdiff(c(reference_level, case_level), unique(as.character(data$condition)))
  if (length(absent) > 0)
    stop("condition level(s) not present in data: ", paste(absent, collapse = ", "),
         ". Check reference_level / case_level against the condition column.")

  # glmer coerces a character `condition` to a factor with ALPHABETICALLY
  # ordered levels, which is not necessarily reference_level first. coef_name
  # below is built as paste0("condition", case_level), so whenever case_level
  # sorted first the model's coefficient was named after reference_level
  # instead, the lookup missed, the error was swallowed by tryCatch, and every
  # site returned NA -- a completely null result indistinguishable from
  # universal convergence failure. Set the contrast explicitly.
  data[, condition := factor(as.character(condition),
                             levels = c(reference_level, case_level))]

  if (!"unedited" %in% names(data))   data[, unedited   := total - edited]
  if (!"edit_ratio" %in% names(data)) data[, edit_ratio := edited / total]

  if (verbose) {
    cat("Samples:   ", paste(sort(unique(data$sample)), collapse = " "), "\n")
    cat("Conditions:", reference_level, "=", sum(data$condition == reference_level),
        "obs |", case_level, "=", sum(data$condition == case_level), "obs\n")
    cat("Sites:     ", uniqueN(data$site), "\n")
    cat("Tests:     ", paste(test, collapse = ", "), "\n\n")
  }

  coef_name <- paste0("condition", case_level)

  glmm_formula <- NULL
  if ("glmm" %in% test) {
    glmm_formula <- stats::as.formula(
      paste("cbind(edited, unedited) ~ condition +", random_effects)
    )
    re_vars <- setdiff(all.vars(glmm_formula), c("edited", "unedited", "condition"))
    missing_re_vars <- setdiff(re_vars, names(data))
    if (length(missing_re_vars) > 0)
      stop("random_effects references column(s) not found in data: ",
           paste(missing_re_vars, collapse = ", "))
  }

  # ── Per-site test functions (closures over reference_level/case_level/min_obs) ──

  run_glmm <- function(d) {
    if (nrow(d) < min_obs || uniqueN(d$condition) < 2L) return(NULL)
    p <- tryCatch({
      m <- withCallingHandlers(
        lme4::glmer(glmm_formula,
                    data = d, family = binomial,
                    control = lme4::glmerControl(optimizer = "bobyqa",
                                                 optCtrl = list(maxfun = 2e5))),
        warning = function(w) invokeRestart("muffleWarning"))
      summary(m)$coefficients[coef_name, "Pr(>|z|)"]
    }, error = function(e) NA_real_)
    data.table(glmm_pvalue = p)
  }

  run_fisher <- function(d) {
    if (uniqueN(d$condition) < 2L) return(NULL)
    mat <- matrix(c(
      sum(d[condition == reference_level, edited]),
      sum(d[condition == case_level,      edited]),
      sum(d[condition == reference_level, unedited]),
      sum(d[condition == case_level,      unedited])
    ), nrow = 2, dimnames = list(c(reference_level, case_level), c("edited", "unedited")))
    p <- tryCatch(fisher.test(mat)$p.value, error = function(e) NA_real_)
    data.table(fisher_pvalue = p)
  }

  run_wilcox <- function(d) {
    if (uniqueN(d$condition) < 2L) return(NULL)
    p <- tryCatch(
      wilcox.test(edit_ratio ~ condition, data = d, exact = FALSE)$p.value,
      error = function(e) NA_real_)
    data.table(wilcox_pvalue = p)
  }

  data_list <- split(data, by = "site")
  cores <- if (.Platform$OS.type == "windows") 1L else n_cores

  run_one_test <- function(fn) {
    rbindlist(Filter(Negate(is.null), parallel::mclapply(data_list, fn, mc.cores = cores)),
              idcol = "site")
  }

  results <- data.table(site = names(data_list))

  if ("glmm" %in% test) {
    if (verbose) cat("[", format(Sys.time()), "] Running GLMM on", length(data_list), "sites...\n")
    t0 <- Sys.time()
    stage <- run_one_test(run_glmm)
    stage[, GLMM_FDR := p.adjust(glmm_pvalue, method = "BH")]
    stage[, GLMM_sig := !is.na(GLMM_FDR) & GLMM_FDR < fdr_threshold]
    results <- merge(results, stage, by = "site", all.x = TRUE)
    if (verbose) cat("  Done in", round(difftime(Sys.time(), t0, units = "mins"), 2),
                     "min. FDR<", fdr_threshold, ":", sum(stage$GLMM_sig, na.rm = TRUE), "\n\n")
  }

  if ("fisher" %in% test) {
    if (verbose) cat("[", format(Sys.time()), "] Running Fisher's exact on", length(data_list), "sites...\n")
    t0 <- Sys.time()
    stage <- run_one_test(run_fisher)
    stage[, Fisher_FDR := p.adjust(fisher_pvalue, method = "BH")]
    stage[, Fisher_sig := !is.na(Fisher_FDR) & Fisher_FDR < fdr_threshold]
    results <- merge(results, stage, by = "site", all.x = TRUE)
    if (verbose) cat("  Done in", round(difftime(Sys.time(), t0, units = "mins"), 2),
                     "min. FDR<", fdr_threshold, ":", sum(stage$Fisher_sig, na.rm = TRUE), "\n\n")
  }

  if ("wilcoxon" %in% test) {
    if (verbose) cat("[", format(Sys.time()), "] Running Wilcoxon on", length(data_list), "sites...\n")
    t0 <- Sys.time()
    stage <- run_one_test(run_wilcox)
    stage[, Wilcox_FDR := p.adjust(wilcox_pvalue, method = "BH")]
    stage[, Wilcox_sig := !is.na(Wilcox_FDR) & Wilcox_FDR < fdr_threshold]
    results <- merge(results, stage, by = "site", all.x = TRUE)
    if (verbose) cat("  Done in", round(difftime(Sys.time(), t0, units = "mins"), 2),
                     "min. FDR<", fdr_threshold, ":", sum(stage$Wilcox_sig, na.rm = TRUE), "\n\n")
  }

  if (!is.null(out_path)) fwrite(results, out_path, sep = "\t")

  sig_cols <- c(glmm = "GLMM_sig", fisher = "Fisher_sig", wilcoxon = "Wilcox_sig")[test]
  summary_dt <- data.table(
    metric = c("sites_tested", unname(sig_cols)),
    value  = c(nrow(results), vapply(sig_cols, function(cn) sum(results[[cn]], na.rm = TRUE), numeric(1)))
  )
  if (!is.null(summary_path)) fwrite(summary_dt, summary_path, sep = "\t")
  if (verbose) { cat("Summary:\n"); print(summary_dt) }

  results
}
