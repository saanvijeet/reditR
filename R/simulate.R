#' Simulate editing data with known ground truth
#'
#' Generates synthetic per-site editing counts under a mixed-effects generative
#' model with known effect sizes. Useful for power analysis and for validating
#' that \code{\link{differential_editing}} recovers planted effects.
#'
#' @param n_null Number of null (no-effect) sites. Default 800.
#' @param n_effects Named numeric vector where names are effect sizes (as
#'   character strings) and values are the number of sites at that effect size.
#'   Default \code{c("0.05" = 100, "0.10" = 100, "0.20" = 100)}.
#' @param n_per_condition Number of samples per condition. Default 6.
#' @param baseline_rate Baseline editing rate in the reference condition.
#'   Default 0.10.
#' @param mean_coverage Mean read coverage per site per sample drawn from
#'   \code{rnbinom(mu = mean_coverage, size = 5)}, floored at 10. Default 30.
#' @param sample_re_sd Standard deviation of the per-sample random effect on
#'   the logit scale. Default 0.3.
#' @param condition_labels Length-2 character vector naming the reference and
#'   case conditions, in that order. Default \code{c("control", "case")}. The
#'   simulated effect is applied to the second label, and these are the values
#'   to pass as \code{reference_level} and \code{case_level} when handing the
#'   result to \code{\link{differential_editing}}.
#' @param seed Random seed for reproducibility. Default 42.
#'
#' @return A named list with three \code{data.table} elements:
#'   \describe{
#'     \item{\code{editing}}{Per-site, per-sample read counts.}
#'     \item{\code{metadata}}{Sample metadata with \code{sample} and
#'       \code{condition} columns.}
#'     \item{\code{truth}}{Per-site ground truth with \code{site} and
#'       \code{true_effect} columns.}
#'   }
#'
#' @examples
#' sim <- simulate_editing_data(n_null = 50, n_effects = c("0.10" = 20), seed = 1)
#' str(sim)
#'
#' @export
simulate_editing_data <- function(n_null          = 800L,
                                   n_effects       = c("0.05" = 100, "0.10" = 100, "0.20" = 100),
                                   n_per_condition = 6L,
                                   baseline_rate   = 0.10,
                                   mean_coverage   = 30L,
                                   sample_re_sd    = 0.3,
                                   condition_labels = c("control", "case"),
                                   seed            = 42L) {
  # Convert the condition labels to character values so they can be used
  # consistently when creating the sample metadata and simulation.
  condition_labels <- as.character(condition_labels)
  # Check that exactly two different condition labels were supplied.
  # The first label is treated as the reference condition and the second
  # label as the condition where the simulated effect is introduced.
  if (length(condition_labels) != 2L || anyNA(condition_labels) ||
      condition_labels[1] == condition_labels[2])
    stop("condition_labels must be two distinct, non-missing values ",
         "(reference first, case second).")
         # Store the reference and case labels separately to make the simulation
  # code easier to read and avoid repeatedly indexing condition_labels.
  ref_label  <- condition_labels[1]
  case_label <- condition_labels[2]

  # Set the random seed so that the same inputs produce the same simulated data.
  # This makes the simulation reproducible and allows results to be checked later.
  set.seed(seed)

  # Create metadata describing which samples belong to each condition.
  # The reference condition is listed first so that the factor levels match
  # the intended reference-versus-case comparison.
  metadata <- data.table(
    sample    = c(paste0("ctrl_", seq_len(n_per_condition)),
                  paste0("case_", seq_len(n_per_condition))),
    condition = factor(
      c(rep(ref_label, n_per_condition), rep(case_label, n_per_condition)),
      levels = condition_labels
    )
  )

  # Create the ground-truth table for sites with no real condition effect.
  # These sites should not be detected as differentially edited and are
  # therefore used to measure the false-positive rate.
  null_sites <- data.table(
    site        = paste0("site_null_", seq_len(n_null)),
    true_effect = 0
  )
  # Create sites with known positive effects.
  # The effect size will later be added to the baseline editing rate
  # in the second (case) condition..
  effect_sites <- rbindlist(Map(function(eff_name, n) {
    data.table(
      site        = paste0("site_eff", eff_name, "_", seq_len(n)),
      true_effect = as.numeric(eff_name)
    )
  }, names(n_effects), n_effects))

  # Combine null and effect sites into one table containing the known
  # "ground truth" for every simulated site.
  truth       <- rbind(null_sites, effect_sites)
  # Store the sample names and conditions separately because they are
  # repeatedly used when generating counts for each site.
  all_samples <- metadata$sample
  conditions  <- as.character(metadata$condition)

  # Generate editing counts independently for each simulated site.
  rows <- lapply(seq_len(nrow(truth)), function(i) {
    # Retrieve the known effect size and site name for the current site.
    # The effect size determines whether this site should show a
    # difference between the two conditions.
    eff     <- truth$true_effect[i]
    site_id <- truth$site[i]

     # Generate a random sample-specific effect to represent variation
    # between biological samples. Larger sample_re_sd produces more
    # between-sample variability in editing rates.
    sample_re  <- rnorm(length(all_samples), 0, sample_re_sd)
    # Set the expected editing rate for each sample.
    # The reference condition stays at the baseline rate, while the
    # second condition receives the known simulated effect.
    mean_rates <- ifelse(conditions == case_label,
                         baseline_rate + eff,
                         baseline_rate)
    # Keep the expected rates within a valid probability range.
    # Values cannot be below 0 or above 1, and the small buffer away
    # from exactly 0 or 1 avoids problems on the logit scale.
    mean_rates <- pmax(0.001, pmin(0.999, mean_rates))
    # Convert the expected editing rates to the logit scale, add the
    # sample-specific random effect, and convert back to probabilities.
    # This allows each sample to deviate from the expected rate while
    # ensuring the final probability remains between 0 and 1.
    p          <- plogis(qlogis(mean_rates) + sample_re)
    # Simulate the total number of reads covering the site in each sample.
    # The negative binomial allows coverage to vary between samples,
    # reflecting the variable read depth seen in real RNA-seq data.
    # A minimum of 10 reads prevents unrealistically low coverage.
    total  <- pmax(10L, rnbinom(length(all_samples), mu = mean_coverage, size = 5))

    # Simulate the number of edited reads out of the total reads.
    # The probability p determines the expected proportion of reads
    # that are edited at this site in each sample.
    edited <- rbinom(length(all_samples), total, p)

    # Store the simulated counts and calculate the observed editing ratio.
    data.table(
      site       = site_id,
      sample     = all_samples,
      edited     = edited,
      total      = total,
      edit_ratio = edited / total
    )
  })

  # Combine the simulated data from all sites and return it together with
  # the sample metadata and the known ground truth.
  list(editing = rbindlist(rows), metadata = metadata, truth = truth)
}


#' Validate differential editing results against simulation ground truth
#'
#' Compares the output of \code{\link{differential_editing}} against the
#' ground truth from \code{\link{simulate_editing_data}}. Reports convergence
#' rate (GLMM), false-positive rate on null sites, and power by effect size --
#' independently for whichever test(s) are present in \code{results} (i.e.
#' whichever were requested via \code{differential_editing(test = ...)}).
#' There is no combined verdict to validate, matching
#' \code{differential_editing()}'s independent-tests design: each test is
#' reported on its own.
#'
#' @param results \code{data.table} returned by \code{differential_editing()}.
#' @param truth \code{data.table} with columns \code{site} and
#'   \code{true_effect}, as returned in \code{simulate_editing_data()$truth}.
#' @param fdr_threshold Retained for API consistency; unused -- the
#'   \verb{<Test>_sig} columns in \code{results} already reflect whichever
#'   FDR threshold \code{differential_editing()} was called with.
#'
#' @return An object of class \code{"reditR_validation"} -- a list with
#'   components \code{convergence} (numeric, GLMM only), \code{fpr} (named
#'   list, one element per test present in \code{results}), and \code{power}
#'   (\code{data.table} with one \verb{power_<Test>} column per test present).
#'
#' @examples
#' sim  <- simulate_editing_data(n_null = 30, n_effects = c("0.10" = 10), seed = 1)
#' \dontrun{
#' res  <- differential_editing(tempfile(), tempfile())  # use sim$editing / sim$metadata
#' val  <- validate_against_truth(res, sim$truth)
#' print(val)
#' }
#'
#' @export
validate_against_truth <- function(results, truth, fdr_threshold = 0.05) {
  # Join the analysis results to the known truth for each site.
  # This allows us to determine which sites truly had an effect and
  # which were simulated as null sites.
  merged <- merge(results, truth, by = "site", all.x = TRUE)

  # Calculate the proportion of GLMM sites that successfully produced
  # a p-value. Missing p-values indicate that the model did not produce
  # a valid result, for example because of convergence failure.
  convergence <- if ("glmm_pvalue" %in% names(merged))
    mean(!is.na(merged$glmm_pvalue)) else NA_real_

  # Define the significance columns associated with each statistical test.
  # Only tests that are actually present in the results are evaluated.
  possible_tests <- c(GLMM = "GLMM_sig", Fisher = "Fisher_sig", Wilcoxon = "Wilcox_sig")
  present_tests  <- possible_tests[possible_tests %in% names(merged)]

  # Select only sites that were simulated with no true effect.
  # These sites are used to assess how often the method produces
  # false-positive discoveries.
  null_sites <- merged[true_effect == 0]
  # Calculate the false-positive rate for each available test.
  # This is the proportion of truly null sites that were incorrectly
  # classified as significant.
  fpr <- lapply(present_tests, function(cn) {
    if (nrow(null_sites) == 0) return(NA_real_)
    mean(null_sites[[cn]], na.rm = TRUE)
  })
  names(fpr) <- names(present_tests)

  # Identify the different non-zero effect sizes that were deliberately
  # introduced into the simulated data.
  effect_sizes <- sort(unique(merged$true_effect[merged$true_effect > 0]))
  power <- data.table(effect_size = effect_sizes)
  # Calculate power separately for each statistical test.
  # Keeping the tests separate allows their performance to be compared
  # without combining them into a single overall result.
  for (nm in names(present_tests)) {
    cn <- present_tests[[nm]]
    power[[paste0("power_", nm)]] <- sapply(effect_sizes, function(eff) {
      # Select the sites that were simulated with this particular
      # effect size so their detection rate can be calculated.
      d <- merged[true_effect == eff]
      if (nrow(d) == 0) return(NA_real_)
      # Power is the proportion of truly affected sites that were
      # correctly identified as significant.
      mean(d[[cn]], na.rm = TRUE)
    })
  }
  # Store the validation results in a custom class.
  # This allows the object to have a dedicated print method for
  # displaying the results in a simple summary format.
  structure(
    list(convergence = convergence, fpr = fpr, power = power),
    class = "reditR_validation"
  )
}


#' Print a reditR validation summary
#'
#' @param x An object of class \code{"reditR_validation"}.
#' @param ... Ignored.
#'
#' @return \code{x} invisibly.
#' @export
print.reditR_validation <- function(x, ...) {
  cat("reditR Validation Summary\n")
  cat("=========================\n")
  if (!is.na(x$convergence))
    cat(sprintf("  Convergence rate (GLMM):        %.1f%%\n", 100 * x$convergence))
  cat("  False-positive rate (null sites):\n")
  for (nm in names(x$fpr)) {
    cat(sprintf("    %-10s%.3f\n", paste0(nm, ":"),
               if (is.na(x$fpr[[nm]])) NaN else x$fpr[[nm]]))
  }
  cat("  Power by effect size:\n")
  print(x$power)
  invisible(x)
}
