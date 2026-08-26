# Load testthat to provide functions for writing automated unit tests.
library(testthat)

# Load the package functions being tested.
library(reditR)

# Load data.table for reading and manipulating the test data.
library(data.table)


# Locate the example input files included with the reditR package.
# system.file() makes the tests independent of the user's working directory.
editing_file  <- system.file("extdata", "example_editing.txt",  package = "reditR")
metadata_file <- system.file("extdata", "example_metadata.txt", package = "reditR")


# Test the default behaviour of differential_editing().
# By default, all three statistical tests should be run and their
# corresponding p-values, FDR values and significance indicators returned.
test_that("differential_editing runs all three tests by default with expected column schema", {

  # Run the complete differential editing analysis using the example data.
  res <- differential_editing(editing_file, metadata_file, verbose = FALSE)

  # Check that the result has the expected data.table structure.
  expect_s3_class(res, "data.table")

  # Define the columns that should be produced by the three tests.
  expected <- c("site",
               "glmm_pvalue",   "GLMM_FDR",   "GLMM_sig",
               "fisher_pvalue", "Fisher_FDR", "Fisher_sig",
               "wilcox_pvalue", "Wilcox_FDR", "Wilcox_sig")

  # Check that all expected result columns are present.
  expect_true(all(expected %in% names(res)))

  # DRE was part of an earlier approach and should no longer be returned.
  # This prevents old functionality from accidentally reappearing.
  expect_false("DRE" %in% names(res))
})


# Test that the test= argument correctly controls which statistical
# methods are run and which result columns are produced.
test_that("test= restricts which columns are computed", {

  # Request only the GLMM.
  res_glmm <- differential_editing(editing_file, metadata_file,
                                   test = "glmm", verbose = FALSE)

  # GLMM results should be present.
  expect_true(all(c("glmm_pvalue", "GLMM_sig") %in% names(res_glmm)))

  # Results from the other tests should not be generated.
  expect_false(any(c("fisher_pvalue", "wilcox_pvalue") %in% names(res_glmm)))


  # Request Fisher's exact test and Wilcoxon, but not the GLMM.
  res_wf <- differential_editing(editing_file, metadata_file,
                                 test = c("wilcoxon", "fisher"), verbose = FALSE)

  # The two requested tests should be present.
  expect_true(all(c("wilcox_pvalue", "fisher_pvalue") %in% names(res_wf)))

  # The GLMM should not have been run.
  expect_false("glmm_pvalue" %in% names(res_wf))
})


# Check that a site with a deliberately clear difference between
# conditions is detected by the GLMM.
test_that("chr1_100 (clear effect) is GLMM significant", {
  res <- differential_editing(editing_file, metadata_file, verbose = FALSE)

  # The site should be classified as significant after FDR correction.
  expect_true(res[site == "chr1_100", GLMM_sig])
})


# Check that a site with no deliberate difference is not incorrectly
# classified as significant by the GLMM.
test_that("chr1_200 (no effect) is not GLMM significant", {
  res <- differential_editing(editing_file, metadata_file, verbose = FALSE)

  # This provides a simple check that the test can distinguish
  # an obvious null site from a site with a real effect.
  expect_false(res[site == "chr1_200", GLMM_sig])
})


# Check that each statistical test is analysed independently.
# Fisher and Wilcoxon should be calculated for all eligible sites rather
# than only for sites that were significant in the GLMM.
# FDR correction should also be performed separately for each test.
test_that("each requested test is independently FDR-corrected (no gating, no combined verdict)", {
  res <- differential_editing(editing_file, metadata_file, verbose = FALSE)

  # Wilcoxon and Fisher should be calculated for every site with
  # observations in both conditions, rather than being restricted
  # to GLMM-significant sites.
  expect_equal(sum(!is.na(res$wilcox_pvalue)), nrow(res))
  expect_equal(sum(!is.na(res$fisher_pvalue)), nrow(res))

  # Verify that each test's p-values are independently adjusted using
  # the Benjamini-Hochberg method to control the false discovery rate.
  expect_equal(res$GLMM_FDR,   p.adjust(res$glmm_pvalue,   method = "BH"))
  expect_equal(res$Fisher_FDR, p.adjust(res$fisher_pvalue, method = "BH"))
  expect_equal(res$Wilcox_FDR, p.adjust(res$wilcox_pvalue, method = "BH"))
})


# Check that results can be written to an output file when requested.
test_that("results are written to out_path when provided", {

  # Create a temporary output file so the test does not modify
  # any permanent files in the package or user's system.
  out <- tempfile(fileext = ".txt")

  # Delete the temporary file after the test has finished.
  on.exit(unlink(out))

  # Run the analysis and request that the results are written to the file.
  differential_editing(editing_file, metadata_file,
                       out_path = out, verbose = FALSE)

  # Confirm that the output file was actually created.
  expect_true(file.exists(out))

  # Read the written file back in and check that it contains
  # the expected number of sites.
  written <- fread(out)
  expect_equal(nrow(written), 2L)
})


# Check that the optional summary output is written correctly.
test_that("summary file is written with expected rows when summary_path is provided", {

  # Create a temporary file for the summary.
  sumf <- tempfile(fileext = ".txt")

  # Remove the temporary file when the test finishes.
  on.exit(unlink(sumf))

  # Run the analysis and request a summary file.
  differential_editing(editing_file, metadata_file,
                       summary_path = sumf, verbose = FALSE)

  # Confirm that the summary file was created.
  expect_true(file.exists(sumf))

  # Read the summary back in for checking.
  sumdt <- fread(sumf)

  # A summary should contain metric names and their corresponding values.
  expect_true(all(c("metric", "value") %in% names(sumdt)))

  # Check that the three requested statistical tests are represented
  # in the summary.
  expect_true("GLMM_sig" %in% sumdt$metric)
  expect_true("Fisher_sig" %in% sumdt$metric)
  expect_true("Wilcox_sig" %in% sumdt$metric)

  # The summary should report only the requested tests. Any other metric
  # appearing here would mean an unrequested result had reached the output.
  expect_false("monte_carlo_pvalue" %in% sumdt$metric)
})


# Test that the random-effects argument can represent crossed random effects.
# This is important for pseudobulk data where observations can share both
# a library and a cell-cluster identity.
test_that("random_effects supports crossed terms for pseudobulk data", {

  # Make the simulated pseudobulk data reproducible.
  set.seed(42)

  # Define four libraries, with two libraries from each condition.
  libraries  <- c("lib_c1", "lib_c2", "lib_d1", "lib_d2")
  conditions <- c("control", "control", "dehydrated", "dehydrated")

  # Define three cell clusters that are present within each library.
  clusters <- paste0("cluster_", 1:3)

  # Generate one site's pseudobulk observations across libraries and clusters.
  rows <- do.call(rbind, lapply(seq_along(libraries), function(i) {

    # Get the condition associated with the current library.
    cond <- conditions[i]

    data.frame(
      site       = "chr1_1",
      sample     = paste0(libraries[i], "_", clusters),
      library    = libraries[i],
      cluster_id = clusters,
      condition  = cond,

      # Simulate more editing in the dehydrated condition so that
      # the GLMM has a known difference to detect.
      edited     = if (cond == "dehydrated") rpois(3, 8) + 2 else rpois(3, 2),

      # Keep total coverage constant for this simple test dataset.
      total      = 40,
      stringsAsFactors = FALSE
    )
  }))

  # Write the simulated pseudobulk data to a temporary file because
  # differential_editing() expects a file as input.
  ed_file <- tempfile(fileext = ".txt")
  on.exit(unlink(ed_file))
  write.table(rows, ed_file, sep = "\t", row.names = FALSE, quote = FALSE)

  # Fit a GLMM containing separate random intercepts for library and cluster.
  # These crossed random effects account for the fact that observations
  # can share the same library or cluster.
  res <- differential_editing(
    ed_file,
    test            = "glmm",
    random_effects  = "(1 | library) + (1 | cluster_id)",
    reference_level = "control",
    case_level      = "dehydrated",
    min_obs         = 4L,
    verbose         = FALSE
  )

  # Confirm that the GLMM ran and returned a valid p-value.
  expect_s3_class(res, "data.table")
  expect_true("glmm_pvalue" %in% names(res))
  expect_false(is.na(res[site == "chr1_1", glmm_pvalue]))
})


# Check that a random-effect variable must actually exist in the input data.
# The error should occur before attempting to fit the model.
test_that("random_effects referencing a missing column errors before any fitting", {

  # Create a minimal dataset that deliberately does not contain
  # the random-effect column requested below.
  ed_file <- tempfile(fileext = ".txt")
  on.exit(unlink(ed_file))

  write.table(
    data.frame(site = "chr1_1",
               sample = c("a", "b", "c", "d"),
               condition = c("control", "control", "dehydrated", "dehydrated"),
               edited = c(1, 2, 8, 9),
               total = c(20, 20, 20, 20)),
    ed_file, sep = "\t", row.names = FALSE, quote = FALSE
  )

  # Check that a clear error is produced rather than allowing the model
  # fitting stage to fail later with a less informative message.
  expect_error(
    differential_editing(ed_file, test = "glmm",
                         random_effects  = "(1 | nonexistent_col)",
                         reference_level = "control",
                         case_level = "dehydrated",
                         verbose = FALSE),
    "random_effects references column"
  )
})


# Check that the default random-effect specification is equivalent
# to explicitly requesting a random intercept for sample.
test_that("default random_effects is unchanged and still matches (1 | sample)", {

  # Run the GLMM using the package default.
  res <- differential_editing(editing_file, metadata_file,
                              test = "glmm", verbose = FALSE)

  # Run the same analysis while explicitly specifying the default model.
  res_explicit <- differential_editing(editing_file, metadata_file,
                                       test = "glmm",
                                       random_effects = "(1 | sample)",
                                       verbose = FALSE)

  # The two approaches should produce identical GLMM p-values.
  expect_equal(res$glmm_pvalue, res_explicit$glmm_pvalue)
})


# ---- condition contrast direction -----------------------------------------

# Test that the explicitly supplied reference and case levels are respected
# even when their names would be ordered differently alphabetically.
#
# This is important because the direction of a GLMM coefficient depends on
# which condition is treated as the reference. The statistical significance
# should not change simply because the condition names are alphabetically
# ordered differently.
test_that("reference_level is honoured when the case level sorts first alphabetically", {

  # Read the original metadata.
  meta <- data.table::fread(metadata_file)

  # Rename the conditions so the intended case level ("acase") comes
  # before the reference level ("zref") alphabetically.
  remap <- c(control = "zref", diabetic = "acase")
  meta[, condition := as.character(remap[condition])]

  # Write the modified metadata to a temporary file.
  mf <- tempfile(fileext = ".txt")
  data.table::fwrite(meta, mf, sep = "\t")
  on.exit(unlink(mf))

  # Run the analysis using the original condition labels.
  ref <- differential_editing(editing_file, metadata_file,
                              test = "glmm", verbose = FALSE)

  # Run it again using the renamed conditions and explicitly specify
  # which level should be treated as the reference and case.
  got <- differential_editing(editing_file, mf,
                              test = "glmm",
                              reference_level = "zref",
                              case_level = "acase",
                              verbose = FALSE)

  # The GLMM should still produce valid p-values rather than returning
  # NA because of an incorrect coefficient lookup.
  expect_false(all(is.na(got$glmm_pvalue)))

  # Relabelling the conditions does not change the underlying division
  # of samples into groups, so the significance of the contrast should
  # remain unchanged. Only the direction/sign of the estimated effect changes.
  expect_equal(got$glmm_pvalue, ref$glmm_pvalue)
  expect_equal(got$GLMM_sig,    ref$GLMM_sig)
})


# Check that the non-GLMM tests are also unaffected by the order in
# which the condition names would be sorted.
test_that("all three tests are unaffected by which level sorts first", {

  # Rename the conditions so that the case level sorts first alphabetically.
  meta <- data.table::fread(metadata_file)
  remap <- c(control = "zref", diabetic = "acase")
  meta[, condition := as.character(remap[condition])]

  # Save the renamed metadata to a temporary file.
  mf <- tempfile(fileext = ".txt")
  data.table::fwrite(meta, mf, sep = "\t")
  on.exit(unlink(mf))

  # Run the original analysis.
  ref <- differential_editing(editing_file, metadata_file,
                              verbose = FALSE)

  # Run the same analysis after renaming the conditions and explicitly
  # specifying the reference and case levels.
  got <- differential_editing(editing_file, mf,
                              reference_level = "zref",
                              case_level = "acase",
                              verbose = FALSE)

  # Fisher and Wilcoxon depend on the group membership rather than
  # the alphabetical names of the conditions, so their p-values should
  # remain unchanged.
  expect_equal(got$fisher_pvalue, ref$fisher_pvalue)
  expect_equal(got$wilcox_pvalue, ref$wilcox_pvalue)
})


# Check that invalid condition labels are detected early.
# This prevents the function from silently returning NA results when
# the requested comparison does not exist in the input data.
test_that("a condition level absent from the data errors instead of returning all NA", {

  # Deliberately request a case level that does not occur in the metadata.
  expect_error(
    differential_editing(editing_file, metadata_file,
                         test = "glmm",
                         case_level = "not_a_condition",
                         verbose = FALSE),
    "condition level\\(s\\) not present in data"
  )
})