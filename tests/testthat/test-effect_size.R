library(testthat)
library(reditR)
library(data.table)


# Locate the example input files bundled with the reditR package.
# The example metadata uses "control" and "diabetic" as its condition labels.
# This does not make editing_difference() diabetes-specific: the function
# accepts case_level and reference_level as arguments and can therefore be
# applied to any two condition labels.
editing_file  <- system.file("extdata", "example_editing.txt",  package = "reditR")
metadata_file <- system.file("extdata", "example_metadata.txt", package = "reditR")


# Check that editing_difference() returns the expected data structure
# and the main columns required to describe the difference in editing.
test_that("editing_difference returns expected schema", {
  eff <- editing_difference(editing_file, meta_path = metadata_file)

  # Check that the result is returned as a data.table.
  expect_s3_class(eff, "data.table")

  # Check that the site identifier, condition-specific means and
  # editing difference are all present in the output.
  expect_true(all(c("site", "diabetic_mean", "control_mean",
                     "editing_difference") %in% names(eff)))
})


# The bundled metadata uses "control" and "diabetic", so these are the
# labels expected when the example data are used. The important behaviour
# being tested is that output column names are generated from the supplied
# condition labels, rather than that the function is restricted to diabetes.
test_that("output column names reflect the supplied condition labels", {
  eff <- editing_difference(editing_file, meta_path = metadata_file,
                             case_level = "diabetic",
                             reference_level = "control")

  # Check that each condition is given its own mean editing column.
  expect_true("diabetic_mean" %in% names(eff))
  expect_true("control_mean"  %in% names(eff))
})


# Check that the direction of the editing difference is correct.
# For chr1_100, editing is higher in the diabetic/case group than in the
# control/reference group, so the difference should be positive.
test_that("editing difference sign is positive for chr1_100 (diabetic > control)", {
  eff <- editing_difference(editing_file, meta_path = metadata_file)

  expect_gt(eff[site == "chr1_100", editing_difference], 0)
})


# Check that a site with very similar editing rates between conditions
# produces an editing difference close to zero.
test_that("chr1_200 editing difference is near zero", {
  eff <- editing_difference(editing_file, meta_path = metadata_file)

  # A small absolute difference indicates little change in editing
  # between the two conditions.
  expect_lt(abs(eff[site == "chr1_200", editing_difference]), 0.05)
})


# Check that sites are ranked by the magnitude of their editing difference
# rather than by its direction.
#
# Taking the absolute value means that a large decrease and a large increase
# are both treated as large effects when results are ranked.
test_that("sites are sorted by decreasing absolute editing difference", {
  eff   <- editing_difference(editing_file, meta_path = metadata_file)
  diffs <- abs(eff$editing_difference)

  # Confirm that absolute effect sizes decrease or remain equal
  # down the results table.
  expect_true(all(diff(diffs) <= 0))
})


# Check that results can be written to an output file when requested.
test_that("results are written to out_path when provided", {
  # Create a temporary file so that the test does not modify
  # any permanent files.
  out <- tempfile(fileext = ".txt")

  # Remove the temporary file after the test finishes.
  on.exit(unlink(out))

  # Run the function and request that the results are saved.
  editing_difference(editing_file,
                     meta_path = metadata_file,
                     out_path = out)

  # Confirm that the output file was successfully created.
  expect_true(file.exists(out))
})


# ---- negative-direction coverage -------------------------------------------
#
# The bundled example dataset contains a positive difference and a near-zero
# difference, but does not provide a clear decrease in editing. Therefore,
# it cannot test whether negative effects are handled correctly.
#
# A small synthetic dataset is created below with known editing proportions.
# This dataset is independent of the biological context of the bundled
# example data and is used to test the mathematical behaviour of the function.
#
# The synthetic data contain:
#   - a decrease: 0.60 -> 0.30, giving a difference of -0.30
#   - an increase: 0.20 -> 0.40, giving a difference of +0.20
#   - no change:  0.20 -> 0.20, giving a difference of  0.00
#
# This allows the expected direction and magnitude of the editing difference
# to be calculated in advance and compared directly with the function output.
#
# The editing ratios are identical across samples within each site, so the
# expected mean editing rates and differences are known exactly.
directional_rows <- data.frame(
  site = rep(c("site_down", "site_up", "site_flat"), each = 6),

  # Create three control and three diabetic samples for each site.
  # The labels are only used to match the metadata; they do not determine
  # how editing_difference() calculates the effect.
  sample = rep(c(paste0("ctrl_", 1:3), paste0("diab_", 1:3)), times = 3),

  # Number of edited reads.
  #
  # site_down: 30/50 = 0.60 in control and 15/50 = 0.30 in diabetic
  # site_up:   10/50 = 0.20 in control and 20/50 = 0.40 in diabetic
  # site_flat: 10/50 = 0.20 in both conditions
  edited = c(30L, 30L, 30L, 15L, 15L, 15L,
             10L, 10L, 10L, 20L, 20L, 20L,
             10L, 10L, 10L, 10L, 10L, 10L),

  # Keep the total number of reads constant so that the expected
  # editing proportions are straightforward to calculate.
  total = rep(50L, 18)
)


# Helper function to write the synthetic directional dataset to a
# temporary file in the same format expected by editing_difference().
# Using a helper avoids repeating the file-writing code in each test.
make_directional_file <- function() {
  f <- tempfile(fileext = ".txt")

  # Write the synthetic data as a tab-separated file.
  write.table(directional_rows, f,
              sep = "\t",
              row.names = FALSE,
              quote = FALSE)

  f
}


# Check that a genuine decrease is reported with a negative editing
# difference and that its magnitude and condition-specific means are correct.
test_that("editing difference is negative, and correct in magnitude, for a decrease", {
  f <- make_directional_file()

  # Remove the temporary file after the test finishes.
  on.exit(unlink(f))

  # Calculate editing differences using the synthetic data.
  eff <- editing_difference(f, meta_path = metadata_file)

  # A decrease from 0.60 to 0.30 should produce a negative difference.
  expect_lt(eff[site == "site_down", editing_difference], 0)

  # Check the exact expected difference:
  # 0.30 - 0.60 = -0.30.
  expect_equal(eff[site == "site_down", editing_difference],
               -0.30, tolerance = 1e-8)

  # Check that the control/reference editing mean is correct.
  expect_equal(eff[site == "site_down", control_mean],
               0.60, tolerance = 1e-8)

  # Check that the diabetic/case editing mean is correct.
  expect_equal(eff[site == "site_down", diabetic_mean],
               0.30, tolerance = 1e-8)
})


# Check that the function correctly reports both directions of change.
# A positive value represents an increase, while a negative value represents
# a decrease. A site with equal editing in both conditions should give zero.
test_that("increases and decreases are both reported with the correct sign", {
  f <- make_directional_file()
  on.exit(unlink(f))

  eff <- editing_difference(f, meta_path = metadata_file)

  # The increase from 0.20 to 0.40 should give +0.20.
  expect_equal(eff[site == "site_up", editing_difference],
               +0.20, tolerance = 1e-8)

  # Equal editing rates should give a difference of exactly zero.
  expect_equal(eff[site == "site_flat", editing_difference],
               0.00, tolerance = 1e-8)
})


# Check that sites are ranked by the absolute size of their effect.
#
# This is important because both increases and decreases can represent
# substantial differential editing. A large decrease should therefore rank
# above a smaller increase.
test_that("ordering is by absolute effect, so a large decrease outranks a smaller increase", {
  f <- make_directional_file()
  on.exit(unlink(f))

  eff <- editing_difference(f, meta_path = metadata_file)

  # The expected order is:
  # |−0.30| > |+0.20| > |0.00|.
  #
  # If the function sorted by the raw difference instead of its absolute
  # value, the positive increase would incorrectly rank above the decrease.
  expect_equal(eff$site,
               c("site_down", "site_up", "site_flat"))

  # Confirm that the absolute effect sizes decrease down the table.
  expect_true(all(diff(abs(eff$editing_difference)) <= 0))
})