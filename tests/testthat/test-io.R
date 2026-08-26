library(testthat)
library(reditR)
library(data.table)


# Locate the example editing and metadata files bundled with the package.
# These files provide a small, consistent dataset for testing the input
# functions without relying on external data.
editing_file  <- system.file("extdata", "example_editing.txt",
                             package = "reditR")
metadata_file <- system.file("extdata", "example_metadata.txt",
                              package = "reditR")


# Check that the editing table is read into the expected format and that
# edit_ratio is calculated when it is not already present in the input.
#
# The editing ratio represents the proportion of reads showing editing:
#
#     edit_ratio = edited reads / total reads
#
# Calculating this during input ensures that downstream filtering and
# statistical analysis have a consistent variable available.
test_that("read_editing_table reads correctly and computes edit_ratio if absent", {

  dt <- read_editing_table(editing_file)

  # The function should return a data.table so that the rest of the
  # package can efficiently manipulate the editing data.
  expect_s3_class(dt, "data.table")

  # Check that all columns required for downstream analysis are present.
  expect_true(all(c("site", "sample", "edited", "total", "edit_ratio")
                  %in% names(dt)))

  # Verify the edit_ratio calculation using a known observation:
  # 2 edited reads out of 30 total reads = 2/30.
  expect_equal(
    dt[site == "chr1_100" & sample == "ctrl_1", edit_ratio],
    2 / 30
  )
})


# Check that the function gives a clear error when essential input columns
# are missing. This prevents invalid data from progressing into later stages
# of the pipeline and producing less informative errors.
test_that("read_editing_table errors on missing required columns", {

  # Create a deliberately incomplete editing table containing only some
  # of the columns required by read_editing_table().
  f <- tempfile(fileext = ".txt")

  fwrite(
    data.table(site = "s1", edited = 1L),
    f,
    sep = "\t"
  )

  # Remove the temporary test file after the test finishes.
  on.exit(unlink(f))

  # The function should identify the missing required columns rather than
  # attempting to process an invalid input table.
  expect_error(
    read_editing_table(f),
    "Missing required columns"
  )
})


# Check that a clear error is returned when the input file does not exist.
# Explicit file validation is important because these functions are the
# entry point for data supplied to the downstream analysis pipeline.
test_that("read_editing_table errors on nonexistent file", {

  expect_error(
    read_editing_table("/nonexistent/path.txt"),
    "File not found"
  )
})


# Check that the metadata file is read into the expected format and contains
# the columns required to associate each sample with its experimental condition.
test_that("read_metadata reads correctly", {

  mt <- read_metadata(metadata_file)

  # Metadata should also be returned as a data.table for consistency
  # with the editing data.
  expect_s3_class(mt, "data.table")

  # These columns are required to identify samples and their conditions.
  expect_true(all(c("sample", "condition") %in% names(mt)))

  # The example metadata contains three control and three diabetic samples.
  # Checking the number of rows confirms that all samples were read.
  expect_equal(nrow(mt), 6L)
})


# Check that the specified reference condition is placed first in the
# factor levels.
#
# This is important for statistical modelling because the reference level
# determines how the condition contrast is interpreted. Explicitly setting
# the reference level avoids relying on the default factor ordering.
test_that("read_metadata with reference_level sets factor levels correctly", {

  mt <- read_metadata(
    metadata_file,
    reference_level = "control"
  )

  # The condition column should be converted to a factor for use in
  # downstream statistical models.
  expect_true(is.factor(mt$condition))

  # Confirm that the requested reference condition is the first level.
  expect_equal(levels(mt$condition)[1], "control")
})


# Check that an invalid reference condition produces an informative error.
# Without this validation, an incorrectly specified reference level could
# lead to an incorrect or ambiguous interpretation of the condition contrast
# in downstream analyses.
test_that("read_metadata errors when reference_level not in conditions", {

  expect_error(
    read_metadata(
      metadata_file,
      reference_level = "missing"
    ),
    "not found"
  )
})