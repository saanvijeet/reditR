library(testthat)
library(reditR)
library(data.table)


# Helper function to create a temporary editing file for each test.
# Using temporary files keeps the tests independent and avoids modifying
# the package's example data or creating permanent test files.
make_temp_editing_file <- function(rows) {
  f <- tempfile(fileext = ".txt")
  write.table(rows, f, sep = "\t", row.names = FALSE, quote = FALSE)
  f
}


# Small example dataset used to test the main filtering criteria.
# It contains different coverage and edited-read counts so that each
# filtering rule can be tested independently.
base_rows <- data.frame(
  site   = c("chr1:100", "chr1:100", "chr1:100",
             "chr1:150", "chr1:150", "chr1:150",
             "chr2:900"),
  sample = paste0("S", 1:7),
  edited = c(5L,  3L, 4L,  2L, 1L, 3L, 1L),
  total  = c(30L, 25L, 20L, 15L, 12L, 18L, 5L)
)


# Test that observations with insufficient sequencing coverage are removed.
# A minimum coverage threshold helps avoid making editing estimates from
# sites supported by very few reads, where the observed editing proportion
# would be unreliable.
test_that("coverage filter removes rows below min_coverage", {
  f   <- make_temp_editing_file(base_rows)
  on.exit(unlink(f))

  res <- filter_editing_sites(
    f,
    out_dir = NULL,
    min_coverage = 10,
    min_edited = 1,
    min_groups = 1
  )

  # Every retained observation should have at least 10 total reads.
  expect_true(all(res$all$total >= 10))
})


# Test that observations with too few edited reads are removed.
# Requiring a minimum number of edited reads helps exclude very weak
# editing observations that may be more susceptible to sequencing noise.
test_that("min_edited filter removes rows with too few edited reads", {
  f   <- make_temp_editing_file(base_rows)
  on.exit(unlink(f))

  res <- filter_editing_sites(
    f,
    out_dir = NULL,
    min_coverage = 1,
    min_edited = 3,
    min_groups = 1
  )

  # Every retained observation should contain at least 3 edited reads.
  expect_true(all(res$all$edited >= 3))
})


# Test that sites represented in too few samples are removed.
# Requiring a minimum number of groups/samples prevents a site from being
# retained when there is insufficient replication to compare conditions.
test_that("min_groups filter removes underpresent sites", {

  # Create a site observed in only one sample.
  single <- data.frame(
    site = "chr1:100",
    sample = "S1",
    edited = 5L,
    total = 30L
  )

  f   <- make_temp_editing_file(single)
  on.exit(unlink(f))

  # Require the site to be observed in at least two samples.
  res <- filter_editing_sites(
    f,
    out_dir = NULL,
    min_coverage = 1,
    min_edited = 1,
    min_groups = 2
  )

  # The site should therefore be removed completely.
  expect_equal(nrow(res$all), 0L)
})


# Check that the clustered output is always a subset of the sites that
# passed the initial filtering criteria.
#
# Clustering is a secondary filtering step, so a site should not appear
# in the clustered results if it was already removed from the full results.
test_that("clustered sites are a subset of all filtered sites", {
  f   <- make_temp_editing_file(base_rows)
  on.exit(unlink(f))

  res <- filter_editing_sites(
    f,
    out_dir = NULL,
    min_coverage = 10,
    min_edited = 1,
    min_groups = 1,
    cluster_window = 100
  )

  # Every clustered site must also be present in the complete filtered set.
  expect_true(all(res$clustered$site %in% res$all$site))
})


# Check that filter_editing_sites() creates the expected output files
# when an output directory is supplied.
test_that("output is written when out_dir is provided", {
  f      <- make_temp_editing_file(base_rows)
  tmpdir <- tempdir()
  on.exit(unlink(f))

  filter_editing_sites(
    f,
    out_dir = tmpdir,
    min_coverage = 10,
    min_edited = 1,
    min_groups = 1
  )

  # Both the complete filtered results and the clustered subset
  # should be written to the specified directory.
  expect_true(file.exists(file.path(tmpdir, "filtered_sites_all.txt")))
  expect_true(file.exists(file.path(tmpdir, "filtered_sites_clustered.txt")))
})


# ---- min_edit_ratio -------------------------------------------------------
#
# Test the minimum editing-ratio filter using observations around the
# 1% threshold:
#
#   1/200 = 0.005  -> below threshold
#   2/200 = 0.010  -> exactly at threshold
#   5/30  = 0.167  -> above threshold
#   3/25  = 0.120  -> above threshold
#
# The exact boundary is important because the filter should be inclusive:
# observations with edit_ratio >= min_edit_ratio should be retained.
ratio_rows <- data.frame(
  site   = c("chr1:100", "chr1:100", "chr1:100", "chr1:100"),
  sample = paste0("S", 1:4),
  edited = c(5L, 3L, 1L, 2L),
  total  = c(30L, 25L, 200L, 200L)
)


# Check that the default minimum editing ratio of zero does not remove
# any observations that would otherwise pass the other filters.
test_that("min_edit_ratio defaults to 0 and is a no-op", {
  f <- make_temp_editing_file(ratio_rows)
  on.exit(unlink(f))

  # Run the filter using the default minimum ratio.
  a <- filter_editing_sites(
    f,
    out_dir = NULL,
    min_coverage = 10,
    min_edited = 1,
    min_groups = 1
  )

  # Explicitly setting the ratio to zero should give the same result.
  b <- filter_editing_sites(
    f,
    out_dir = NULL,
    min_coverage = 10,
    min_edited = 1,
    min_groups = 1,
    min_edit_ratio = 0
  )

  # All four observations should remain and both approaches
  # should produce identical results.
  expect_equal(nrow(a$all), 4L)
  expect_equal(a$all, b$all)
})


# Check that observations below the minimum editing ratio are removed
# while an observation exactly at the threshold is retained.
test_that("min_edit_ratio drops observations below the threshold, boundary inclusive", {
  f <- make_temp_editing_file(ratio_rows)
  on.exit(unlink(f))

  res <- filter_editing_sites(
    f,
    out_dir = NULL,
    min_coverage = 10,
    min_edited = 1,
    min_groups = 1,
    min_edit_ratio = 0.01
  )

  # Only the observation below 1% should be removed.
  expect_equal(nrow(res$all), 3L)

  # Every retained observation must meet or exceed the threshold.
  expect_true(all(res$all$edit_ratio >= 0.01))

  # S3 has an editing ratio of 1/200 = 0.005 and should be removed.
  expect_false("S3" %in% res$all$sample)

  # S4 has an editing ratio of 2/200 = 0.010 and is exactly at the
  # threshold, so it should be retained.
  expect_true("S4" %in% res$all$sample)
})


# Check how the editing-ratio filter interacts with the minimum number
# of samples required for a site.
#
# Both observations initially belong to the same site, but applying the
# ratio threshold removes one of them. The site is then represented by
# fewer samples than required by min_groups and should therefore be removed.
test_that("min_edit_ratio interacts with min_groups on the post-ratio count", {

  # Use only S3 and S4:
  # S3 = 1/200 = 0.005
  # S4 = 2/200 = 0.010
  f <- make_temp_editing_file(ratio_rows[3:4, ])
  on.exit(unlink(f))

  res <- filter_editing_sites(
    f,
    out_dir = NULL,
    min_coverage = 10,
    min_edited = 1,
    min_groups = 2,
    min_edit_ratio = 0.02
  )

  # Both observations are below 2%, so the site has no remaining
  # observations and cannot satisfy the minimum group requirement.
  expect_equal(nrow(res$all), 0L)
})


# Check that zero-coverage observations are handled safely when coverage
# filtering is disabled.
#
# With total = 0, the calculation edited / total produces NaN rather than
# a valid editing ratio. The ratio filter should only be applied when a
# positive min_edit_ratio is requested, preventing these rows from being
# unintentionally removed because of an NA/NaN comparison.
test_that("zero-coverage rows are not silently dropped when min_coverage is 0", {

  zero_rows <- data.frame(
    site   = c("chr1:100", "chr1:100"),
    sample = c("S1", "S2"),
    edited = c(0L, 5L),
    total  = c(0L, 30L)
  )

  f <- make_temp_editing_file(zero_rows)
  on.exit(unlink(f))

  res <- filter_editing_sites(
    f,
    out_dir = NULL,
    min_coverage = 0,
    min_edited = 0,
    min_groups = 1
  )

  # Both rows should remain because coverage filtering has been disabled
  # and no positive editing-ratio threshold was requested.
  expect_equal(nrow(res$all), 2L)

  # Confirm that the zero-coverage observation has the expected NaN
  # editing ratio rather than being silently removed.
  expect_true(is.nan(
    res$all$edit_ratio[res$all$sample == "S1"]
  ))
})