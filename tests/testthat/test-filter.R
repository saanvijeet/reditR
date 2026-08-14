library(testthat)
library(reditR)
library(data.table)

# Write a temp editing file and run filter_editing_sites() on it
make_temp_editing_file <- function(rows) {
  f <- tempfile(fileext = ".txt")
  write.table(rows, f, sep = "\t", row.names = FALSE, quote = FALSE)
  f
}

base_rows <- data.frame(
  site   = c("chr1:100", "chr1:100", "chr1:100",
             "chr1:150", "chr1:150", "chr1:150",
             "chr2:900"),
  sample = paste0("S", 1:7),
  edited = c(5L,  3L, 4L,  2L, 1L, 3L, 1L),
  total  = c(30L, 25L, 20L, 15L, 12L, 18L, 5L)
)

test_that("coverage filter removes rows below min_coverage", {
  f    <- make_temp_editing_file(base_rows)
  on.exit(unlink(f))
  res  <- filter_editing_sites(f, out_dir = NULL, min_coverage = 10, min_edited = 1, min_groups = 1)

  expect_true(all(res$all$total >= 10))
})

test_that("min_edited filter removes rows with too few edited reads", {
  f    <- make_temp_editing_file(base_rows)
  on.exit(unlink(f))
  res  <- filter_editing_sites(f, out_dir = NULL, min_coverage = 1, min_edited = 3, min_groups = 1)

  expect_true(all(res$all$edited >= 3))
})

test_that("min_groups filter removes underpresent sites", {
  single <- data.frame(site = "chr1:100", sample = "S1", edited = 5L, total = 30L)
  f      <- make_temp_editing_file(single)
  on.exit(unlink(f))
  res    <- filter_editing_sites(f, out_dir = NULL, min_coverage = 1, min_edited = 1, min_groups = 2)

  expect_equal(nrow(res$all), 0L)
})

test_that("clustered sites are a subset of all filtered sites", {
  f   <- make_temp_editing_file(base_rows)
  on.exit(unlink(f))
  res <- filter_editing_sites(f, out_dir = NULL, min_coverage = 10,
                               min_edited = 1, min_groups = 1, cluster_window = 100)

  expect_true(all(res$clustered$site %in% res$all$site))
})

test_that("output is written when out_dir is provided", {
  f      <- make_temp_editing_file(base_rows)
  tmpdir <- tempdir()
  on.exit(unlink(f))

  filter_editing_sites(f, out_dir = tmpdir, min_coverage = 10, min_edited = 1, min_groups = 1)

  expect_true(file.exists(file.path(tmpdir, "filtered_sites_all.txt")))
  expect_true(file.exists(file.path(tmpdir, "filtered_sites_clustered.txt")))
})

# ---- min_edit_ratio -------------------------------------------------------
# Rows span the 1% boundary: 1/200 = 0.005 (below), 2/200 = 0.01 (exactly on),
# 5/30 and 3/25 (well above). The boundary row must be KEPT, since the filter
# is >= not >.
ratio_rows <- data.frame(
  site   = c("chr1:100", "chr1:100", "chr1:100", "chr1:100"),
  sample = paste0("S", 1:4),
  edited = c(5L,  3L,  1L,   2L),
  total  = c(30L, 25L, 200L, 200L)
)

test_that("min_edit_ratio defaults to 0 and is a no-op", {
  f <- make_temp_editing_file(ratio_rows)
  on.exit(unlink(f))
  a <- filter_editing_sites(f, out_dir = NULL, min_coverage = 10, min_edited = 1, min_groups = 1)
  b <- filter_editing_sites(f, out_dir = NULL, min_coverage = 10, min_edited = 1, min_groups = 1,
                            min_edit_ratio = 0)

  expect_equal(nrow(a$all), 4L)
  expect_equal(a$all, b$all)
})

test_that("min_edit_ratio drops observations below the threshold, boundary inclusive", {
  f <- make_temp_editing_file(ratio_rows)
  on.exit(unlink(f))
  res <- filter_editing_sites(f, out_dir = NULL, min_coverage = 10, min_edited = 1,
                              min_groups = 1, min_edit_ratio = 0.01)

  expect_equal(nrow(res$all), 3L)
  expect_true(all(res$all$edit_ratio >= 0.01))
  expect_false("S3" %in% res$all$sample)   # 1/200 = 0.005, below
  expect_true("S4"  %in% res$all$sample)   # 2/200 = 0.010, exactly on
})

test_that("min_edit_ratio interacts with min_groups on the post-ratio count", {
  # S3 and S4 both fall below 2%, leaving one sample, so the site fails min_groups = 2
  f <- make_temp_editing_file(ratio_rows[3:4, ])
  on.exit(unlink(f))
  res <- filter_editing_sites(f, out_dir = NULL, min_coverage = 10, min_edited = 1,
                              min_groups = 2, min_edit_ratio = 0.02)

  expect_equal(nrow(res$all), 0L)
})

test_that("zero-coverage rows are not silently dropped when min_coverage is 0", {
  # total == 0 gives edit_ratio = NaN. Guarding the ratio filter behind
  # min_edit_ratio > 0 keeps the default path free of NA-driven row loss.
  zero_rows <- data.frame(
    site   = c("chr1:100", "chr1:100"),
    sample = c("S1", "S2"),
    edited = c(0L, 5L),
    total  = c(0L, 30L)
  )
  f <- make_temp_editing_file(zero_rows)
  on.exit(unlink(f))
  res <- filter_editing_sites(f, out_dir = NULL, min_coverage = 0, min_edited = 0, min_groups = 1)

  expect_equal(nrow(res$all), 2L)
  expect_true(is.nan(res$all$edit_ratio[res$all$sample == "S1"]))
})
