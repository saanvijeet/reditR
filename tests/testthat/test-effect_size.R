library(testthat)
library(reditR)
library(data.table)

editing_file  <- system.file("extdata", "example_editing.txt",  package = "reditR")
metadata_file <- system.file("extdata", "example_metadata.txt", package = "reditR")

test_that("editing_difference returns expected schema", {
  eff <- editing_difference(editing_file, meta_path = metadata_file)

  expect_s3_class(eff, "data.table")
  expect_true(all(c("site", "diabetic_mean", "control_mean",
                     "editing_difference") %in% names(eff)))
})

test_that("output column names reflect the supplied condition labels", {
  eff <- editing_difference(editing_file, meta_path = metadata_file,
                             case_level = "diabetic", reference_level = "control")
  expect_true("diabetic_mean" %in% names(eff))
  expect_true("control_mean"  %in% names(eff))
})

test_that("editing difference sign is positive for chr1_100 (diabetic > control)", {
  eff <- editing_difference(editing_file, meta_path = metadata_file)

  expect_gt(eff[site == "chr1_100", editing_difference], 0)
})

test_that("chr1_200 editing difference is near zero", {
  eff <- editing_difference(editing_file, meta_path = metadata_file)

  expect_lt(abs(eff[site == "chr1_200", editing_difference]), 0.05)
})

test_that("sites are sorted by decreasing absolute editing difference", {
  eff   <- editing_difference(editing_file, meta_path = metadata_file)
  diffs <- abs(eff$editing_difference)

  expect_true(all(diff(diffs) <= 0))
})

test_that("results are written to out_path when provided", {
  out <- tempfile(fileext = ".txt")
  on.exit(unlink(out))

  editing_difference(editing_file, meta_path = metadata_file, out_path = out)

  expect_true(file.exists(out))
})

# ---- negative-direction coverage -------------------------------------------
# The bundled fixture contains only an increase (chr1_100, +0.17) and a flat
# site (chr1_200, -0.001). Every delta in it is therefore >= 0, which leaves
# two things untested: whether a genuine decrease is reported with the right
# sign and magnitude, and whether the ordering is really by ABSOLUTE effect --
# with no negative values present, sorting by raw delta and sorting by |delta|
# are indistinguishable. Both matter, because the directional results this
# package was written for are predominantly decreases.
#
# Ratios are identical across samples within a site, so the mean of per-sample
# ratios equals the pooled ratio and the expected delta is exact.
directional_rows <- data.frame(
  site   = rep(c("site_down", "site_up", "site_flat"), each = 6),
  sample = rep(c(paste0("ctrl_", 1:3), paste0("diab_", 1:3)), times = 3),
  edited = c(30L, 30L, 30L, 15L, 15L, 15L,    # down: 0.60 -> 0.30  = -0.30
             10L, 10L, 10L, 20L, 20L, 20L,    # up:   0.20 -> 0.40  = +0.20
             10L, 10L, 10L, 10L, 10L, 10L),   # flat: 0.20 -> 0.20  =  0.00
  total  = rep(50L, 18)
)

make_directional_file <- function() {
  f <- tempfile(fileext = ".txt")
  write.table(directional_rows, f, sep = "\t", row.names = FALSE, quote = FALSE)
  f
}

test_that("editing difference is negative, and correct in magnitude, for a decrease", {
  f <- make_directional_file()
  on.exit(unlink(f))
  eff <- editing_difference(f, meta_path = metadata_file)

  expect_lt(eff[site == "site_down", editing_difference], 0)
  expect_equal(eff[site == "site_down", editing_difference], -0.30, tolerance = 1e-8)
  expect_equal(eff[site == "site_down", control_mean],       0.60, tolerance = 1e-8)
  expect_equal(eff[site == "site_down", diabetic_mean],      0.30, tolerance = 1e-8)
})

test_that("increases and decreases are both reported with the correct sign", {
  f <- make_directional_file()
  on.exit(unlink(f))
  eff <- editing_difference(f, meta_path = metadata_file)

  expect_equal(eff[site == "site_up",   editing_difference], +0.20, tolerance = 1e-8)
  expect_equal(eff[site == "site_flat", editing_difference],  0.00, tolerance = 1e-8)
})

test_that("ordering is by absolute effect, so a large decrease outranks a smaller increase", {
  # The decisive case the bundled fixture cannot provide: site_down (-0.30)
  # must sort ABOVE site_up (+0.20). Sorting by raw delta would invert this.
  f <- make_directional_file()
  on.exit(unlink(f))
  eff <- editing_difference(f, meta_path = metadata_file)

  expect_equal(eff$site, c("site_down", "site_up", "site_flat"))
  expect_true(all(diff(abs(eff$editing_difference)) <= 0))
})
