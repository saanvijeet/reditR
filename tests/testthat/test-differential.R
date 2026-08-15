library(testthat)
library(reditR)
library(data.table)

editing_file  <- system.file("extdata", "example_editing.txt",  package = "reditR")
metadata_file <- system.file("extdata", "example_metadata.txt", package = "reditR")

test_that("differential_editing runs all three tests by default with expected column schema", {
  res <- differential_editing(editing_file, metadata_file, verbose = FALSE)

  expect_s3_class(res, "data.table")
  expected <- c("site",
               "glmm_pvalue",   "GLMM_FDR",   "GLMM_sig",
               "fisher_pvalue", "Fisher_FDR", "Fisher_sig",
               "wilcox_pvalue", "Wilcox_FDR", "Wilcox_sig")
  expect_true(all(expected %in% names(res)))
  expect_false("DRE" %in% names(res))
})

test_that("test= restricts which columns are computed", {
  res_glmm <- differential_editing(editing_file, metadata_file,
                                   test = "glmm", verbose = FALSE)
  expect_true(all(c("glmm_pvalue", "GLMM_sig") %in% names(res_glmm)))
  expect_false(any(c("fisher_pvalue", "wilcox_pvalue") %in% names(res_glmm)))

  res_wf <- differential_editing(editing_file, metadata_file,
                                 test = c("wilcoxon", "fisher"), verbose = FALSE)
  expect_true(all(c("wilcox_pvalue", "fisher_pvalue") %in% names(res_wf)))
  expect_false("glmm_pvalue" %in% names(res_wf))
})

test_that("chr1_100 (clear effect) is GLMM significant", {
  res <- differential_editing(editing_file, metadata_file, verbose = FALSE)
  expect_true(res[site == "chr1_100", GLMM_sig])
})

test_that("chr1_200 (no effect) is not GLMM significant", {
  res <- differential_editing(editing_file, metadata_file, verbose = FALSE)
  expect_false(res[site == "chr1_200", GLMM_sig])
})

test_that("each requested test is independently FDR-corrected (no gating, no combined verdict)", {
  res <- differential_editing(editing_file, metadata_file, verbose = FALSE)

  # Wilcoxon and Fisher must be computed for every site with >=2 conditions,
  # not just GLMM-significant ones (i.e. not gated on GLMM_sig).
  expect_equal(sum(!is.na(res$wilcox_pvalue)), nrow(res))
  expect_equal(sum(!is.na(res$fisher_pvalue)), nrow(res))

  # FDR correction is independent per test.
  expect_equal(res$GLMM_FDR,   p.adjust(res$glmm_pvalue,   method = "BH"))
  expect_equal(res$Fisher_FDR, p.adjust(res$fisher_pvalue, method = "BH"))
  expect_equal(res$Wilcox_FDR, p.adjust(res$wilcox_pvalue, method = "BH"))
})

test_that("results are written to out_path when provided", {
  out <- tempfile(fileext = ".txt")
  on.exit(unlink(out))

  differential_editing(editing_file, metadata_file, out_path = out, verbose = FALSE)

  expect_true(file.exists(out))
  written <- fread(out)
  expect_equal(nrow(written), 2L)
})

test_that("summary file is written with expected rows when summary_path is provided", {
  sumf <- tempfile(fileext = ".txt")
  on.exit(unlink(sumf))

  differential_editing(editing_file, metadata_file, summary_path = sumf, verbose = FALSE)

  expect_true(file.exists(sumf))
  sumdt <- fread(sumf)
  expect_true(all(c("metric", "value") %in% names(sumdt)))
  expect_true("GLMM_sig" %in% sumdt$metric)
  expect_true("Fisher_sig" %in% sumdt$metric)
  expect_true("Wilcox_sig" %in% sumdt$metric)
  expect_false("monte_carlo_pvalue" %in% sumdt$metric)
})

test_that("random_effects supports crossed terms for pseudobulk data", {
  set.seed(42)
  libraries  <- c("lib_c1", "lib_c2", "lib_d1", "lib_d2")
  conditions <- c("control", "control", "dehydrated", "dehydrated")
  clusters   <- paste0("cluster_", 1:3)

  rows <- do.call(rbind, lapply(seq_along(libraries), function(i) {
    cond <- conditions[i]
    data.frame(
      site       = "chr1_1",
      sample     = paste0(libraries[i], "_", clusters),
      library    = libraries[i],
      cluster_id = clusters,
      condition  = cond,
      edited     = if (cond == "dehydrated") rpois(3, 8) + 2 else rpois(3, 2),
      total      = 40,
      stringsAsFactors = FALSE
    )
  }))

  ed_file <- tempfile(fileext = ".txt")
  on.exit(unlink(ed_file))
  write.table(rows, ed_file, sep = "\t", row.names = FALSE, quote = FALSE)

  res <- differential_editing(
    ed_file,
    test            = "glmm",
    random_effects  = "(1 | library) + (1 | cluster_id)",
    reference_level = "control",
    case_level      = "dehydrated",
    min_obs         = 4L,
    verbose         = FALSE
  )

  expect_s3_class(res, "data.table")
  expect_true("glmm_pvalue" %in% names(res))
  expect_false(is.na(res[site == "chr1_1", glmm_pvalue]))
})

test_that("random_effects referencing a missing column errors before any fitting", {
  ed_file <- tempfile(fileext = ".txt")
  on.exit(unlink(ed_file))
  write.table(
    data.frame(site = "chr1_1", sample = c("a", "b", "c", "d"),
               condition = c("control", "control", "dehydrated", "dehydrated"),
               edited = c(1, 2, 8, 9), total = c(20, 20, 20, 20)),
    ed_file, sep = "\t", row.names = FALSE, quote = FALSE
  )

  expect_error(
    differential_editing(ed_file, test = "glmm",
                         random_effects  = "(1 | nonexistent_col)",
                         reference_level = "control", case_level = "dehydrated",
                         verbose = FALSE),
    "random_effects references column"
  )
})

test_that("default random_effects is unchanged and still matches (1 | sample)", {
  res <- differential_editing(editing_file, metadata_file,
                              test = "glmm", verbose = FALSE)
  res_explicit <- differential_editing(editing_file, metadata_file,
                                       test = "glmm",
                                       random_effects = "(1 | sample)",
                                       verbose = FALSE)
  expect_equal(res$glmm_pvalue, res_explicit$glmm_pvalue)
})

# ---- condition contrast direction -----------------------------------------

test_that("reference_level is honoured when the case level sorts first alphabetically", {
  # glmer coerces a character `condition` to a factor with alphabetically
  # ordered levels. Without an explicit relevel, coef_name (condition<case>)
  # names a coefficient that does not exist whenever case_level sorts before
  # reference_level, the lookup error is swallowed by tryCatch, and every site
  # returns NA -- indistinguishable from universal convergence failure.
  meta  <- data.table::fread(metadata_file)
  remap <- c(control = "zref", diabetic = "acase")   # case now sorts first
  meta[, condition := as.character(remap[condition])]
  mf <- tempfile(fileext = ".txt")
  data.table::fwrite(meta, mf, sep = "\t")
  on.exit(unlink(mf))

  ref <- differential_editing(editing_file, metadata_file, test = "glmm", verbose = FALSE)
  got <- differential_editing(editing_file, mf, test = "glmm",
                              reference_level = "zref", case_level = "acase",
                              verbose = FALSE)

  expect_false(all(is.na(got$glmm_pvalue)))
  # Relabelling the same partition cannot change the Wald p-value for the
  # condition contrast -- only the sign of the estimate.
  expect_equal(got$glmm_pvalue, ref$glmm_pvalue)
  expect_equal(got$GLMM_sig,    ref$GLMM_sig)
})

test_that("all three tests are unaffected by which level sorts first", {
  meta  <- data.table::fread(metadata_file)
  remap <- c(control = "zref", diabetic = "acase")
  meta[, condition := as.character(remap[condition])]
  mf <- tempfile(fileext = ".txt")
  data.table::fwrite(meta, mf, sep = "\t")
  on.exit(unlink(mf))

  ref <- differential_editing(editing_file, metadata_file, verbose = FALSE)
  got <- differential_editing(editing_file, mf, reference_level = "zref",
                              case_level = "acase", verbose = FALSE)

  expect_equal(got$fisher_pvalue, ref$fisher_pvalue)
  expect_equal(got$wilcox_pvalue, ref$wilcox_pvalue)
})

test_that("a condition level absent from the data errors instead of returning all NA", {
  expect_error(
    differential_editing(editing_file, metadata_file, test = "glmm",
                         case_level = "not_a_condition", verbose = FALSE),
    "condition level\\(s\\) not present in data"
  )
})
