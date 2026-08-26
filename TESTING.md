# reditR Testing Log

## Overview

This document records the full testing workflow for the `reditR` package: the tools used, each error encountered during installation and testing, and the fix applied.

---

## Test Suite

Tests are written using **testthat 3rd edition** and live in `tests/testthat/`. There are 4 test files covering all major functions:

| File | Functions tested | Blocks | Expectations |
|---|---|---|---|
| `test-io.R` | `read_editing_table`, `read_metadata` | 6 | 11 |
| `test-filter.R` | `filter_editing_sites` | 9 | 15 |
| `test-differential.R` | `differential_editing` | 13 | 33 |
| `test-effect_size.R` | `editing_difference` | 9 | 16 |
| **Total** | | **37** | **75** (0 failures) |

Two of the seven exported functions, `simulate_editing_data()` and
`validate_against_truth()`, have no tests. They are exercised by the
calibration scripts under `inst/scripts/calibration/` rather than by the suite.

(`permutation_pvalue` was removed from `test-differential.R` during the 0.2.0
rewrite: it was declared as an export and tested, but never implemented
anywhere in the package, so that test would have errored against a clean
install.)

---

## Step 1: Installing devtools, roxygen2, testthat

### Problem

`devtools`, `roxygen2`, and `testthat` were not present in the `rna_editing` conda environment. Running the standard install command:

```r
install.packages(c("devtools", "roxygen2", "testthat"),
                 repos = "https://cloud.r-project.org", Ncpus = 4)
```

failed for 18 packages because the conda environment was missing system C libraries required to compile from source:

- **`systemfonts`** → `cannot find -lfreetype` (libfreetype not installed)
- **`httpuv`** → `fatal error: zlib.h: No such file or directory` (zlib headers not installed)

These are system-level dependencies that `install.packages()` cannot resolve on its own when compiling from source.

### Fix

Used conda-forge instead, which distributes **pre-compiled binaries** that bundle all system dependencies:

```bash
conda install -n rna_editing -c conda-forge r-devtools r-roxygen2 r-testthat -y
```

After install, `devtools` failed to load because `usethis` and then `pkgload`, `pkgbuild`, `callr`, and several other dependencies were also missing. These were installed in a follow-up:

```bash
conda install -n rna_editing -c conda-forge r-usethis r-pkgload r-pkgbuild \
  r-remotes r-sessioninfo r-desc r-rprojroot r-ellipsis r-withr r-fs \
  r-pkgdown r-callr r-processx -y
```

After this, `library(devtools)` loaded successfully.

---

## Step 2: Generating Documentation with roxygen2

### Command

```r
roxygen2::roxygenize('/rds/general/user/sj1825/home/msc_prj/reditR')
```

### Problem

roxygen2 emitted two warnings:

```
✖ reditR-package.R:7: @importFrom must be only 1 line long, not 2.
✖ reditR-package.R:10: @importFrom must be only 1 line long, not 2.
```

The `@importFrom` tags in `reditR-package.R` were wrapped across two lines:

```r
#' @importFrom data.table ":=" ".N" as.data.table data.table fread fwrite
#'   is.data.table rbindlist setnames uniqueN
#' @importFrom stats binomial p.adjust plogis qlogis rbinom rnbinom rnorm
#'   wilcox.test
```

### Fix

Collapsed each `@importFrom` tag onto a single line:

```r
#' @importFrom data.table ":=" ".N" as.data.table data.table fread fwrite is.data.table rbindlist setnames uniqueN
#' @importFrom stats binomial p.adjust plogis qlogis rbinom rnbinom rnorm wilcox.test
```

After this, `roxygenize()` completed cleanly and generated all 10 Rd man pages.

---

## Step 3: Installing the Package

### Command

```r
devtools::install('/rds/general/user/sj1825/home/msc_prj/reditR')
```

Completed with no errors:

```
* DONE (reditR)
```

---

## Step 4: Running Tests

### Command

```r
devtools::test('/rds/general/user/sj1825/home/msc_prj/reditR')
```

### Initial Result

```
[ FAIL 6 | WARN 0 | SKIP 0 | PASS 34 ]
```

All 6 failures were in `test-effect_size.R`:

```
ERROR: 'test-effect_size.R:9:3'
Error: object 'edit_ratio' not found
Backtrace:
  └─reditR::editing_difference(editing_file, metadata_file)
    └─...[] at reditR/R/effect_size.R:61:3
```

### Root Cause

`editing_difference()` reads the input file with `fread()` and then immediately uses `edit_ratio` as a bare column name in a `data.table` aggregation:

```r
editing_summary <- data[, .(
  case_mean = mean(edit_ratio[condition == case_level], na.rm = TRUE),
  ...
), by = site]
```

The example data file (`inst/extdata/example_editing.txt`) only has columns `site`, `sample`, `edited`, and `total`; it does not include a pre-computed `edit_ratio` column. The function never computed it before trying to use it.

### Fix

Added a guard to compute `edit_ratio` from `edited / total` if it is absent, matching the behaviour of `read_editing_table()`:

```r
# Compute edit_ratio if absent
if (!"edit_ratio" %in% names(data))
  data[, edit_ratio := edited / total]
```

This was inserted in `R/effect_size.R` immediately after `fread()` and before the merge with metadata.

### Result After Fix

```
[ FAIL 0 | WARN 0 | SKIP 0 | PASS 42 ]
```

All 42 tests pass.

---

## Step 5: R CMD Check

### Command

```r
devtools::check('/rds/general/user/sj1825/home/msc_prj/reditR', cran = FALSE)
```

### Initial Result

```
0 errors | 1 warning | 1 note
```

### Warning: Non-ASCII Characters

```
checking code files for non-ASCII characters ... WARNING
Found the following files with non-ASCII characters:
  R/differential.R
  R/simulate.R
```

The `∩` (intersection) symbol had been used in two places to represent the logical AND of GLMM and Wilcoxon results:

- `R/differential.R` line 126, in a comment
- `R/differential.R` line 155, inside a `cat()` call
- `R/simulate.R` line 179, inside a `cat()` / `sprintf()` call

**Fix:** Replaced `∩` with `&` in all three locations.

---

### Note: Undefined Global Functions

```
checking R code for possible problems ... NOTE
editing_difference: no visible global function definition for '.'
filter_editing_sites: no visible global function definition for 'read.table'
filter_editing_sites: no visible global function definition for 'write.table'
read_editing_table: no visible global function definition for '.'
read_metadata: no visible global function definition for '.'
Undefined global functions or variables:
  . read.table write.table
```

Two separate sub-issues:

**`.` (data.table list alias)**
The `.()` shorthand used in `data.table` aggregations is not a formally exported symbol from the `data.table` package, so `R CMD check` cannot resolve it through the normal `@importFrom` mechanism. Attempting to add `@importFrom data.table "."` caused roxygen2 to emit `Excluding unknown export from data.table: "."`.

**Fix:** Added `"."` to the `utils::globalVariables()` declaration in `R/reditR-package.R`.

**`read.table` / `write.table`**
`filter_editing_sites()` (copied from `filtered_sites.R`) calls these base R utility functions. They are not in `stats` and need to be explicitly imported from the `utils` package.

**Fix:** Added the following `@importFrom` tag to `R/reditR-package.R`:

```r
#' @importFrom utils read.table write.table
```

---

### Final Result

```
0 errors | 0 warnings | 0 notes
```

---

## Final Package State (0.1.0)

| Check | Result |
|---|---|
| `devtools::document()` | OK, 10 Rd files generated |
| `devtools::install()` | OK, `* DONE (reditR)` |
| `devtools::test()` | 42/42 pass |
| `devtools::check()` | 0 errors, 0 warnings, 0 notes |

---

## Step 6: 0.2.0 rewrite, three independent tests, no Monte Carlo

`differential_editing()` was redesigned from a GLMM-primary /
Wilcoxon-robustness-flag / Monte-Carlo-overlap framework to three
independent, user-selectable tests (`test = c("glmm","fisher","wilcoxon")`,
each run on every site, each independently BH-corrected, no combined `DRE`
verdict). Full rationale in `NEWS.md`. Package also moved out of the
`msc_prj` monorepo into its own standalone repository at this point.

### Discovered in the process

**`fisher.test` was an undeclared dependency.** Fisher's exact test is new
in 0.2.0 and calls `fisher.test()` unqualified from within the package
namespace, which requires an explicit import (base/recommended packages are
not automatically visible inside a package namespace the way they are in a
user session). Fixed by adding `fisher.test` to the `@importFrom stats` tag
in `R/reditR-package.R` and re-running `roxygen2::roxygenise()`.

**`validate_against_truth()`/`print.reditR_validation()` hard-coded the old
`DRE` column.** Since `DRE` no longer exists, these degraded silently
(reporting `NaN` for the DRE row) rather than erroring. Rewritten to detect
whichever `<Test>_sig` columns are actually present in the `differential_editing()`
output and report false-positive rate / power per test found, generically.

**`permutation_pvalue()` was already broken pre-0.2.0.** It was declared in
`NAMESPACE`/`man/permutation_pvalue.Rd`/`tests/testthat/test-differential.R`,
but never defined anywhere in `R/`. Calling it would have errored on a clean
install. Re-running `roxygen2::roxygenise()` correctly dropped the export
and orphaned man page once nothing in `R/` backed it. Not reimplemented as
part of this rewrite -- flagged as a separate, pre-existing gap.

### Result after fix

```
devtools::test(".")
[ FAIL 0 | WARN 0 | SKIP 0 | PASS 47 ]
```

Also smoke-tested against a real (non-toy) dataset with all three tests
requested and `n_cores = 4`, and the results matched an independent
from-scratch script exactly. The counts themselves are reported in the
dissertation rather than here.

## Package state at the 0.2.0 rewrite

| Check | Result |
|---|---|
| `devtools::document()` | OK, Rd files regenerated, `permutation_pvalue.Rd` correctly dropped |
| `devtools::test()` | 47/47 pass |
| Real-data smoke test | Matched an independent from-scratch script exactly |

## Changes since

This document records the setup and verification work carried out up to the
0.2.0 rewrite. Three releases have followed, each described in `NEWS.md`:

| Version | Change |
|---|---|
| 0.2.1 | `random_effects` argument, allowing crossed or nested terms for pseudobulk designs |
| 0.2.2 | Explicit condition contrast, `min_edit_ratio`, paired-end handling in `bulk_bam_to_fastq.sh`, and added test coverage for the negative direction in `editing_difference()` |
| 0.2.3 | `condition_labels` argument for `simulate_editing_data()` |

The test-suite table at the top of this document reflects the current suite,
not the 47 expectations recorded above for 0.2.0.
