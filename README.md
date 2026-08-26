# reditR

> Differential RNA editing analysis for bulk and single-cell pseudobulk data.

reditR provides statistical tools for detecting differentially RNA-edited (DRE)
sites from SPRINT-derived read-count tables. It is the statistical downstream
layer of the SPRINT pipeline. 

## Installation

```r
# install.packages("devtools")
devtools::install_github("saanvijeet/reditR")
```

## Quick example

```r
library(reditR)

# 1. Read data
ed <- read_editing_table("filtered_sites_clustered.txt")
mt <- read_metadata("sample_metadata.txt", reference_level = "control")

# 2. Filter (or use filter_editing_sites() for SPRINT output)
# filter_editing_sites("all_samples_editing.txt", out_dir = "results/")

# 3. Run differential analysis -- pick the test(s) that suit your design
res <- differential_editing(
  data_path = "filtered_sites_clustered.txt",
  meta_path = "sample_metadata.txt",
  test      = c("glmm", "fisher", "wilcoxon"),
  out_path  = "DRE_results.txt"
)

# 4. Effect sizes - computes the mean editing ratio for each site and condition
eff <- editing_difference("filtered_sites_clustered.txt", meta_path = "sample_metadata.txt")
head(eff)
```

## Pipeline integration

The SPRINT output parsing script is shipped at:

```r
system.file("scripts", "extracting_read_counts.sh", package = "reditR")
```

Run it on your HPC to produce `all_samples_editing.txt`, then call
`filter_editing_sites()` from within R.

## Methodology

reditR offers **three independent significance tests**, each run on every
site and independently Benjamini-Hochberg FDR corrected. Pick
whichever test(s) suit your experimental design with the `test` argument
(any subset of `c("glmm", "fisher", "wilcoxon")`), and compare the
`<Test>_sig`/`<Test>_FDR` columns you asked for directly.

| Test | Model | Notes |
|---|---|---|
| **GLMM** | `cbind(edited, unedited) ~ condition + (1\|sample)` | Primary test (Srivastava et al. 2017). Can be severely anti-conservative when between-replicate variance is high relative to the number of independent control replicates. |
| **Fisher** | Exact test on pooled per-sample counts | Independent of GLMM. Anti-conservative whenever real between-replicate variance is present, regardless of dataset -- and unlike GLMM, adding replicates does not fix it. |
| **Wilcoxon** | Rank-sum on per-sample editing ratios | The most consistently well-calibrated test across every design checked, at the cost of lower power and a coarse p-value floor at small sample sizes. |

**Which test should I use?** This isn't a rule of thumb — it comes from
running permutation and parametric null-calibration simulations across three
independent datasets (see the package vignette for methodology). A design
with few, unbalanced, or confounded replicates should lean on Wilcoxon; a
well-powered, balanced design can reasonably use GLMM; Fisher should be
treated cautiously whenever real between-replicate variance is plausible.
Requesting more than one test lets you compare them directly instead of
trusting either in isolation:

```r
res[GLMM_sig == TRUE & Wilcox_sig == TRUE]   # your own AND, if you want one
```

## Simulation and power analysis

```r
sim <- simulate_editing_data(n_null = 800, n_effects = c("0.10" = 100, "0.20" = 100))
res <- differential_editing(..., case_level = "case")   # run on sim$editing / sim$metadata
val <- validate_against_truth(res, sim$truth)
print(val)   # false-positive rate and power, reported per test requested
```

## Citation

TBD
