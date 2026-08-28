# reditR

> Differential RNA editing analysis for bulk and single-cell pseudobulk data.

reditR provides statistical tools for detecting differentially RNA-edited (DRE)
sites from per-site, per-sample read-count tables, such as those derived from
SPRINT output. It is the statistical downstream layer of an editing-detection
pipeline: the package API begins at the count table, and site detection itself
is performed externally.

## Requirements

R >= 4.0.0.

| Type | Packages |
|---|---|
| Imports | `data.table`, `dplyr`, `lme4`, `parallel`, `stats`, `tidyr` |
| Suggests | `testthat` (>= 3.0.0), `knitr`, `rmarkdown` |

Installing with `devtools::install_github()` will pull the imports
automatically. **SPRINT is not an R package and is not installed by reditR.** It
is an external tool that must be installed separately if you intend to run the
detection stage; reditR does not call it.

## Installation

```r
# install.packages("devtools")
devtools::install_github("saanvijeet/reditR")
```

## Input formats

Two tab-separated files are required. Both are read with
`data.table::fread()`, so any delimiter it auto-detects will work, but
tab-separated is what the shipped pipeline scripts produce.

### Editing count table

One row per site per sample. Required columns:

| Column | Type | Description |
|---|---|---|
| `site` | character | Site identifier, unique per genomic position |
| `sample` | character | Sample identifier, matching the metadata file |
| `edited` | integer | Number of reads supporting the edited base |
| `total` | integer | Total reads covering the site in that sample |

```
site        sample    edited    total
chr1_100    ctrl_1    2         30
chr1_100    ctrl_2    1         25
chr1_100    diab_1    8         35
chr1_200    ctrl_1    1         28
```

An `edit_ratio` column may be present but is not required: it is computed as
`edited / total` where absent. Additional columns are permitted and are used
where the analysis needs them, for example grouping columns named in
`random_effects`. Note that `read_editing_table()` returns only
`site`, `sample`, `edited`, `total` and `edit_ratio`, so if you need extra
columns downstream, pass the file path directly to the analysis functions
rather than routing it through the reader.

### Sample metadata

One row per sample. Required columns:

| Column | Type | Description |
|---|---|---|
| `sample` | character | Sample identifier, matching the count table |
| `condition` | character | Experimental group label |

```
sample    condition
ctrl_1    control
ctrl_2    control
diab_1    diabetic
diab_2    diabetic
```

Only the two conditions named by `reference_level` and `case_level` are
analysed; other rows are dropped, so a multi-arm experiment is run as a series
of pairwise comparisons. Both defaults must be checked against your own labels:
`reference_level` defaults to `"control"` and `case_level` to `"diabetic"`. If
either is absent from the data, the analysis stops with an error naming the
levels that are available.

`reference_level` sets the baseline of the comparison. The reported effect is
case relative to reference, so reversing the two reverses the sign. If the
metadata file already contains a `condition` column merged into the count table,
`meta_path` may be omitted.

## Functions

**Input**

| Function | Description |
|---|---|
| `read_editing_table()` | Read and validate a count table, computing `edit_ratio` if absent |
| `read_metadata()` | Read and validate a metadata file, optionally setting the reference level as the first factor level |

**Filtering**

| Function | Description |
|---|---|
| `filter_editing_sites()` | Apply per-observation coverage, edited-read and editing-ratio thresholds, a minimum-sample requirement, and identify genomically clustered sites |

**Testing**

| Function | Description |
|---|---|
| `differential_editing()` | Run any subset of three significance tests per site, each independently FDR-corrected |

**Effect sizes**

| Function | Description |
|---|---|
| `editing_difference()` | Per-site mean editing ratio in each condition and their signed difference, ordered by absolute effect |

**Simulation and validation**

| Function | Description |
|---|---|
| `simulate_editing_data()` | Generate synthetic count data with known planted effects |
| `validate_against_truth()` | Score analysis output against simulated truth, reporting convergence, false-positive rate and power per test |

## Minimal reproducible example

The package bundles a small example dataset so the workflow can be run without
your own data. It contains 2 sites across 6 samples, 3 control and 3 diabetic,
and uses the default condition labels.

```r
library(reditR)

ed <- system.file("extdata", "example_editing.txt",  package = "reditR")
mt <- system.file("extdata", "example_metadata.txt", package = "reditR")

# Inspect the inputs
read_editing_table(ed)
read_metadata(mt, reference_level = "control")

# Differential testing
res <- differential_editing(ed, mt,
                            test    = c("glmm", "fisher", "wilcoxon"),
                            verbose = FALSE)
res

# Effect sizes
editing_difference(ed, meta_path = mt)
```

`res` contains one row per site with, for each requested test, a raw p-value, a
BH-adjusted value and a logical call: `glmm_pvalue` / `GLMM_FDR` / `GLMM_sig`,
and correspondingly `fisher_pvalue` / `Fisher_FDR` / `Fisher_sig` and
`wilcox_pvalue` / `Wilcox_FDR` / `Wilcox_sig`. No combined verdict column is
produced; intersect the columns yourself if you want one.

`editing_difference()` returns `site`, a mean column per condition named from
the labels (here `control_mean` and `diabetic_mean`), and `editing_difference`,
the case minus reference difference.

The example is small enough that the GLMM produces a singular fit warning. This
is expected at this size and does not indicate an error.

## Simulation

```r
sim <- simulate_editing_data(n_null = 40, n_effects = c("0.20" = 10),
                             n_per_condition = 4, seed = 1)

d <- tempfile(); m <- tempfile()
data.table::fwrite(sim$editing,  d, sep = "\t")
data.table::fwrite(sim$metadata, m, sep = "\t")

res <- differential_editing(d, m, test = c("glmm", "fisher"),
                            case_level = "case", verbose = FALSE)
validate_against_truth(res, sim$truth)
```

`simulate_editing_data()` labels its arms `control` and `case` by default, set
by `condition_labels`, so `case_level = "case"` must be passed explicitly since
`differential_editing()` defaults to `"diabetic"`.

## Methodology

reditR offers three independent significance tests, each run on every site and
independently Benjamini-Hochberg FDR corrected. Select them with the `test`
argument, which accepts any subset of `c("glmm", "fisher", "wilcoxon")`.

| Test | Model | What it uses |
|---|---|---|
| **GLMM** | `cbind(edited, unedited) ~ condition + <random_effects>` | Binomial mixed model on read counts. The fixed effect is `condition`; the random-effects term is supplied by the caller through the `random_effects` argument. The unit of replication is whatever grouping variable that term names. Fitted with `lme4::glmer()`. |
| **Fisher** | Exact test on a 2x2 table of pooled counts | Edited and unedited reads summed across all samples within each condition. Sample identity is not represented, so the unit of replication is the read. |
| **Wilcoxon** | Rank-sum on per-sample editing ratios | Ranks `edited / total` per sample and compares the two condition groups. The unit of replication is the sample; read depth does not enter beyond forming the ratio. |

The `random_effects` argument takes a character string giving everything after
`condition +` in the formula. The default is `"(1 | sample)"`, appropriate for
bulk data where each sample is one biological replicate:

```r
# Bulk: one random intercept per sample
differential_editing(ed, mt, random_effects = "(1 | sample)")

# Pseudobulk: crossed terms, where each row is a library-cluster unit
differential_editing(ed, mt,
                     random_effects = "(1 | library) + (1 | cluster_id)")
```

Any grouping column named in `random_effects` must be present in the count
table. The function checks this before fitting and stops with an error naming
any missing column.

`min_obs` sets how many observations a site needs before the GLMM is fitted to
it at all. The default is 4. A site with fewer rows than this is skipped and
returns no GLMM p-value, because there is too little data to estimate an
intercept, a condition effect and a variance component from. Fisher and
Wilcoxon are unaffected by this setting, so the three tests can end up scored
on different numbers of sites. If you report how many sites the GLMM called
significant, give that count against the number of sites it actually fitted
rather than the number tested.

## Pipeline integration

The steps upstream of the count table are shell and Python scripts, run outside
R. They are shipped with the package for reference and are located with
`system.file("scripts", "<name>", package = "reditR")`.

| Script | Language | Role |
|---|---|---|
| `build_splitter_annotation.sh` | shell | Build the barcode-to-cluster annotation table required by the BAM splitter, from CellRanger output |
| `bam_extract_barcode_reads_commandline_chr_V2.py` | Python | Split a multiplexed single-cell BAM into per-cell-type BAMs using that annotation |
| `scrna_preprocessing.sh` | shell | Single-cell driver: splits the BAM, converts each per-cell-type BAM to FASTQ, then runs SPRINT per cell type |
| `bulk_bam_to_fastq.sh` | shell | Convert a bulk BAM to FASTQ for SPRINT input. Only needed when starting from BAM rather than FASTQ |
| `extracting_read_counts.sh` | shell | Parse `SPRINT_identified_regular.res` files into the combined `all_samples_editing.txt` count table |

`inst/scripts/test_data/mouse/` additionally holds the two R scripts that built
the barcode-to-cluster annotations for the mouse dataset analysed in the
dissertation. They are records of that specific analysis rather than portable
utilities.

`inst/scripts/calibration/` holds the scripts and outputs for the calibration
and validation analyses, documented in their own README.

The R workflow starts once `all_samples_editing.txt` exists:

```r
filter_editing_sites("all_samples_editing.txt", out_dir = "results/")
```

## Testing

```r
devtools::test()
```

The suite comprises 37 test blocks and 75 expectations across four files,
covering the readers, the filtering thresholds and their boundary behaviour,
the three significance tests and their independent FDR correction, the
random-effects specification, condition-level handling, and effect-size sign
and ordering. `R CMD check` runs it automatically.

## Citation

Saanvi Jiteendra, MSc Dissertation, Genomic Medicine, Imperial College London.
