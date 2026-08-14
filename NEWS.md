# reditR 0.2.2

* `filter_editing_sites()` gains a `min_edit_ratio` argument (default `0`,
  unchanged behaviour) applying a minimum per-observation editing ratio
  after the coverage and edited-read filters. Applied only when
  `min_edit_ratio > 0`: `edit_ratio` is `NaN` where `total == 0` and
  `NaN >= 0` is `NA`, so an unconditional filter would silently discard
  zero-coverage rows whenever `min_coverage` is `0`.

* `filter_editing_sites()` documentation now states that only
  quality-control filters are supported, and that design-dependent
  filters -- requiring a site in both arms, for example -- belong with
  the caller, since they depend on the contrast rather than on data
  quality and do not generalise across two-arm, multi-arm and paired
  designs. The two-line idiom is given in `@details`.

* `inst/scripts/bulk_bam_to_fastq.sh` now handles paired-end BAMs. It
  previously passed only `-fq` to `bamToFastq`, which writes both mates
  into one file with no mate information, so a paired-end BAM would be
  silently flattened and then processed by SPRINT as single-end. Layout
  is detected from the first alignment record's FLAG (overridable with
  `PAIRED=0`/`PAIRED=1`), and paired input is split into `<out>_1.fastq`
  and `<out>_2.fastq`. The template is not covered by the test suite,
  which has no BAM fixtures; this is now stated in the script header.

* `inst/scripts/scrna_preprocessing.sh` no longer describes SPRINT's `-c`
  as a "minimum coverage cutoff". `-c` trims the first N bp of each read.
  The comment is corrected and now warns against confusing it with the
  cluster-size defaults `-csrg`/`-cshp`, which happen to share the value 5.

* New `inst/scripts/test_data/mouse/`: the two annotation builders that
  produced the analysed mouse (E-MTAB-8145) clusters, copied verbatim,
  with a README recording that they apply no minimum cluster-size filter
  -- unlike `inst/scripts/build_splitter_annotation.sh`, which is a
  single-library smoke test and not the production path.

# reditR 0.2.1

* `differential_editing()` gains a `random_effects` argument controlling
  the random-effects side of the GLMM formula (default `"(1 | sample)"`,
  unchanged behaviour). Lets pseudobulk callers pass crossed or nested
  terms, e.g. `"(1 | library) + (1 | cluster_id)"`, using whatever
  grouping columns are present in the data -- a single `sample`-level
  term pseudoreplicates when multiple pseudobulk units share a library
  and/or cluster identity. Errors up front (before fitting anything) if
  `random_effects` references a column not present in the data, rather
  than silently returning `NA` for every site.

# reditR 0.2.0

* **Breaking change:** `differential_editing()` redesigned around three
  independent, user-selectable tests rather than a single GLMM-primary /
  Wilcoxon-robustness-flag framework:
  - New `test` argument (any subset of `c("glmm", "fisher", "wilcoxon")`,
    default all three). Only the requested test(s) are computed.
  - All requested tests now run on **every** site, independently of one
    another -- Wilcoxon is no longer gated on GLMM significance.
  - Each requested test gets its own independent Benjamini-Hochberg FDR
    correction (`GLMM_FDR`/`GLMM_sig`, `Fisher_FDR`/`Fisher_sig`,
    `Wilcox_FDR`/`Wilcox_sig`).
  - Fisher's exact test is now a first-class option (previously absent).
  - Removed: the global Monte Carlo permutation overlap test (and its
    `n_perm`/`seed` arguments and `monte_carlo_pvalue` attribute), the
    combined `DRE` verdict column, and the raw-p-value Wilcoxon
    "robustness flag" convention. Requesting more than one test lets you
    compare columns directly, or build your own combined criterion.
  - Guidance on which test to pick for a given experimental design (based
    on null-calibration simulations across three independent datasets) is
    now documented in `?differential_editing` and the package vignette.
* `validate_against_truth()`/`print.reditR_validation()` updated to report
  false-positive rate and power independently per test present in the
  `differential_editing()` output, instead of hard-coding GLMM + DRE.
* Removed the `permutation_pvalue()` export: it was declared in `NAMESPACE`
  and documented, but never actually implemented in this package -- calling
  it would have errored. No functional change (it was already unusable).

# reditR 0.1.0

* Initial release.
* Core functions: `read_editing_table()`, `read_metadata()`,
  `filter_editing_sites()`, `differential_editing()`,
  `editing_difference()`, `simulate_editing_data()`,
  `validate_against_truth()`.
* Implements the GLMM + per-site robustness-flag framework for
  small-sample differential editing analysis.
