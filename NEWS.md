# reditR 0.2.3

* `simulate_editing_data()` gains a `condition_labels` argument, replacing the
  hardcoded `"control"` / `"diabetic"` arm names. The default is now
  `c("control", "case")`. The labels are arbitrary names for the two arms and
  do not enter the generative model, so simulated counts are unchanged for a
  given seed -- only the `condition` column of the returned metadata differs.
  The previous names were residue from the package having been written against
  the cardiomyocyte dataset first, and made every calibration script pass
  `case_level = "diabetic"` to simulations that had nothing to do with
  diabetes.

  **Note the resulting asymmetry:** `differential_editing()` and
  `editing_difference()` still default to `case_level = "diabetic"`, so
  handing simulated data straight to either with default arguments now errors.
  It errors up front with a clear message rather than silently returning all
  `NA` (see 0.2.2), but the case level must be supplied. Changing those
  defaults is a breaking change and is deferred.

  Callers updated: `inst/scripts/calibration/reditr_shared_effect_sim.R` and
  `qc_permutation_null.R`, the vignette, and the README example.
  `reditr_parametric_null.R` needed no change -- it already derives the case
  level from the simulated metadata rather than assuming it.

# reditR 0.2.2

* Test coverage for the negative direction in `editing_difference()`. The
  bundled fixture contains only an increase (chr1_100, +0.17) and a flat site
  (chr1_200, -0.001), so every delta in it is >= 0. Two things were therefore
  untested: whether a genuine decrease is reported with the correct sign and
  magnitude, and whether the output ordering is really by ABSOLUTE effect --
  with no negative values present, `order(-abs(delta))` and `order(-delta)`
  produce identical output, so the existing sort assertion could not
  distinguish them. Three blocks added, using a self-contained temporary
  fixture spanning -0.30, +0.20 and 0.00; the ordering assertion was verified
  to discriminate between the two sort keys. No change to any function -- the
  behaviour was already correct, but it was not pinned.

* New `inst/scripts/calibration/`: the scripts used to assess whether
  `differential_editing()` controls its false discovery rate: the parametric
  sweep over between-sample variance and sample size, QC on the simulator, the
  shared-effect and cross-site-correlation experiments, and the four
  permutation nulls (one per dataset, since the labelling space depends on the
  design) with their PBS array wrappers. Previously these existed only in the
  analysis tree, so the calibration results could not be reproduced from the
  package. They are records of specific analyses rather than portable
  utilities: paths are hardcoded and the permutation scripts assume a
  particular group structure. `qc_permutation_null.R` validates the
  permutation procedure itself and should be re-run if that code is adapted.

* **Bug fix:** `differential_editing()` now sets the `condition` contrast
  explicitly with `reference_level` as the base. Previously `condition`
  reached `glmer()` as a character vector, which R coerces to a factor with
  alphabetically ordered levels -- not necessarily `reference_level` first.
  The Wald p-value is looked up by the coefficient name
  `paste0("condition", case_level)`, so whenever `case_level` sorted before
  `reference_level` the fitted coefficient carried the other name, the lookup
  failed, the error was swallowed by the per-site `tryCatch()`, and **every
  site silently returned `NA`** -- a wholly null result indistinguishable
  from universal convergence failure. All condition pairs used in the
  bundled examples happen to sort safely (`control` < `diabetic`), so
  existing results are unaffected, but e.g.
  `reference_level = "treated", case_level = "control"` was silently broken.

* `differential_editing()` now errors up front if `reference_level` or
  `case_level` is absent from the data, instead of failing later with
  `object 'glmm_pvalue' not found` once every site has been skipped.

* `filter_editing_sites()` gains a `min_edit_ratio` argument (default `0`,
  unchanged behaviour) applying a minimum per-observation editing ratio
  after the coverage and edited-read filters. Applied only when
  `min_edit_ratio > 0`: `edit_ratio` is `NaN` where `total == 0` and
  `NaN >= 0` is `NA`, so an unconditional filter would silently discard
  zero-coverage rows whenever `min_coverage` is `0`.

* `filter_editing_sites()` documentation now states that only
  quality-control filters are supported, and that design-dependent
  filters, such as requiring a site in both arms, belong with the caller,
  since they depend on the contrast rather than on data quality and do not
  generalise across two-arm, multi-arm and paired designs.

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
