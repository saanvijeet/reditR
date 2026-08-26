# Calibration analyses

These scripts assess the statistical behaviour of the reditR differential-editing
workflow under controlled and data-derived null conditions.

Differential editing is tested one site at a time across many thousands of
sites, and the resulting p-values are corrected for multiple testing. That
correction is meaningful only if the underlying tests are calibrated: if
p-values are not approximately uniform under the null, the false discovery rate
is not controlled at its nominal level and the number of sites reported cannot
be interpreted. Editing datasets make this a substantive concern, because they
combine few biological replicates with high read depth, so read-level precision
can exceed the precision the number of independent samples supports.

Significance throughout is assessed at a Benjamini–Hochberg adjusted p-value
below 0.05, the default of `differential_editing()`, with each requested test
corrected independently.

## Overview

Two complementary approaches are used, because neither is sufficient alone.

**Parametric simulation** generates data with known ground truth, so
false-positive rate and power can be measured against sites planted with and
without an effect, under explicitly varied sample size, between-sample variance
and effect size. It assesses the complete `differential_editing()` workflow, but
only under the simulator's assumptions: counts are binomial conditional on a
per-sample random intercept, and sites are independent. A clean result therefore
excludes specific mechanisms rather than establishing general correctness.

**Permutation testing** reassigns condition labels within the observed
datasets, so coverage, editing-rate distributions, site-to-site structure and
sample-level variability are those of the real data and no generative model is
assumed. The question is whether the observed comparison yields more discoveries
than relabelled data produce. Resolution is bounded by the design: the smallest
attainable empirical p-value is one over the number of valid labellings plus
one, which for small designs may exceed conventional thresholds however large
the observed effect. Permutation cannot measure power, since which sites carry
an effect is unknown; simulation can, but only under its own model.

## Scripts

| Script | Dataset / design | Analysis |
|---|---|---|
| `reditr_parametric_null.R` | Simulated, known truth | Grid over variance, sample size and seed; false-positive rate and power |
| `qc_reditr_simulations.R` | Simulated vs all observed datasets | Parameter recovery; dispersion comparison |
| `reditr_shared_effect_sim.R` | Simulated | Per-site vs per-sample random-effect generation |
| `measure_cross_site_correlation.R` | Simulated vs observed | Site-to-site correlation of per-sample deviations |
| `qc_permutation_null.R` | Simulated, known truth | Controls on the permutation procedure |
| `reditr_permutation_null_diabetes.R` | Cardiomyocyte, 6 vs 6 | Site-level null, exhaustive |
| `reditr_permutation_null_endo2.R` | Endothelial knockdown, 3 vs 3 per arm | Site-level null, exhaustive |
| `reditr_permutation_null_endo2_genes.R` | Endothelial knockdown, 3 vs 3 per arm | Gene-level null; p-value diagnostics |
| `reditr_permutation_null_mouse.R` | Mouse pseudobulk, 3 vs 10 libraries | Site-level null, crossed random effects |
| `perm_null_diabetes_array.pbs` | — | PBS array driver, cardiomyocyte |
| `perm_null_mouse_array.pbs` | — | PBS array driver, mouse |

## Parametric simulation

`reditr_parametric_null.R` builds a factorial grid over between-sample
random-effect standard deviation on the logit scale, samples per condition and
replicate seed. Each cell simulates a dataset with `simulate_editing_data()`
holding true-null sites plus equal numbers carrying three known effect sizes on
the probability scale, then runs `differential_editing()` for the GLMM, Fisher's
exact test and the Wilcoxon rank-sum test. Results are joined to the truth
table, and false-positive rate, power overall and at each effect size, and the
proportion of sites returning a p-value — for the GLMM, the convergence rate —
are computed per test.

Null and effect sites share a dataset so both quantities come from the same run:
comparing tests on power alone is uninformative when their false-positive rates
differ. Sample size is swept because the tests do not respond to it in the same
direction.

`qc_reditr_simulations.R` verifies that the simulator produces what it was asked
for, comparing recovered sample and site counts, baseline editing rate,
coverage, logit-scale variance and realised effect sizes against the requested
values. It then compares dispersion — Pearson chi-square over its degrees of
freedom, per site and condition — between simulated and observed data, since the
parametric arm is informative only if the simulation is not an easier problem
than reality. As dispersion is undefined for a site and condition seen in one
sample, the proportion of estimable groups is reported alongside it.

### Shared sample-effect simulation

`simulate_editing_data()` draws its random intercept inside the per-site loop,
so a sample running high at one site has no tendency to do so at the next. Real
data need not behave this way: a library with a systematic offset would deviate
coherently across many sites at once. `reditr_shared_effect_sim.R` isolates this
distinction, changing exactly one thing — whether the random effect is drawn
once per site or once per sample and reused across sites — while holding
baseline rates, coverage, sample sizes and effect sizes constant. Both modes are
generated at every setting so the comparison is like-for-like, and a coherent
per-sample offset statistic confirms the manipulation produced the intended
structure.

`measure_cross_site_correlation.R` measures the correlation between sites of
their per-sample deviation vectors. Because sites retained by
`filter_editing_sites()` may be genomically clustered and share reads,
correlations are reported separately for site pairs on different chromosomes,
which cannot.

## Permutation analyses

Condition labels are reassigned with group sizes preserved, and the complete
`differential_editing()` call is refitted for each labelling using identical
filters, formula, optimiser and FDR threshold. Only the labels change. The
empirical p-value is

```
p = (#{null >= observed} + 1) / (n_nulls + 1)
```

so the smallest attainable value is `1 / (n_nulls + 1)`, a property of the
design that cannot be lowered by additional computation.

Three conventions apply: the observed labelling is excluded from the null; for
balanced two-group comparisons a labelling and its complement describe the same
partition and are counted once; and the permutation unit is the biological
replicate, which is not always the row of the analysis table. Permuting rows
instead would create labellings the experiment could not have produced, and
would understate the null.

`qc_permutation_null.R` validates this procedure rather than any dataset, on
simulated data with known truth. Under a negative control with no planted effect
the observed labelling is one draw from the same null as the rest, so its
empirical p-value should be uniform; under a positive control it should rank
first. Together these show whether a null result reflects absence of signal
rather than absence of power. Re-run it if the permutation code is adapted.

## Dataset-specific permutation designs

### Cardiomyocyte (diabetes)

Twelve samples in a balanced two-group design, permuted at sample level. As the
arms are equal, complementary assignments are collapsed; excluding the observed
labelling leaves `C(12,6)/2 - 1 = 461` valid nulls, few enough to enumerate
exhaustively. All three tests are evaluated, with `(1 | sample)` as the
random-effects term.

### Endothelial knockdown (endothelial2)

Two arms, siADAR1 and siADAR2, each compared against the same scrambled-control
samples with three per group. They are analysed separately because they are
distinct experimental contrasts rather than levels of one factor. With equal arms
and complements collapsed, `C(6,3)/2 - 1 = 9` valid nulls remain, so enumeration
is exhaustive by construction. Analysis is at site level with `(1 | sample)`; a
companion script performs the gene-level equivalent.

### Mouse dehydration

Single-cell derived pseudobulk data, in which each analysis unit is a
library–cluster combination but the biological replicate is the sequencing
library. Labels are therefore permuted at library level, so every pseudobulk unit
inherits its library's label and no library is split between arms. The arms are
unequal, so complements are distinct labellings and are not collapsed; excluding
the observed labelling leaves `C(13,3) - 1 = 285` valid nulls. The model uses
crossed random effects, `(1 | library) + (1 | cluster_id)`, since cluster labels
recur across libraries.

An important limitation applies here: condition is confounded with dehydration
timepoint and partly with compartment balance, so relabelling destroys structure
the real design carries. This null is consequently not directly comparable with
the other two.

## Gene-level aggregation

`reditr_permutation_null_endo2_genes.R` repeats the endothelial permutations,
recording results at both site and gene level. Sites are mapped to annotated
gene symbols, and a gene counts as significant if at least one of its sites is —
the same rule used in the corresponding observed analysis, so permuted and
observed counts follow identical logic. A gene-level null must be generated
explicitly rather than derived from site counts, because the mapping depends on
how significant sites distribute across genes.

It also records p-value distribution diagnostics: the median p-value and a
Kolmogorov–Smirnov test against a uniform distribution. Counting significant
sites summarises the distribution at a single threshold, whereas uniformity under
the null is what FDR control asserts.

## Reproducibility

These scripts are records of specific analyses rather than portable utilities.
Input and output locations are set near the top of each and must be changed
before use elsewhere; the permutation scripts also assume a particular group
structure and sample naming.

Where supported, behaviour is controlled by environment variables:

| Variable | Used by | Controls |
|---|---|---|
| `SIGMAS`, `NPC`, `N_REPS`, `N_NULL` | `reditr_parametric_null.R` | Simulation grid and dataset size |
| `N_CORES` | Most scripts | Cores used for per-site fitting |
| `PERM_FROM`, `PERM_TO` | Cardiomyocyte and mouse permutation scripts | Contiguous range of permutations to run |
| `NMAX` | `measure_cross_site_correlation.R` | Sites sampled per dataset |
| `N_REP`, `N_SITES`, `SIGMA` | `qc_permutation_null.R` | Control-experiment size and variance |

The two larger permutation analyses run as PBS array jobs, with
`perm_null_diabetes_array.pbs` and `perm_null_mouse_array.pbs` assigning each
task a contiguous range. Chunked runs write suffixed output so concurrent tasks
cannot overwrite one another; these must be concatenated and deduplicated by
permutation index before scoring.

Result tables are provided under `outputs/`, with a separate README describing
their contents.

## Relationship to the thesis

These scripts provide the computational implementation of the calibration and
validation analyses described in the thesis Methods. The corresponding findings
are reported and interpreted in the Results and Discussion, and are not
reproduced here.
