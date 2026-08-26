# Null-calibration scripts

The scripts used to assess whether `differential_editing()` controls its false
discovery rate. They are deposited so the published calibration results can be
reproduced and adapted.

**These are records of specific analyses, not portable utilities.** Every one
hardcodes absolute paths and, for the permutation scripts, a particular
experimental design. Read the header of any script before running it.

## Two arms, testing different things

| Arm | Data | Question | Weakness |
|---|---|---|---|
| Parametric | simulated, known truth | How does the false-positive rate respond as between-sample variance grows? | Circular — validates the package against the model it assumes |
| Permutation | real, condition labels shuffled | Does the observed result separate from a null built from the actual data? | Only probes the variance that dataset carries; the achievable p has a hard floor |

The parametric arm *excludes mechanisms*; the permutation arm *detects the
failure*. Where they disagree, the disagreement is the finding.

## Parametric

| Script | Purpose |
|---|---|
| `reditr_parametric_null.R` | The sweep. `sample_re_sd` × `n_per_condition` × 3 seeds = 45 cells, each 800 true-null plus 300 known-effect sites, via `simulate_editing_data()` → `differential_editing()` → `validate_against_truth()`. Grid is set by environment variables (`SIGMAS`, `NPC`, `N_REPS`, `N_NULL`, `N_CORES`), so this is the most reusable script here. |
| `qc_reditr_simulations.R` | QC on the generator: parameter recovery, and dispersion (Pearson chi-square / df) of simulated against all three real datasets. Reports `pct_estimable` beside each dispersion, because dispersion does not exist for a site × condition group of one observation and small designs are mostly such groups. |
| `reditr_shared_effect_sim.R` | Tests whether a per-sample random intercept shared *across* sites — which `simulate_editing_data()` does not generate — reproduces the observed inflation. It does not. |
| `measure_cross_site_correlation.R` | Measures site-to-site correlation in simulated versus real data, reported separately for different-chromosome pairs so read-sharing between nearby sites is excluded. |

## Permutation

One per dataset, because the labelling space depends on the design. All refit
the complete `differential_editing()` call with identical filters, formula,
optimiser and FDR threshold; only the condition labels change.

The empirical p is `(#{null >= observed} + 1) / (n_nulls + 1)`, so **the floor
is `1/(n_nulls + 1)` and is a property of the design, not of the analysis.**
Getting this wrong is easy and consequential: two of these datasets were
originally run with 10 sampled labellings, pinning both at p = 0.091, and
exhaustive enumeration later moved one of them to p = 0.0043.

| Script | Design | Nulls | Floor |
|---|---|---|---|
| `reditr_permutation_null_mouse.R` | 3 vs 10 libraries, unequal arms so no complement symmetry | 285 = C(13,3) − 1 | 0.0035 |
| `reditr_permutation_null_diabetes.R` | 6 vs 6, complements collapse | 461 = C(12,6)/2 − 1 | 0.00216 |
| `reditr_permutation_null_endo2.R` | 3 vs 3, site-level | 9 = C(6,3)/2 − 1 | 0.100 |
| `reditr_permutation_null_endo2_genes.R` | 3 vs 3, gene-level aggregation plus p-value-distribution diagnostics | 9 | 0.100 |

`perm_null_mouse_array.pbs` and `perm_null_diabetes_array.pbs` run the two
large ones as PBS arrays; both scripts accept `PERM_FROM` / `PERM_TO` to select
a contiguous chunk, and write chunk-suffixed output so concurrent tasks cannot
overwrite each other. Concatenate and **deduplicate by permutation index**
before scoring.

See `outputs/` for the results these scripts produced, with a README mapping
each file to what it supports.

`qc_permutation_null.R` validates the permutation *procedure* rather than any
dataset, using positive and negative controls on simulated data with known
truth. Run it if you adapt the permutation code: it is what establishes that a
null result means "no signal" rather than "the test has no power".

## Adapting these to a new dataset

`reditr_parametric_null.R` needs only its output path changed. The permutation
scripts need more: the labelling enumeration assumes a specific group structure,
and the random-effects specification is dataset-dependent — the mouse uses
`(1 | library) + (1 | cluster_id)` rather than the package default
`(1 | sample)`, because a pseudobulk unit is unique per row and `(1 | sample)`
is near-degenerate there.

## Things worth knowing before citing the outputs

**Observed values are not stored in the diabetes and mouse permutation output
files.** They were only printed to console. Scoring against a row of one of
those files yields a permuted value, not the observed one — the observed counts
are in `outputs/observed_values.txt`. The endothelial files do carry theirs, as
the row with `null == 0`.

**The GLMM is reproducible to roughly ±0.5% on the mouse data.** Where
permutations were inadvertently repeated, Fisher agreed exactly on all of them
while the GLMM differed by a median of 4 significant sites out of ~780. This is
`glmer` convergence jitter at borderline sites. It does not move the empirical
p, but exact site counts should not be quoted to the last digit.

**The mouse permutation null is not like-for-like with the other two.**
Condition is completely confounded with timepoint and partly with compartment
balance, so relabelling destroys structure the real design has.
