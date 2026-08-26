# Calibration outputs

The results produced by the scripts one directory up, deposited so the
published calibration numbers can be checked without re-running them. The
permutation nulls in particular took days of cluster time; re-deriving them to
verify a figure is not reasonable, so the tables are shipped instead.

These are **outputs of specific runs**, not example data. They will not
regenerate identically if the scripts are re-run with different inputs.

## Parametric arm

| File | Produced by | What it supports |
|---|---|---|
| `reditr_parametric_null.txt` | `reditr_parametric_null.R` | The 45-cell sweep: one row per (sigma, n_per_condition, seed, test), giving FPR and power at each planted effect size. False discovery proportion is derived from these as `(fpr*800)/(fpr*800 + power_all*300)`. |
| `qc_reditr_simulations.txt` | `qc_reditr_simulations.R` | Two tables. Parameter recovery across six configurations, then the dispersion comparison of simulated against real data (Pearson chi-square / df), with `pct_estimable` and `n_median` alongside. |
| `reditr_shared_effect_sim.txt` | `reditr_shared_effect_sim.R` | The refuted hypothesis: a per-sample intercept shared across sites does not reproduce the observed inflation. |
| `cross_site_correlation.txt` | `measure_cross_site_correlation.R` | Site-to-site correlation, simulated versus real, reported separately for different-chromosome pairs. |
| `realised_fdp_reditr.txt` | `reditr_parametric_null.R` | FDP with confidence bounds, per cell. |

## Permutation arm

One file per dataset. Each row is one null labelling; the columns are the
number of sites (or genes) called significant under that labelling.

| File | Design | Nulls | Floor |
|---|---|---|---|
| `perm_null_diabetes.txt` | 6 vs 6, complements collapsed | 461 | 1/462 = 0.00216 |
| `perm_null_mouse.txt` | 3 vs 10 libraries, unequal arms | 285 | 1/286 = 0.00350 |
| `perm_null_endo2_siADAR1.txt` | 3 vs 3, site level | 9 | 1/10 = 0.100 |
| `perm_null_endo2_siADAR2.txt` | 3 vs 3, site level | 9 | 0.100 |
| `perm_null_endo2_genes_siADAR1.txt` | 3 vs 3, gene level, plus p-value distribution diagnostics | 9 | 0.100 |
| `perm_null_endo2_genes_siADAR2.txt` | 3 vs 3, gene level | 9 | 0.100 |

### Scoring them

The empirical p is

```
p = (#{null >= observed} + 1) / (n_nulls + 1)
```

**The observed values are not in the diabetes and mouse files.** They are in
`observed_values.txt`. The endothelial files do carry theirs, as the row with
`null == 0`. Scoring a diabetes or mouse file against one of its own rows
yields a permuted value, not the observed one.

Checked against the deposited files:

| Dataset | Observed | Null range | p |
|---|---|---|---|
| diabetes, GLMM | 226 | 43–253 | 0.0043 |
| mouse, GLMM | 851 | 279–994 | 0.0769 |
| endothelial siADAR1, GLMM | 824 | 13–515 | 0.100 (at the floor) |

Note the mouse permuted maximum, 994, exceeds the observed 851.

## Two things to know before citing these

**`perm_null_mouse.txt` is a consolidation, not a raw output.** The mouse null
ran as a PBS array writing 159 chunk files, and array tasks that re-ran after
walltime failures produced 560 rows for 285 permutations. This file is those
rows deduplicated by permutation index, keeping the first occurrence.

Among the 275 repeated indices, Fisher agreed exactly every time, while the
**GLMM disagreed on 134 of them** — median spread 0 sites, maximum 18, out of
roughly 780. That is `glmer` convergence jitter at borderline sites. It does
not move the empirical p, but mouse GLMM site counts should not be quoted to
the last digit.

**The mouse null is not like-for-like with the other two.** Condition is
completely confounded with timepoint and partly with compartment balance, so
relabelling destroys structure the real design has.

## Not deposited

`reditr_parametric_null_raw.rds` (1.4 MB) holds the per-site results so the
sweep can be re-scored under a different metric without re-simulating. It is a
working convenience rather than a result, and a binary format that is fragile
across R versions. The numbers cited in the write-up all come from
`reditr_parametric_null.txt`.

Differential editing results for the three datasets are also not deposited
here. They are results of applying the package to particular datasets rather
than properties of the package, and belong with the data deposition.
