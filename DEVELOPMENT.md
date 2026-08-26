# reditR Developer Notes

Design decisions and worked examples for changes that have been considered but
not implemented. Each entry records enough detail that the change can be made
without re-deriving it, and states why it was left out.

---

## Adding clinical covariates to the GLMM

**Status:** not implemented.

The fixed-effects side of the GLMM formula is `condition` alone
(`R/differential.R`, in the `if ("glmm" %in% test)` block):

```r
glmm_formula <- stats::as.formula(
  paste("cbind(edited, unedited) ~ condition +", random_effects)
)
```

`random_effects` (added in 0.2.1) opened up the random side, so a caller can
pass crossed or nested grouping terms for pseudobulk data. There is no
equivalent for the fixed side, so a clinical or technical variable (age, sex,
HbA1c, batch) cannot currently be adjusted for.

### What the change would look like

Four edits, of which only the first two touch the model.

**1. Take the covariates as an argument.** Default `NULL`, so every existing
call behaves identically:

```r
differential_editing <- function(data_path,
                                 ...,
                                 random_effects = "(1 | sample)",
                                 covariates     = NULL,
                                 ...)
```

**2. Paste them into the fixed side, next to `condition`:**

```r
fixed <- paste(c("condition", covariates), collapse = " + ")
glmm_formula <- stats::as.formula(
  paste("cbind(edited, unedited) ~", fixed, "+", random_effects)
)
```

With `covariates = c("age", "sex")` this builds

```
cbind(edited, unedited) ~ condition + age + sex + (1 | sample)
```

and with `covariates = NULL` it collapses to the current formula exactly, so
the default path is unchanged.

**3. Check the columns exist before fitting anything.** Without this a typo
returns `NA` at every site. This is the failure mode of the `reference_level`
bug present in 0.2.1 and fixed in 0.2.2, where the error was swallowed by the
per-site `tryCatch()` and the result was indistinguishable from universal
convergence failure. The existing
random-effects check already covers covariate names once they are in the
formula; it only needs generalising from "random effects" to "formula":

```r
fixed_vars <- setdiff(all.vars(glmm_formula),
                      c("edited", "unedited", "condition"))
missing    <- setdiff(fixed_vars, names(data))
if (length(missing) > 0)
  stop("formula references column(s) not found in data: ",
       paste(missing, collapse = ", "))
```

**4. Nothing in `run_glmm()` changes.** It fits whatever formula it is given,
and the p-value is looked up by coefficient name:

```r
coef_name <- paste0("condition", case_level)
summary(m)$coefficients[coef_name, "Pr(>|z|)"]
```

Additional fixed terms shift the condition coefficient's row position without
affecting a name-based lookup. The reported p-value stays the condition
effect, now adjusted for the covariates.

### Requirements on the covariates

- **Must be a column of the data** passed to `differential_editing()`, not of
  the metadata file. A clinical variable is a property of the donor, so it has
  to be merged in alongside `condition`; the existing merge on `sample` will
  carry it if it is present in `meta_path`.
- **Must be constant within a sample.** A per-row covariate would be modelling
  something other than a donor characteristic.
- **Must not be collinear with `condition`.** A variable that only takes one
  value in each arm (a treatment given only to cases, say) cannot be separated
  from the condition effect and produces a rank-deficient fit.
- **Categorical covariates should be factors.** A character column is coerced
  with alphabetical level ordering, the same trap that made the explicit
  `condition` contrast necessary in 0.2.2.

### Why it is not implemented

Degrees of freedom. The model is fitted **per site**, so each covariate costs
a degree of freedom at every site, and the datasets this package was developed
on carry 6 to 13 samples. A per-site model already estimating an intercept, a
condition effect and a variance component has very little room for more.

The null-calibration work is relevant here. The Wald statistic is already
unreliable at some sites in these datasets: a small number returned p-values
between 1e-119 and true double-precision underflow, where a likelihood ratio
test on the same fitted models returned p between 0.59 and 1.00. Those sites
were *not* the singular ones, so the failure is in the Wald approximation to a
non-quadratic likelihood surface rather than in degenerate variance
estimation. Adding fixed terms to a model already at the edge of
identifiability would be expected to make that worse, not better.

So the argument against is not that the change is difficult, since it is four
edits and two of them are trivial, but that on a design of this size the
adjusted model is unlikely to be estimable site by site. The change is worth
making for callers with larger sample sizes; it should ship with a note in
`?differential_editing` recording the above, and ideally with the calibration
sweep re-run on a design carrying a covariate before it is recommended.

### If implemented, also

- Add a `covariates` entry to the roxygen block and re-run `roxygenise()`.
- Add a test asserting `covariates = NULL` reproduces the current formula
  byte-for-byte, and one asserting a missing covariate column errors up front
  rather than returning all-`NA`.
- Add a `NEWS.md` entry under the next version.
