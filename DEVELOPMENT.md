# reditR Developer Notes

Implementation guides for extensions to the package. 

---

## Adding clinical covariates to the GLMM

The fixed-effects side of the GLMM formula is `condition` alone
(`R/differential.R`, in the `if ("glmm" %in% test)` block):

```r
glmm_formula <- stats::as.formula(
  paste("cbind(edited, unedited) ~ condition +", random_effects)
)
```

`random_effects` (added in 0.2.1) opens up the random side, so a caller can
pass crossed or nested grouping terms for pseudobulk data. The fixed side has
no equivalent, so adjusting for a clinical or technical variable (age, sex,
HbA1c, batch) means the four edits below. Only the first two touch the model.

### 1. Take the covariates as an argument

Default `NULL`, so every existing call behaves identically:

```r
differential_editing <- function(data_path,
                                 ...,
                                 random_effects = "(1 | sample)",
                                 covariates     = NULL,
                                 ...)
```

### 2. Paste them into the fixed side, next to `condition`

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

### 3. Check the columns exist before fitting anything

Without this a typo returns `NA` at every site. That is the failure mode of the
`reference_level` bug present in 0.2.1 and fixed in 0.2.2, where the error was
swallowed by the per-site `tryCatch()` and the result was indistinguishable
from universal convergence failure. The existing random-effects check already
covers covariate names once they are in the formula; it only needs generalising
from "random effects" to "formula":

```r
fixed_vars <- setdiff(all.vars(glmm_formula),
                      c("edited", "unedited", "condition"))
missing    <- setdiff(fixed_vars, names(data))
if (length(missing) > 0)
  stop("formula references column(s) not found in data: ",
       paste(missing, collapse = ", "))
```

### 4. Leave `run_glmm()` alone

It fits whatever formula it is given, and the p-value is looked up by
coefficient name:

```r
coef_name <- paste0("condition", case_level)
summary(m)$coefficients[coef_name, "Pr(>|z|)"]
```

Additional fixed terms shift the condition coefficient's row position without
affecting a name-based lookup. The reported p-value stays the condition effect,
now adjusted for the covariates.

### Requirements on the covariates

- **Must be a column of the data** passed to `differential_editing()`, not of
  the metadata file. A clinical variable is a property of the donor, so it has
  to be merged in alongside `condition`; the existing merge on `sample` will
  carry it if it is present in `meta_path`.
- **Must be constant within a sample.** A per-row covariate would be modelling
  something other than a donor characteristic.
- **Must not be collinear with `condition`.** A variable taking one value in
  each arm (a treatment given only to cases, say) cannot be separated from the
  condition effect and produces a rank-deficient fit.
- **Categorical covariates should be factors.** A character column is coerced
  with alphabetical level ordering, the same trap that made the explicit
  `condition` contrast necessary in 0.2.2.

### Sizing the model

A separate model is fitted at every site, and each one has only your samples to
work with. With six samples, the model is already estimating three things from
six numbers: a baseline editing rate, the difference between the two
conditions, and how much samples vary among themselves. Every covariate you add
is one more thing estimated from those same six numbers.

Eventually there is not enough data to go round. When that happens `glmer`
either fails to fit, and the site returns `NA`, or it returns something that
looks like a p-value but should not be trusted. The second case is the
dangerous one, because nothing in the output tells you it happened.

The datasets this package was built on have 6 to 13 samples, which is small
enough that even one covariate may not fit. Larger designs have more room, and
that is where this argument is worth using.

Whatever the size, do not judge it by eye. Run the calibration sweep in
`inst/scripts/calibration/` on a simulated design that matches yours and
carries a covariate. That measures whether the p-values behave on a design of
that size, which you cannot tell from the results themselves.

### Also do

- Add a `covariates` entry to the roxygen block and re-run `roxygenise()`.
- Record the sizing caveat above in `?differential_editing`.
- Add a test asserting `covariates = NULL` reproduces the current formula
  byte-for-byte, and one asserting a missing covariate column errors up front
  rather than returning all-`NA`.
