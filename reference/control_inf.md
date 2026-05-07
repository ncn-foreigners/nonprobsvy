# Control Parameters for Inference

`control_inf` constructs a `list` with all necessary control parameters
for statistical inference.

## Usage

``` r
control_inf(
  var_method = c("analytic", "bootstrap"),
  rep_type = c("subbootstrap", "auto", "JK1", "JKn", "BRR", "bootstrap", "mrbbootstrap",
    "Fay"),
  vars_selection = FALSE,
  vars_combine = FALSE,
  bias_correction = FALSE,
  num_boot = 500,
  alpha = 0.05,
  cores = 1,
  keep_boot = TRUE,
  nn_exact_se = FALSE
)
```

## Arguments

- var_method:

  the variance method (default `"analytic"`).

- rep_type:

  the replication type for weights in the bootstrap method for variance
  estimation passed to
  [`survey::as.svrepdesign()`](https://rdrr.io/pkg/survey/man/as.svrepdesign.html).
  Default is `"subbootstrap"`.

- vars_selection:

  default `FALSE`; if `TRUE`, then the variables selection model is
  used.

- vars_combine:

  whether variables should be combined after variable selection for
  doubly robust estimators (default `FALSE`)

- bias_correction:

  default `FALSE`; if `TRUE`, then the bias minimization estimation used
  during model fitting.

- num_boot:

  the number of iteration for bootstrap algorithms.

- alpha:

  significance level (default 0.05).

- cores:

  the number of cores in parallel computing (default 1).

- keep_boot:

  a logical value indicating whether statistics from bootstrap should be
  kept (default `TRUE`)

- nn_exact_se:

  a logical value indicating whether to compute the exact standard error
  estimate for `nn` or `pmm` estimator. The variance estimator for
  estimation based on `nn` or `pmm` can be decomposed into three parts,
  with the third computed using covariance between imputed values for
  units in the probability sample using predictive matches from the
  non-probability sample. In most situations this term is negligible and
  is very computationally expensive so by default it is set to `FALSE`,
  but the recommended option is to set this value to `TRUE` before
  submitting the final results.

## Value

A `list` with selected parameters.

## See also

[`nonprob()`](https://ncn-foreigners.github.io/nonprobsvy/reference/nonprob.md)
– for fitting procedure with non-probability samples.
