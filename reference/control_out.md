# Control Parameters for Outcome Model

`control_out` constructs a list with all necessary control parameters
for outcome model.

## Usage

``` r
control_out(
  epsilon = 1e-08,
  maxit = 100,
  trace = FALSE,
  k = 5,
  penalty = c("SCAD", "lasso", "MCP"),
  a_SCAD = 3.7,
  a_MCP = 3,
  lambda_min = 0.001,
  nlambda = 100,
  nfolds = 10,
  treetype = c("kd", "rp", "ball"),
  searchtype = c("standard", "priority"),
  pmm_match_type = 1,
  pmm_weights = c("none", "dist"),
  pmm_k_choice = c("none", "min_var"),
  pmm_reg_engine = c("glm", "loess"),
  npar_loess = stats::loess.control(surface = "direct", trace.hat = "approximate")
)
```

## Arguments

- epsilon:

  Tolerance for fitting algorithms. Default is `1e-6`.

- maxit:

  Maximum number of iterations.

- trace:

  logical value. If `TRUE` trace steps of the fitting algorithms.
  Default is `FALSE`.

- k:

  The k parameter in the
  [`RANN::nn2()`](https://jefferislab.github.io/RANN/reference/nn2.html)
  function. Default is 5.

- penalty:

  penalty algorithm for variable selection. Default is `SCAD`

- a_SCAD:

  The tuning parameter of the SCAD penalty for outcome model. Default is
  3.7.

- a_MCP:

  The tuning parameter of the MCP penalty for outcome model. Default is
  3.

- lambda_min:

  The smallest value for lambda, as a fraction of lambda.max. Default is
  .001.

- nlambda:

  The number of lambda values. Default is 100.

- nfolds:

  The number of folds during cross-validation for variables selection
  model.

- treetype:

  Type of tree for nearest neighbour imputation (for the NN and PMM
  estimator) passed to
  [`RANN::nn2()`](https://jefferislab.github.io/RANN/reference/nn2.html)
  function.

- searchtype:

  Type of search for nearest neighbour imputation (for the NN and PMM
  estimator) passed to
  [`RANN::nn2()`](https://jefferislab.github.io/RANN/reference/nn2.html)
  function.

- pmm_match_type:

  (Only for the PMM Estimator) Indicates how to select 'closest' unit
  from non-probability sample for each unit in probability sample.
  Either `1` (default) or `2` where `2` is matching by minimizing
  distance between \\\hat{y}\_{i}\\ for \\i \in S\_{A}\\ and \\y\_{j}\\
  for \\j \in S\_{B}\\ and `1` is matching by minimizing distance
  between \\\hat{y}\_{i}\\ for \\i \in S\_{A}\\ and \\\hat{y}\_{i}\\ for
  \\i \in S\_{A}\\.

- pmm_weights:

  (Only for the PMM Estimator) Indicate how to weight `k` nearest
  neighbours in \\S\_{B}\\ to create imputed value for units in
  \\S\_{A}\\. The default value `"none"` indicates that mean of `k`
  nearest \\y\\'s from \\S\_{B}\\ should be used whereas `"prop_dist"`
  results in weighted mean of these `k` values where weights are
  inversely proportional to distance between matched values.

- pmm_k_choice:

  (Only for the PMM Estimator) Character value indicating how `k`
  hyper-parameter should be chosen, by default `"none"` meaning `k`
  provided in `control_outcome` argument will be used. For now the only
  other option `"min_var"` means that `k` will be chosen by minimizing
  estimated variance of estimator for mean. Parameter `k` provided in
  this control list will be chosen as starting point.

- pmm_reg_engine:

  (Only for the PMM Estimator) whether to use parametric (`"glm"`) or
  non-parametric (`"loess"`) regression model for the outcome. The
  default is `"glm"`.

- npar_loess:

  control parameters for the
  [stats::loess](https://rdrr.io/r/stats/loess.html) via the
  [stats::loess.control](https://rdrr.io/r/stats/loess.control.html)
  function.

## Value

List with selected parameters.

## See also

[`nonprob()`](https://ncn-foreigners.github.io/nonprobsvy/reference/nonprob.md)
– for fitting procedure with non-probability samples.
