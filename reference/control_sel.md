# Control Parameters for the Selection Model

`control_sel` constructs a list with all necessary control parameters
for selection model.

## Usage

``` r
control_sel(
  est_method = c("mle", "gee"),
  gee_h_fun = 1,
  optimizer = c("maxLik", "optim"),
  maxlik_method = c("NR", "BFGS", "NM"),
  optim_method = c("BFGS", "Nelder-Mead"),
  epsilon = 1e-04,
  maxit = 500,
  trace = FALSE,
  penalty = c("SCAD", "lasso", "MCP"),
  a_SCAD = 3.7,
  a_MCP = 3,
  lambda = -1,
  lambda_min = 0.001,
  nlambda = 50,
  nfolds = 10,
  print_level = 0,
  start_type = c("zero", "mle", "naive"),
  nleqslv_method = c("Broyden", "Newton"),
  nleqslv_global = c("dbldog", "pwldog", "cline", "qline", "gline", "hook", "none"),
  nleqslv_xscalm = c("fixed", "auto"),
  dependence = FALSE,
  key = NULL
)
```

## Arguments

- est_method:

  Method of estimation for propensity score model (`"mle"` or `"gee"`;
  default is `"mle"`).

- gee_h_fun:

  Smooth function for the generalized estimating equations (GEE) method.

- optimizer:

  (for the `est_method="mle"` only) optimization function for maximum
  likelihood estimation.

- maxlik_method:

  (for the `est_method="mle"` only) maximisation method that will be
  passed to
  [`maxLik::maxLik()`](https://rdrr.io/pkg/maxLik/man/maxLik.html)
  function. Default is `NR`.

- optim_method:

  (for the `est_method="mle"` only) maximisation method that will be
  passed to [`stats::optim()`](https://rdrr.io/r/stats/optim.html)
  function. Default is `BFGS`.

- epsilon:

  Tolerance for fitting algorithms by default `1e-6`.

- maxit:

  Maximum number of iterations.

- trace:

  logical value. If `TRUE` trace steps of the fitting algorithms.
  Default is `FALSE`

- penalty:

  The penalization function used during variables selection.

- a_SCAD:

  The tuning parameter of the SCAD penalty for selection model. Default
  is 3.7.

- a_MCP:

  The tuning parameter of the MCP penalty for selection model. Default
  is 3.

- lambda:

  A user-specified \\\lambda\\ value during variable selection model
  fitting.

- lambda_min:

  The smallest value for lambda, as a fraction of `lambda.max`. Default
  is .001.

- nlambda:

  The number of `lambda` values. Default is 50.

- nfolds:

  The number of folds for cross validation. Default is 10.

- print_level:

  this argument determines the level of printing which is done during
  the optimization (for propensity score model) process.

- start_type:

  - Type of method for start points for model fitting taking the
    following values

    - if `zero` then start is a vector of zeros (default for all
      methods).

    - if `mle` (for `est_method="gee"` only) starting parameters are
      taken from the result of the `est_method="mle"` method.

- nleqslv_method:

  (for the `est_method="gee"` only) The method that will be passed to
  [`nleqslv::nleqslv()`](https://bertcarnell.github.io/nleqslv/reference/nleqslv.html)
  function.

- nleqslv_global:

  (for the `est_method="gee"` only) The global strategy that will be
  passed to
  [`nleqslv::nleqslv()`](https://bertcarnell.github.io/nleqslv/reference/nleqslv.html)
  function.

- nleqslv_xscalm:

  (for the `est_method="gee"` only) The type of x scaling that will be
  passed to
  [`nleqslv::nleqslv()`](https://bertcarnell.github.io/nleqslv/reference/nleqslv.html)
  function.

- dependence:

  logical value (default `TRUE`) informing whether samples overlap (NOT
  YET IMPLEMENTED, FOR FUTURE DEVELOPMENT).

- key:

  binary key variable allowing to identify the overlap (NOT YET
  IMPLEMENTED, FOR FUTURE DEVELOPMENT).

## Value

List with selected parameters.

## Details

Smooth function (`gee_h_fun`) for the generalized estimating equations
(GEE) method taking the following values

- if `1` then \\\boldsymbol{h}\left(\boldsymbol{x},
  \boldsymbol{\theta}\right) = \frac{\pi(\boldsymbol{x},
  \boldsymbol{\theta})}{\boldsymbol{x}}\\,

- if `2` then \\ \boldsymbol{h}\left(\boldsymbol{x},
  \boldsymbol{\theta}\right) = \boldsymbol{x}\\

## See also

[`nonprob()`](https://ncn-foreigners.github.io/nonprobsvy/reference/nonprob.md)
– for fitting procedure with non-probability samples.
