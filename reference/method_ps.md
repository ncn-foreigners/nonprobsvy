# Propensity Score Model Functions

Function to specify the propensity score (PS) model for the inverse
probability weighting estimator. This function provides basic functions
logistic regression with a given link function (currently we support
`logit`, `probit` and `cloglog`) with additional information about the
analytic variance estimator of the mean.

This is a function returns a `list` of functions that refer to specific
estimation methods and variance estimators when whether the IPW alone or
the DR estimator is applied. The export of this function is mainly
because the functions are used in the variable selection algorithms.

Functions starting with `make_log_like`, `make_gradient` and
`make_hessian` refer to the maximum likelihood estimation as described
in the Chen et al. (2020) paper. These functions take into account
different link functions defined through the `link` argument.

Functions `make_link`, `make_link_inv`, `make_link_der`,
`make_link_inv_der`, `make_link_inv_rev`, and `make_link_inv_rev_der`
refer to specific link functions and are used in the estimation process.

Functions `variance_covariance1` and `variance_covariance2` refer to the
variance estimator of the IPW estimator as defined by Chen et al.
(2020).

Functions `b_vec_ipw`, `b_vec_dr` and `t_vec` are specific functions
defined in the Chen et al. (2020) that are used in the variance
estimator of the IPW or the DR.

Finally, `var_nonprob` is the non-probability component of the DR
estimator as defined by Chen et al. (2020).

## Usage

``` r
method_ps(link = c("logit", "probit", "cloglog"), ...)
```

## Arguments

- link:

  link for the PS model

- ...:

  Additional, optional arguments.

## Value

A `list` of functions and elements for a specific link function with the
following entries:

- make_log_like:

  log-likelihood function for a specific link function

- make_gradient:

  gradient of the loglik

- make_hessian:

  hessian of the loglik

- make_link:

  link function

- make_link_inv:

  inverse link function

- make_link_der:

  first derivative of the link function

- make_link_inv_der:

  first derivative of the the inverse link function

- make_link_inv_rev:

  defines 1/inv_link

- make_link_inv_rev_der:

  first derivative of 1/inv_link

- variance_covariance1:

  for the IPW estimator: variance component for the non-probability
  sample

- variance_covariance2:

  for the IPW estimator: variance component for the probability sample

- b_vec_ipw:

  for the IPW estimator: the \\b\\ function as defined in the Chen et
  al. (2020, sec. 3.2, eq. (9)-(10); sec 4.1)

- b_vec_dr:

  for the DR estimator: the \\b\\ function as defined in the Chen et al.
  (2020, sec. 3.3., eq. (14); sec 4.1)

- t_vec:

  for the DR estimator: the \\b\\ function as defined in the Chen et al.
  (2020, sec. 3.3., eq. (14); sec 4.1)

- var_nonprob:

  for the DR estimator: non-probability component of the variance for DR
  estimator

- link:

  name of the selected link function for the PS model (character)

- model:

  model type (character)

## Author

Łukasz Chrostowski, Maciej Beręsewicz

## Examples

``` r
# Printing information on the model selected

method_ps()
#> Propensity score model with logit link

# extracting specific field

method_ps("cloglog")$make_gradient
#> function (X_nons, X_rand, weights, weights_rand, ...) 
#> {
#>     function(theta) {
#>         eta1 <- as.matrix(X_nons) %*% theta
#>         eta2 <- as.matrix(X_rand) %*% theta
#>         invLink1 <- inv_link(eta1)
#>         invLink2 <- inv_link(eta2)
#>         exp_eta1 <- stable_cloglog_rate(eta1)
#>         exp_eta2 <- stable_cloglog_rate(eta2)
#>         t(crossprod(X_nons, weights * exp_eta1/invLink1) - crossprod(X_rand, 
#>             weights_rand * exp_eta2))
#>     }
#> }
#> <bytecode: 0x55e22e35a440>
#> <environment: 0x55e229be58d0>
```
