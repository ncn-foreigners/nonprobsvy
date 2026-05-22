# Mass Imputation Using the Generalized Linear Model Method

Model for the outcome for the mass imputation estimator using
generalized linear models via the
[`stats::glm`](https://rdrr.io/r/stats/glm.html) function. Estimation of
the mean is done using \\S\_{\mathrm{P}}\\ probability sample or known
population totals.

## Usage

``` r
method_glm(
  y_nons,
  X_nons,
  X_rand,
  svydesign,
  weights = NULL,
  family_outcome = "gaussian",
  start_outcome = NULL,
  vars_selection = FALSE,
  pop_totals = NULL,
  pop_size = NULL,
  control_outcome = control_out(),
  control_inference = control_inf(),
  verbose = FALSE,
  se = TRUE
)
```

## Arguments

- y_nons:

  target variable from non-probability sample

- X_nons:

  a `model.matrix` with auxiliary variables from non-probability sample

- X_rand:

  a `model.matrix` with auxiliary variables from non-probability sample

- svydesign:

  a svydesign object

- weights:

  case / frequency weights from non-probability sample

- family_outcome:

  family for the glm model

- start_outcome:

  start parameters (default `NULL`)

- vars_selection:

  whether variable selection should be conducted

- pop_totals:

  population totals from the `nonprob` function

- pop_size:

  population size from the `nonprob` function

- control_outcome:

  controls passed by the `control_out` function

- control_inference:

  controls passed by the `control_inf` function (currently not used, for
  further development)

- verbose:

  parameter passed from the main `nonprob` function

- se:

  whether standard errors should be calculated

## Value

an `nonprob_method` class which is a `list` with the following entries

- model_fitted:

  fitted model either a `glm.fit`, `cv.ncvreg`, or `ncvreg` object

- y_nons_pred:

  predicted values for the non-probablity sample

- y_rand_pred:

  predicted values for the probability sample or population totals

- coefficients:

  coefficients for the model (if available)

- svydesign:

  an updated `surveydesign2` object (new column `y_hat_MI` is added)

- y_mi_hat:

  estimated population mean for the target variable

- vars_selection:

  whether variable selection was performed

- var_prob:

  variance for the probability sample component (if available)

- var_nonprob:

  variance for the non-probability sampl component

- var_total:

  total variance, if possible it should be `var_prob+var_nonprob` if
  not, just a scalar

- model:

  model type (character `"glm"`)

- family:

  family type (character `"glm"`)

## Details

Analytical variance

The variance of the mean is estimated based on the following approach

\(a\) non-probability part (\\S\_{\mathrm{NP}}\\ with size
\\n\_{\mathrm{NP}}\\; denoted as `var_nonprob` in the result)

\$\$ \hat{V}\_1 =
\frac{1}{n\_{\mathrm{NP}}^2}\sum\_{i=1}^{n\_{\mathrm{NP}}} \hat{e}\_i^2
\left\lbrace \boldsymbol{h}(\boldsymbol{x}\_i;
\hat{\boldsymbol{\beta}})^\prime\hat{\boldsymbol{c}}\right\rbrace^2,
\$\$

where \\\hat{e}\_i = y_i - m(\boldsymbol{x}\_i;
\hat{\boldsymbol{\beta}})\\ and
\$\$\widehat{\boldsymbol{c}}=\left\lbrace n\_{\mathrm{NP}}^{-1} \sum\_{i
\in S\_{\mathrm{NP}}} \dot{\boldsymbol{m}}\left(\boldsymbol{x}\_i ;
\boldsymbol{\beta}^\*\right) \boldsymbol{h}\left(\boldsymbol{x}\_i ;
\boldsymbol{\beta}^\*\right)^{\prime}\right\rbrace^{-1} N^{-1} \sum\_{i
\in S\_{\mathrm{P}}} w_i \dot{\boldsymbol{m}}\left(\boldsymbol{x}\_i ;
\boldsymbol{\beta}^\*\right).\$\$

Under the linear regression model
\\\boldsymbol{h}\left(\boldsymbol{x}\_i ;
\widehat{\boldsymbol{\beta}}\right)=\boldsymbol{x}\_i\\ and
\\\widehat{\boldsymbol{c}}=\left(n\_{\mathrm{NP}}^{-1} \sum\_{i \in
S\_{\mathrm{NP}}} \boldsymbol{x}\_i
\boldsymbol{x}\_i^{\prime}\right)^{-1} N^{-1} \sum\_{i \in
S\_{\mathrm{P}}} w_i \boldsymbol{x}\_i .\\

\(b\) probability part (\\S\_{\mathrm{P}}\\ with size
\\n\_{\mathrm{P}}\\; denoted as `var_prob` in the result)

This part uses functionalities of the `{survey}` package and the
variance is estimated using the following equation:

\$\$ \hat{V}\_2=\frac{1}{N^2} \sum\_{i=1}^{n\_{\mathrm{P}}}
\sum\_{j=1}^{n\_{\mathrm{P}}} \frac{\pi\_{i j}-\pi_i \pi_j}{\pi\_{i j}}
\frac{m(\boldsymbol{x}\_i; \hat{\boldsymbol{\beta}})}{\pi_i}
\frac{m(\boldsymbol{x}\_i; \hat{\boldsymbol{\beta}})}{\pi_j}. \$\$

Note that \\\hat{V}\_2\\ in principle can be estimated in various ways
depending on the type of the design and whether population size is known
or not.

Furthermore, if only population totals/means are known and assumed to be
fixed we set \\\hat{V}\_2=0\\.

Information on the case when `svydesign` is not available:

1.  variance is estimated only for the non-probability part with
    \\\hat{V}\_1\\ defined above.

2.  point estimator of \\\hat{\mu}\_y\\ for linear regression is
    estimated using \\\mu_x^\prime\hat{\boldsymbol{\beta}}\\ where
    \\\mu_x\\ is the vector of population means

3.  for non-linear functions such as logistic or Poisson regression we
    use a simplification, i.e. we report point estimate as
    \\\exp(\mu_x^\prime\hat{\boldsymbol{\beta}})\\ for Poisson and
    \\\frac{\exp(\mu_x^\prime\hat{\boldsymbol{\beta}})}{1+\exp(\mu_x^\prime\hat{\boldsymbol{\beta}})}\\
    for logistic regression.

## References

Kim, J. K., Park, S., Chen, Y., & Wu, C. (2021). Combining
non-probability and probability survey samples through mass imputation.
Journal of the Royal Statistical Society Series A: Statistics in
Society, 184(3), 941-963.

## Examples

``` r

data(admin)
data(jvs)
jvs_svy <- svydesign(ids = ~ 1,  weights = ~ weight, strata = ~ size + nace + region, data = jvs)

res_glm <- method_glm(y_nons = admin$single_shift,
                      X_nons = model.matrix(~ region + private + nace + size, admin),
                      X_rand = model.matrix(~ region + private + nace + size, jvs),
                      svydesign = jvs_svy)

res_glm
#> Mass imputation model (GLM approach). Estimated mean: 0.7039 (se: 0.0115)
```
