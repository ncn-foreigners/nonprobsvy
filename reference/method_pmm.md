# Mass Imputation Using Predictive Mean Matching Method

Model for the outcome for the mass imputation estimator. The
implementation is currently based on
[RANN::nn2](https://jefferislab.github.io/RANN/reference/nn2.html)
function and thus it uses Euclidean distance for matching units from
\\S_A\\ (non-probability) to \\S_B\\ (probability) based on predicted
values from model \\\boldsymbol{x}\_i\\ based either on `method_glm` or
`method_npar`. Estimation of the mean is done using \\S_B\\ sample.
Matching ties are randomized by the nearest-neighbour step before donor
values are aggregated.

This implementation extends Yang et al. (2021) approach as described in
Chlebicki et al. (2025), namely:

- pmm_weights:

  if k\>1 weighted aggregation of the mean for a given unit is used. We
  use distance matrix returned by
  [RANN::nn2](https://jefferislab.github.io/RANN/reference/nn2.html)
  function (`pmm_weights` from the
  [`control_out()`](https://ncn-foreigners.github.io/nonprobsvy/reference/control_out.md)
  function)

- nn_exact_se:

  if the non-probability sample is small we recommend using a
  mini-bootstrap approach to estimate variance from the non-probability
  sample (`nn_exact_se` from the
  [`control_inf()`](https://ncn-foreigners.github.io/nonprobsvy/reference/control_inf.md)
  function). If non-constant pseudo-weights are supplied, bootstrap
  samples are drawn with probabilities proportional to inverse weights
  and the resampled weights are used in each refitted outcome model.

- pmm_k_choice:

  the main `nonprob` function allows for dynamic selection of `k`
  neighbours based on the variance minimization procedure
  (`pmm_k_choice` from the
  [`control_out()`](https://ncn-foreigners.github.io/nonprobsvy/reference/control_out.md)
  function)

## Usage

``` r
method_pmm(
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

  case / frequency weights from non-probability sample. If
  `nn_exact_se=TRUE`, non-constant weights also define mini-bootstrap
  sampling probabilities proportional to their inverses and are
  resampled for each bootstrap refit.

- family_outcome:

  family for the glm model

- start_outcome:

  start parameters

- vars_selection:

  whether variable selection should be conducted

- pop_totals:

  a place holder (not used in `method_pmm`)

- pop_size:

  population size from the `nonprob` function

- control_outcome:

  controls passed by the `control_out` function

- control_inference:

  controls passed by the `control_inf` function

- verbose:

  parameter passed from the main `nonprob` function

- se:

  whether standard errors should be calculated

## Value

an `nonprob_method` class which is a `list` with the following entries

- model_fitted:

  fitted model either an `glm.fit` or `cv.ncvreg` object

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

- model:

  model type (character `"pmm"`)

- family:

  depends on the method selected for estimating E(Y\|X)

## Details

Matching

In the package we support two types of matching:

1.  \\\hat{y} - \hat{y}\\ matching (default;
    `control_out(pmm_match_type = 1)`).

2.  \\\hat{y} - y\\ matching (`control_out(pmm_match_type = 2)`).

Analytical variance

The variance of the mean is estimated based on the following approach
(a) non-probability part (\\S_A\\ with size \\n_A\\; denoted as
`var_nonprob` in the result) is currently estimated using the
non-parametric mini-bootstrap estimator proposed by Chlebicki et al.
(2025, Algorithm 2). It is not proved to be consistent but with good
finite population properties. This bootstrap can be applied using
`control_inference(nn_exact_se=TRUE)` and can be summarized as follows:

1.  Sample \\n_A\\ units from \\S_A\\ with replacement to create
    \\S_A'\\. If non-constant pseudo-weights are supplied through
    `weights`, sampling probabilities are proportional to their
    inverses; equal weights use uniform resampling.

2.  Estimate regression model
    \\\mathbb{E}\[Y\|\boldsymbol{X}\]=m(\boldsymbol{X}, \cdot)\\ based
    on \\S\_{A}'\\ from step 1, using the resampled weights when
    supplied.

3.  Compute \\\hat{\nu}'(i,t)\\ for \\t=1,\dots,k, i\in S\_{B}\\ using
    estimated \\m(\boldsymbol{x}', \cdot)\\ and
    \\\left\lbrace(y\_{j},\boldsymbol{x}\_{j})\| j\in
    S\_{A}'\right\rbrace\\.

4.  Compute
    \\\displaystyle\frac{1}{k}\sum\_{t=1}^{k}y\_{\hat{\nu}'(i)}\\ using
    \\Y\\ values from \\S\_{A}'\\.

5.  Repeat steps 1-4 \\M\\ times (we set (hard-coded) \\M=50\\ in our
    code).

6.  Estimate \\\hat{V}\_1=\text{var}({\hat{\boldsymbol{\mu}}})\\
    obtained from simulations and save it as `var_nonprob`.

\(b\) probability part (\\S_B\\ with size \\n_B\\; denoted as `var_prob`
in the result)

This part uses functionalities of the `{survey}` package and the
variance is estimated using the following equation:

\$\$ \hat{V}\_2=\frac{1}{N^2} \sum\_{i=1}^{n_B} \sum\_{j=1}^{n_B}
\frac{\pi\_{i j}-\pi_i \pi_j}{\pi\_{i j}} \frac{m(\boldsymbol{x}\_i;
\hat{\boldsymbol{\beta}})}{\pi_i} \frac{m(\boldsymbol{x}\_i;
\hat{\boldsymbol{\beta}})}{\pi_j}. \$\$

Note that \\\hat{V}\_2\\ in principle can be estimated in various ways
depending on the type of the design and whether population size is known
or not.

## Examples

``` r

data(admin)
data(jvs)
jvs_svy <- svydesign(ids = ~ 1,  weights = ~ weight, strata = ~ size + nace + region, data = jvs)

res_pmm <- method_pmm(y_nons = admin$single_shift,
                      X_nons = model.matrix(~ region + private + nace + size, admin),
                      X_rand = model.matrix(~ region + private + nace + size, jvs),
                      svydesign = jvs_svy)

res_pmm
#> Mass imputation model (PMM approach). Estimated mean: 0.6969 (se: 0.0146)
```
