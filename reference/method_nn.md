# Mass Imputation Using Nearest Neighbours Matching Method

Mass imputation using nearest neighbours approach as described in Yang
et al. (2021). The implementation is currently based on
[RANN::nn2](https://jefferislab.github.io/RANN/reference/nn2.html)
function and thus it uses Euclidean distance for matching units from
\\S_A\\ (non-probability) to \\S_B\\ (probability). Matching ties are
randomized before donor values are aggregated, so tied nearest
neighbours are not selected only by input row order. Estimation of the
mean is done using \\S_B\\ sample.

## Usage

``` r
method_nn(
  y_nons,
  X_nons,
  X_rand,
  svydesign,
  weights = NULL,
  family_outcome = NULL,
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
  sampling probabilities proportional to their inverses.

- family_outcome:

  a placeholder (not used in `method_nn`)

- start_outcome:

  a placeholder (not used in `method_nn`)

- vars_selection:

  whether variable selection should be conducted

- pop_totals:

  a placeholder (not used in `method_nn`)

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

  [`RANN::nn2`](https://jefferislab.github.io/RANN/reference/nn2.html)
  object

- y_nons_pred:

  predicted values for the non-probablity sample (query to itself)

- y_rand_pred:

  predicted values for the probability sample

- coefficients:

  coefficients for the model (if available)

- svydesign:

  an updated `surveydesign2` object (new column `y_hat_MI` is added)

- y_mi_hat:

  estimated population mean for the target variable

- vars_selection:

  whether variable selection was performed (not implemented, for further
  development)

- var_prob:

  variance for the probability sample component (if available)

- var_nonprob:

  variance for the non-probability sample component

- var_tot:

  total variance, if possible it should be `var_prob+var_nonprob` if
  not, just a scalar

- model:

  model type (character `"nn"`)

- family:

  placeholder for the `NN approach` information

## Details

Analytical variance

The variance of the mean is estimated based on the following approach

\(a\) non-probability part (\\S_A\\ with size \\n_A\\; denoted as
`var_nonprob` in the result)

This may be estimated using

\$\$ \hat{V}\_1 =
\frac{1}{N^2}\sum\_{i=1}^{S_A}\frac{1-\hat{\pi}\_B(\boldsymbol{x}\_i)}{\hat{\pi}\_B(\boldsymbol{x}\_i)}\hat{\sigma}^2(\boldsymbol{x}\_i),
\$\$

where \\\hat{\pi}\_B(\boldsymbol{x}\_i)\\ is an estimator of propensity
scores which we currently estimate using \\n_A/N\\ (constant) and
\\\hat{\sigma}^2(\boldsymbol{x}\_i)\\ is estimated using based on the
average of \\(y_i - y_i^\*)^2\\.

Chlebicki et al. (2025, Algorithm 2) proposed non-parametric
mini-bootstrap estimator (without assuming that it is consistent) but
with good finite population properties. This bootstrap can be applied
using `control_inference(nn_exact_se=TRUE)` and can be summarized as
follows:

1.  Sample \\n_A\\ units from \\S_A\\ with replacement to create
    \\S_A'\\. If non-constant pseudo-weights are supplied through
    `weights`, sampling probabilities are proportional to their
    inverses; equal weights use uniform resampling.

2.  Match units from \\S_B\\ to \\S_A'\\ to obtain predictions
    \\y^\*\\=\\{k}^{-1}\sum\_{k}y_k\\.

3.  Estimate \\\hat{\mu}=\frac{1}{N} \sum\_{i \in S_B} d_i y_i^\*\\.

4.  Repeat steps 1-3 \\M\\ times (we set \\M=50\\ in our simulations;
    this is hard-coded).

5.  Estimate \\\hat{V}\_1=\text{var}({\hat{\boldsymbol{\mu}}})\\
    obtained from simulations and save it as `var_nonprob`.

\(b\) probability part (\\S_B\\ with size \\n_B\\; denoted as `var_prob`
in the result)

This part uses functionalities of the `{survey}` package and the
variance is estimated using the following equation:

\$\$ \hat{V}\_2=\frac{1}{N^2} \sum\_{i=1}^n \sum\_{j=1}^n \frac{\pi\_{i
j}-\pi_i \pi_j}{\pi\_{i j}} \frac{y_i^\*}{\pi_i} \frac{y_j^\*}{\pi_j},
\$\$

where \\y^\*\_i\\ and \\y_j^\*\\ are values imputed imputed as an
average of \\k\\-nearest neighbour, i.e. \\{k}^{-1}\sum\_{k}y_k\\. Note
that \\\hat{V}\_2\\ in principle can be estimated in various ways
depending on the type of the design and whether population size is known
or not.

## References

Yang, S., Kim, J. K., & Hwang, Y. (2021). Integration of data from
probability surveys and big found data for finite population inference
using mass imputation. Survey Methodology, June 2021 29 Vol. 47, No. 1,
pp. 29-58

Chlebicki, P., Chrostowski, Ł., & Beręsewicz, M. (2025). Data
integration of non-probability and probability samples with predictive
mean matching. arXiv preprint arXiv:2403.13750.

## Examples

``` r

data(admin)
data(jvs)
jvs_svy <- svydesign(ids = ~ 1,  weights = ~ weight, strata = ~ size + nace + region, data = jvs)

res_nn <- method_nn(y_nons = admin$single_shift,
                    X_nons = model.matrix(~ region + private + nace + size, admin),
                    X_rand = model.matrix(~ region + private + nace + size, jvs),
                    svydesign = jvs_svy)

res_nn
#> Mass imputation model (NN approach). Estimated mean: 0.6800 (se: 0.0157)
```
