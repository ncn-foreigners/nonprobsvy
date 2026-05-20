# Mass Imputation Using Non-Parametric Model Method

Model for the outcome for the mass imputation estimator using loess via
[`stats::loess`](https://rdrr.io/r/stats/loess.html). Estimation of the
mean is done using the \\S_B\\ probability sample.

## Usage

``` r
method_npar(
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

  case / frequency weights from non-probability sample (default NULL)

- family_outcome:

  family for the glm model

- start_outcome:

  a place holder (not used in `method_npar`)

- vars_selection:

  whether variable selection should be conducted

- pop_totals:

  a place holder (not used in `method_npar`)

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

  fitted model object returned by
  [`stats::loess`](https://rdrr.io/r/stats/loess.html)

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

  model type (character `"npar"`)

## Details

Analytical variance

The variance of the mean is estimated based on the following approach

\(a\) non-probability part (\\S_A\\ with size \\n_A\\; denoted as
`var_nonprob` in the result)

\$\$ \hat{V}\_1 = \frac{1}{N^2} \sum\_{i=1}^{n_A}
\left\lbrace\hat{g}\_B(\boldsymbol{x}\_i)\right\rbrace^{2} \hat{e}\_i^2,
\$\$

where \\\hat{e}\_i=y_i - \hat{m}(x_i)\\ is the residual and
\\\hat{g}\_B(\boldsymbol{x}\_i) = \left\lbrace \pi_B(\boldsymbol{x}\_i)
\right\rbrace^{-1}\\ can be estimated various ways. In the package we
estimate \\\hat{g}\_B(\boldsymbol{x}\_i)\\ using
\\\pi_B(\boldsymbol{x}\_i)=E(R \| \boldsymbol{x})\\ as suggested by Chen
et al. (2022, p. 6). In particular, we currently support this using
stats::loess`with`"gaussian"\` family.

\(b\) probability part (\\S_B\\ with size \\n_B\\; denoted as `var_prob`
in the result)

This part uses functionalities of the `{survey}` package and the
variance is estimated using the following equation:

\$\$ \hat{V}\_2=\frac{1}{N^2} \sum\_{i=1}^{n_B} \sum\_{j=1}^{n_B}
\frac{\pi\_{i j}-\pi_i \pi_j}{\pi\_{i j}} \frac{\hat{m}(x_i)}{\pi_i}
\frac{\hat{m}(x_j)}{\pi_j}. \$\$

Note that \\\hat{V}\_2\\ in principle can be estimated in various ways
depending on the type of the design and whether population size is known
or not.

## References

Chen, S., Yang, S., & Kim, J. K. (2022). Nonparametric mass imputation
for data integration. Journal of Survey Statistics and Methodology,
10(1), 1-24.

## Examples

``` r

set.seed(123123123)
N <- 10000
n_a <- 500
n_b <- 1000
n_b1 <- 0.7*n_b
n_b2 <- 0.3*n_b
x1 <- rnorm(N, 2, 1)
x2 <- rnorm(N, 2, 1)
y1 <- rnorm(N, 0.3 + 2*x1+ 2*x2, 1)
y2 <- rnorm(N, 0.3 + 0.5*x1^2+ 0.5*x2^2, 1)
strata <- x1 <= 2
pop <- data.frame(x1, x2, y1, y2, strata)
sample_a <- pop[sample(1:N, n_a),]
sample_a$w_a <- N/n_a
sample_a_svy <- svydesign(ids=~1, weights=~w_a, data=sample_a)
pop1 <- subset(pop, strata == TRUE)
pop2 <- subset(pop, strata == FALSE)
sample_b <- rbind(pop1[sample(1:nrow(pop1), n_b1), ],
                  pop2[sample(1:nrow(pop2), n_b2), ])
res_y_npar <- nonprob(outcome = y1 + y2 ~ x1 + x2,
                      data = sample_b,
                      svydesign = sample_a_svy,
                      method_outcome = "npar")
#> Warning: Variable(s) y1, y2 are present in the probability sample (`svydesign`) but are ignored: `nonprob()` uses the outcome only from the non-probability sample; the probability sample contributes covariates and design weights.
res_y_npar
#> A nonprob object
#>  - estimator type: mass imputation
#>  - method: npar (gaussian)
#>  - auxiliary variables source: survey
#>  - vars selection: false
#>  - variance estimator: analytic
#>  - population size fixed: false
#>  - naive (uncorrected) estimators:
#>    - variable y1: 7.5715
#>    - variable y2: 4.6272
#>  - selected estimators:
#>    - variable y1: 8.3408 (se=0.1221, ci=(8.1015, 8.5802))
#>    - variable y2: 5.3156 (se=0.1295, ci=(5.0618, 5.5695))
```
