
<!-- README.md is generated from README.Rmd. Please edit that file -->

# `nonprobsvy`: an R package for modern statistical inference methods based on non-probability samples <img src="man/figures/logo.png" align="right" width="150"/>

<!-- badges: start -->

[![R-CMD-check](https://github.com/ncn-foreigners/nonprobsvy/workflows/R-CMD-check/badge.svg)](https://github.com/ncn-foreigners/nonprobsvy/actions)
[![revdep-check](https://github.com/ncn-foreigners/nonprobsvy/actions/workflows/revdep-check.yaml/badge.svg)](https://github.com/ncn-foreigners/nonprobsvy/actions/workflows/revdep-check.yaml)
[![Codecov test
coverage](https://codecov.io/gh/ncn-foreigners/nonprobsvy/branch/main/graph/badge.svg)](https://app.codecov.io/gh/ncn-foreigners/nonprobsvy?branch=main)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.10280114.svg)](https://doi.org/10.5281/zenodo.10280114)
[![CRAN
status](https://www.r-pkg.org/badges/version/nonprobsvy)](https://CRAN.R-project.org/package=nonprobsvy)
[![CRAN
downloads](https://cranlogs.r-pkg.org/badges/grand-total/nonprobsvy)](https://cran.r-project.org/package=nonprobsvy)
[![CRAN
downloads](https://cranlogs.r-pkg.org/badges/nonprobsvy)](https://cran.r-project.org/package=nonprobsvy)
[![Mentioned in Awesome Official
Statistics](https://awesome.re/mentioned-badge.svg)](https://github.com/SNStatComp/awesome-official-statistics-software)

<!-- badges: end -->

## Basic information

The goal of this package is to provide R users access to modern methods
for non-probability samples when auxiliary information from the
population or probability sample is available:

- inverse probability weighting estimators with possible calibration
  constraints ([Chen et al. 2020](#ref-chen2020)),
- mass imputation estimators based on nearest neighbours ([Yang et al.
  2021](#ref-yang2021)), predictive mean matching ([Chlebicki et al.
  2025](#ref-chlebicki2025)), non-parametric ([Chen et al.
  2022](#ref-chen2022nonparametric)) and regression imputation ([Kim et
  al. 2021](#ref-kim2021)),
- doubly robust estimators ([Chen et al. 2020](#ref-chen2020)) with bias
  minimization ([Yang et al. 2020](#ref-yang2020)).

The package allows for:

- variable section in high-dimensional space using SCAD ([Yang et al.
  2020](#ref-yang2020)), Lasso and MCP penalty (via the `ncvreg`,
  `Rcpp`, `RcppArmadillo` packages),
- estimation of variance using analytical and bootstrap approach (see Wu
  ([2023](#ref-wu2023))),
- integration with the `survey` and `srvyr` packages when probability
  sample is available ([Lumley 2004](#ref-Lumley2004),
  [2023](#ref-Lumley2023); [Freedman Ellis and Schneider
  2024](#ref-srvyr2024)),
- different links for selection (`logit`, `probit` and `cloglog`) and
  outcome (`gaussian`, `binomial` and `poisson`) variables.

Details on the use of the package can be found:

- please cite the following paper: Chrostowski, Ł., Chlebicki, P., &
  Beręsewicz, M(2026). nonprobsvy: An R package for modern methods for
  non-probability surveys. Journal of Statistical Software, 117(2),
  1-37. <https://doi.org/10.18637/jss.v117.i02>
- in example codes that reproduce papers available on github in the
  repository [software
  tutorials](https://github.com/ncn-foreigners/software-tutorials).

## Installation

You can install the recent version of `nonprobsvy` package from main
branch [Github](https://github.com/ncn-foreigners/nonprobsvy) with:

``` r
pak::pkg_install("ncn-foreigners/nonprobsvy")
```

or install the stable version from
[CRAN](https://CRAN.R-project.org/package=nonprobsvy)

``` r
install.packages("nonprobsvy")
```

or development version from the `dev` branch

``` r
pak::pkg_install("ncn-foreigners/nonprobsvy@dev")
```

## Basic idea

Consider the following setting where two samples are available:
non-probability (denoted as $S_{\text{NP}}$) and probability (denoted as
$S_{\text{P}}$) where a set of auxiliary variables (denoted as
$\boldsymbol{X}$) is available for both sources, the target variable $Y$
is observed in the non-probability sample, and design or calibrated
weights ($\boldsymbol{d}$ or $\boldsymbol{w}$) are observed in the
probability sample.

| Sample |  | Auxiliary variables $\boldsymbol{X}$ | Target variable $Y$ | Design ($\boldsymbol{d}$) or calibrated ($\boldsymbol{w}$) weights |
|----|---:|:--:|:--:|:--:|
| $S_{\text{NP}}$ (non-probability) | 1 | $\checkmark$ | $\checkmark$ | ? |
|  | … | $\checkmark$ | $\checkmark$ | ? |
|  | $n_{\text{NP}}$ | $\checkmark$ | $\checkmark$ | ? |
| $S_{\text{P}}$ (probability) | $n_{\text{NP}}+1$ | $\checkmark$ | ? | $\checkmark$ |
|  | … | $\checkmark$ | ? | $\checkmark$ |
|  | $n_{\text{NP}}+n_{\text{P}}$ | $\checkmark$ | ? | $\checkmark$ |

The current implementation does not use target-variable values from the
probability sample. Data structures where $Y$ is observed in both
samples, or where overlapping units must be linked across samples, are
not currently implemented. The `dependence` and `key` arguments in
`control_sel()` are reserved for future overlap handling and currently
raise a “not yet implemented” error if set.

## Basic functionalities

Suppose $Y$ is the target variable, $\boldsymbol{X}$ is a matrix of
auxiliary variables, $R$ is the inclusion indicator. Then, if we are
interested in estimating the mean $\bar{\tau}_Y$ or the sum $\tau_Y$ of
the of the target variable given the observed data set
$(y_k, \boldsymbol{x}_k, I_{\text{NP}, k})$, we can approach this
problem with the possible scenarios:

- unit-level data is available for the non-probability sample
  $S_{\text{NP}}$, i.e. $(y_k,\boldsymbol{x}_k)$ is available for all
  units $k \in S_{\text{NP}}$, and population-level data is available
  for $\boldsymbol{x}_1,\ldots,\boldsymbol{x}_p$, denoted as
  $\tau_{x_{1}},\tau_{x_{2}},\ldots,\tau_{x_{p}}$ and population size
  $N$ is known. We can also consider situations where population data
  are estimated (e.g. on the basis of a survey to which we do not have
  access),
- unit-level data is available for the non-probability sample
  $S_{\text{NP}}$ and the probability sample $S_{\text{P}}$, i.e.
  $(\boldsymbol{x}_k,I_{\text{NP}, k})$ is determined by the data:
  $I_{\text{NP}, k}=1$ if $k \in S_{\text{NP}}$ otherwise
  $I_{\text{NP}, k}=0$, $y_k$ is observed only for sample
  $S_{\text{NP}}$ and $\boldsymbol{x}_k$ is observed in both
  $S_{\text{NP}}$ and $S_{\text{P}}$.

Supported target-variable types depend on the estimator family:

| Estimator family | Supported target variable `Y` |
|----|----|
| IPW | Numeric targets whose population mean is meaningful, including continuous, count, and 0/1 binary variables. No outcome model is fitted. |
| Mass imputation with `method_outcome = "glm"` | Continuous, count, or binary variables through `family_outcome = "gaussian"`, `"poisson"`, or `"binomial"`. |
| Mass imputation with `method_outcome = "nn"`, `"pmm"`, or `"npar"` | Numeric targets; categorical, ordinal, survival, and other structured outcomes are not supported. |
| Doubly robust | GLM outcome models only; use `family_outcome = "gaussian"`, `"poisson"`, or `"binomial"`. |

The compact examples below use the built-in `admin` non-probability
sample and `jvs` probability sample.

``` r
library(survey)
library(nonprobsvy)
data(admin)
data(jvs)

prob <- svydesign(
  ids = ~1,
  weights = ~weight,
  strata = ~size + nace + region,
  data = jvs
)
pop_totals <- colSums(model.matrix(~region + private + nace + size, jvs) * jvs$weight)
```

### When unit-level data is available for non-probability survey only

<table class='table'>

<tr>

<th>

Estimator
</th>

<th>

Example code
</th>

<tr>

<tr>

<td>

Mass imputation based on regression imputation
</td>

<td>

``` r
nonprob(
  outcome = single_shift ~ region + private + nace + size,
  data = admin,
  pop_totals = pop_totals,
  method_outcome = "glm",
  family_outcome = "binomial",
  se = FALSE
)
```

</td>

<tr>

<tr>

<td>

Inverse probability weighting
</td>

<td>

``` r
nonprob(
  selection = ~region + private + nace + size,
  target = ~single_shift,
  data = admin,
  pop_totals = pop_totals,
  method_selection = "logit",
  se = FALSE
)
```

</td>

<tr>

<tr>

<td>

Inverse probability weighting with calibration constraint
</td>

<td>

``` r
nonprob(
  selection = ~region + private + nace + size,
  target = ~single_shift,
  data = admin,
  pop_totals = pop_totals,
  method_selection = "logit",
  control_selection = control_sel(est_method = "gee", gee_h_fun = 1),
  se = FALSE
)
```

</td>

<tr>

<tr>

<td>

Doubly robust estimator
</td>

<td>

``` r
nonprob(
  selection = ~region + private + nace + size,
  outcome = single_shift ~ region + private + nace + size,
  data = admin,
  pop_totals = pop_totals,
  method_outcome = "glm", 
  family_outcome = "binomial",
  se = FALSE
)
```

</td>

<tr>

</table>

### When unit-level data are available for both surveys

<table class='table'>

<tr>

<th>

Estimator
</th>

<th>

Example code
</th>

<tr>

<tr>

<td>

Mass imputation based on regression imputation
</td>

<td>

``` r
nonprob(
  outcome = single_shift ~ region + private + nace + size,
  data = admin,
  svydesign = prob,
  method_outcome = "glm", 
  family_outcome = "binomial",
  se = FALSE
)
```

</td>

<tr>

<tr>

<td>

Mass imputation based on nearest neighbour imputation
</td>

<td>

``` r
nonprob(
  outcome = single_shift ~ region + private + nace + size,
  data = admin,
  svydesign = prob,
  method_outcome = "nn", 
  control_outcome = control_out(k = 2),
  se = FALSE
)
```

</td>

<tr>

<tr>

<td>

Mass imputation based on predictive mean matching
</td>

<td>

``` r
nonprob(
  outcome = single_shift ~ region + private + nace + size,
  data = admin,
  svydesign = prob,
  method_outcome = "pmm", 
  se = FALSE
)
```

</td>

<tr>

<tr>

<td>

Mass imputation based on regression imputation with variable selection
(LASSO)
</td>

<td>

``` r
nonprob(
  outcome = single_shift ~ region + private + nace + size,
  data = admin,
  svydesign = prob,
  method_outcome = "glm",
  family_outcome = "binomial",
  control_outcome = control_out(penalty = "lasso"), 
  control_inference = control_inf(vars_selection = TRUE),
  se = FALSE
)
```

</td>

<tr>

<tr>

<td>

Inverse probability weighting
</td>

<td>

``` r
nonprob(
  selection = ~region + private + nace + size,
  target = ~single_shift,
  data = admin,
  svydesign = prob,
  method_selection = "logit",
  se = FALSE
)
```

</td>

<tr>

<tr>

<td>

Inverse probability weighting with calibration constraint
</td>

<td>

``` r
nonprob(
  selection = ~region + private + nace + size,
  target = ~single_shift,
  data = admin,
  svydesign = prob,
  method_selection = "logit", 
  control_selection = control_sel(est_method = "gee", gee_h_fun = 1),
  se = FALSE
)
```

</td>

<tr>

<tr>

<td>

Inverse probability weighting with calibration constraint with variable
selection (SCAD)
</td>

<td>

``` r
nonprob(
  selection = ~region + private + nace + size,
  target = ~single_shift,
  data = admin,
  svydesign = prob,
  method_selection = "logit",
  control_selection = control_sel(penalty = "SCAD"),
  control_inference = control_inf(vars_selection = TRUE),
  se = FALSE
)
```

</td>

<tr>

<tr>

<td>

Doubly robust estimator
</td>

<td>

``` r
nonprob(
  selection = ~region + private + nace + size,
  outcome = single_shift ~ region + private + nace + size,
  data = admin,
  svydesign = prob,
  method_outcome = "glm", 
  family_outcome = "binomial",
  se = FALSE
)
```

</td>

<tr>

<tr>

<td>

Doubly robust estimator with variable selection (SCAD) and bias
minimization
</td>

<td>

``` r
nonprob(
  selection = ~region + private + nace + size,
  outcome = single_shift ~ region + private + nace + size,
  data = admin,
  svydesign = prob,
  method_outcome = "glm", 
  family_outcome = "binomial",
  control_inference = control_inf(
    vars_selection = TRUE, 
    vars_combine = TRUE,
    bias_correction = TRUE
  ),
  se = FALSE
)
```

</td>

<tr>

</table>

## Examples

Simulate example data from the following paper: Kim, Jae Kwang, and
Zhonglei Wang. “Sampling techniques for big data analysis.”
International Statistical Review 87 (2019): S177-S191 \[section 5.2\]

``` r
library(survey)
library(nonprobsvy)

set.seed(1234567890)
N <- 1e6 ## 1000000
n <- 1000
x1 <- rnorm(n = N, mean = 1, sd = 1)
x2 <- rexp(n = N, rate = 1)
epsilon <- rnorm(n = N) # rnorm(N)
y1 <- 1 + x1 + x2 + epsilon
y2 <- 0.5*(x1 - 0.5)^2 + x2 + epsilon
p1 <- exp(x2)/(1+exp(x2))
p2 <- exp(-0.5+0.5*(x2-2)^2)/(1+exp(-0.5+0.5*(x2-2)^2))
flag_bd1 <- rbinom(n = N, size = 1, prob = p1)
flag_srs <- as.numeric(1:N %in% sample(1:N, size = n))
base_w_srs <- N/n
population <- data.frame(x1,x2,y1,y2,p1,p2,base_w_srs, flag_bd1, flag_srs, pop_size = N)
base_w_bd <- N/sum(population$flag_bd1)
```

Declare `svydesign` object with `survey` package

``` r
sample_prob <- svydesign(ids= ~1, weights = ~ base_w_srs, 
                         data = subset(population, flag_srs == 1),
                         fpc = ~ pop_size)

sample_prob
#> Independent Sampling design
#> svydesign(ids = ~1, weights = ~base_w_srs, data = subset(population, 
#>     flag_srs == 1), fpc = ~pop_size)
```

or with the `srvyr` package

``` r
sample_prob <- srvyr::as_survey_design(.data = subset(population, flag_srs == 1),
                                       weights = base_w_srs)

sample_prob
```

``` r
Independent Sampling design (with replacement)
Called via srvyr
Sampling variables:
Data variables: 
  - x1 (dbl), x2 (dbl), y1 (dbl), y2 (dbl), p1 (dbl), p2 (dbl), base_w_srs (dbl), flag_bd1 (int), flag_srs (dbl)
```

Estimate population mean of `y1` based on doubly robust estimator using
IPW with calibration constraints and we specify that auxiliary variables
should not be combined for the inference.

``` r
result_dr <- nonprob(
  selection = ~ x2,
  outcome = y1 + y2 ~ x1 + x2,
  data = subset(population, flag_bd1 == 1),
  svydesign = sample_prob
)
```

Results

``` r
result_dr
#> A nonprob object
#>  - estimator type: doubly robust
#>  - method: glm (gaussian)
#>  - IPW point estimator: Hajek (denominator: 1025063)
#>  - auxiliary variables source: survey
#>  - vars selection: false
#>  - variance estimator: analytic
#>  - population size fixed: false
#>  - naive (uncorrected) estimators:
#>    - variable y1: 3.1817
#>    - variable y2: 1.8087
#>  - selected estimators:
#>    - variable y1: 2.9500 (se=0.0414, ci=(2.8689, 3.0312))
#>    - variable y2: 1.5762 (se=0.0313, ci=(1.5150, 1.6375))
```

Mass imputation estimator

``` r
result_mi <- nonprob(
  outcome = y1 + y2 ~ x1 + x2,
  data = subset(population, flag_bd1 == 1),
  svydesign = sample_prob
)
```

Results

``` r
result_mi
#> A nonprob object
#>  - estimator type: mass imputation
#>  - method: glm (gaussian)
#>  - auxiliary variables source: survey
#>  - vars selection: false
#>  - variance estimator: analytic
#>  - population size fixed: false
#>  - naive (uncorrected) estimators:
#>    - variable y1: 3.1817
#>    - variable y2: 1.8087
#>  - selected estimators:
#>    - variable y1: 2.9498 (se=0.0420, ci=(2.8675, 3.0321))
#>    - variable y2: 1.5760 (se=0.0326, ci=(1.5122, 1.6398))
```

Inverse probability weighting estimator

For IPW, MLE without a fixed `pop_size`, `pop_totals`, or `pop_means`
uses the Hajek-type estimator. Supplying a fixed population size uses
the Horvitz-Thompson-type estimator, while IPW-GEE with a reference
survey uses `sum(weights(svydesign))` as the denominator.

``` r
result_ipw <- nonprob(
  selection = ~ x2,
  target = ~y1+y2,
  data = subset(population, flag_bd1 == 1),
  svydesign = sample_prob)
```

Results

``` r
result_ipw
#> A nonprob object
#>  - estimator type: IPW (Hajek, denominator: 1025063)
#>  - method: logit (mle)
#>  - auxiliary variables source: survey
#>  - vars selection: false
#>  - variance estimator: analytic
#>  - population size fixed: false
#>  - naive (uncorrected) estimators:
#>    - variable y1: 3.1817
#>    - variable y2: 1.8087
#>  - selected estimators:
#>    - variable y1: 2.9248 (se=0.0500, ci=(2.8269, 3.0227))
#>    - variable y2: 1.5517 (se=0.0499, ci=(1.4539, 1.6496))
```

## Funding

Work on this package is supported by the National Science Centre, OPUS
20 grant no. 2020/39/B/HS4/00941.

## References (selected)

<div id="refs" class="references csl-bib-body hanging-indent">

<div id="ref-chen2022nonparametric" class="csl-entry">

Chen, Sixia, Shu Yang, and Jae Kwang Kim. 2022. “Nonparametric Mass
Imputation for Data Integration.” *Journal of Survey Statistics and
Methodology* 10 (1): 1–24.

</div>

<div id="ref-chen2020" class="csl-entry">

Chen, Yilin, Pengfei Li, and Changbao Wu. 2020. “Doubly Robust Inference
with Nonprobability Survey Samples.” *Journal of the American
Statistical Association* 115 (532): 2011–21.
<https://doi.org/10.1080/01621459.2019.1677241>.

</div>

<div id="ref-chlebicki2025" class="csl-entry">

Chlebicki, Piotr, Łukasz Chrostowski, and Maciej Beręsewicz. 2025. *Data
Integration of Non-Probability and Probability Samples with Predictive
Mean Matching*. <https://arxiv.org/abs/2403.13750>.

</div>

<div id="ref-srvyr2024" class="csl-entry">

Freedman Ellis, Greg, and Ben Schneider. 2024. *Srvyr: ’Dplyr’-Like
Syntax for Summary Statistics of Survey Data*.
<https://CRAN.R-project.org/package=srvyr>.

</div>

<div id="ref-kim2021" class="csl-entry">

Kim, Jae Kwang, Seho Park, Yilin Chen, and Changbao Wu. 2021. “Combining
Non-Probability and Probability Survey Samples Through Mass Imputation.”
*Journal of the Royal Statistical Society Series A: Statistics in
Society* 184 (3): 941–63. <https://doi.org/10.1111/rssa.12696>.

</div>

<div id="ref-Lumley2004" class="csl-entry">

Lumley, Thomas. 2004. “Analysis of Complex Survey Samples.” *Journal of
Statistical Software* 9 (1): 1–19.

</div>

<div id="ref-Lumley2023" class="csl-entry">

Lumley, Thomas. 2023. *Survey: Analysis of Complex Survey Samples*.

</div>

<div id="ref-wu2023" class="csl-entry">

Wu, Changbao. 2023. “Statistical Inference with Non-Probability Survey
Samples.” *Survey Methodology* 48 (2): 283–311.
<https://www150.statcan.gc.ca/n1/pub/12-001-x/2022002/article/00002-eng.htm>.

</div>

<div id="ref-yang2021" class="csl-entry">

Yang, Shu, Jae Kwang Kim, and Youngdeok Hwang. 2021. “Integration of
Data from Probability Surveys and Big Found Data for Finite Population
Inference Using Mass Imputation.” *Survey Methodology* 47 (1): 29–58.
<https://www150.statcan.gc.ca/n1/pub/12-001-x/2021001/article/00004-eng.htm>.

</div>

<div id="ref-yang2020" class="csl-entry">

Yang, Shu, Jae Kwang Kim, and Rui Song. 2020. “Doubly Robust Inference
When Combining Probability and Non-Probability Samples with High
Dimensional Data.” *Journal of the Royal Statistical Society Series B:
Statistical Methodology* 82 (2): 445–65.
<https://doi.org/10.1111/rssb.12354>.

</div>

</div>
