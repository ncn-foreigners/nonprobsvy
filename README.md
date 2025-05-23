
<!-- README.md is generated from README.Rmd. Please edit that file -->

# `nonprobsvy`: an R package for modern statistical inference methods based on non-probability samples <img src="man/figures/logo.png" align="right" width="150"/>

<!-- badges: start -->

[![R-CMD-check](https://github.com/ncn-foreigners/nonprobsvy/workflows/R-CMD-check/badge.svg)](https://github.com/ncn-foreigners/nonprobsvy/actions)
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
  constraints ([Y. Chen, Li, and Wu 2020](#ref-chen2020)),
- mass imputation estimators based on nearest neighbours ([Yang, Kim,
  and Hwang 2021](#ref-yang2021)), predictive mean matching ([Chlebicki,
  Chrostowski, and Beręsewicz 2025](#ref-chlebicki2025)), non-parametric
  ([S. Chen, Yang, and Kim 2022](#ref-chen2022nonparametric)) and
  regression imputation ([Kim et al. 2021](#ref-kim2021)),
- doubly robust estimators ([Y. Chen, Li, and Wu 2020](#ref-chen2020))
  with bias minimization ([Yang, Kim, and Song 2020](#ref-yang2020)).

The package allows for:

- variable section in high-dimensional space using SCAD ([Yang, Kim, and
  Song 2020](#ref-yang2020)), Lasso and MCP penalty (via the `ncvreg`,
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

- in the draft (and not proofread) version of the book [Modern inference
  methods for non-probability samples with
  R](https://ncn-foreigners.github.io/nonprobsvy-book/),
- in example codes that reproduce papers available on github in the
  repository [software
  tutorials](https://github.com/ncn-foreigners/software-tutorials).

## Installation

You can install the recent version of `nonprobsvy` package from main
branch [Github](https://github.com/ncn-foreigners/nonprobsvy) with:

``` r
remotes::install_github("ncn-foreigners/nonprobsvy")
```

or install the stable version from
[CRAN](https://CRAN.R-project.org/package=nonprobsvy)

``` r
install.packages("nonprobsvy")
```

or development version from the `dev` branch

``` r
remotes::install_github("ncn-foreigners/nonprobsvy@dev")
```

## Basic idea

Consider the following setting where two samples are available:
non-probability (denoted as $S_A$ ) and probability (denoted as $S_B$)
where set of auxiliary variables (denoted as $\boldsymbol{X}$) is
available for both sources while $Y$ and $\boldsymbol{d}$ (or
$\boldsymbol{w}$) is present only in probability sample.

| Sample |  | Auxiliary variables $\boldsymbol{X}$ | Target variable $Y$ | Design ($\boldsymbol{d}$) or calibrated ($\boldsymbol{w}$) weights |
|----|---:|:--:|:--:|:--:|
| $S_A$ (non-probability) | 1 | $\checkmark$ | $\checkmark$ | ? |
|  | … | $\checkmark$ | $\checkmark$ | ? |
|  | $n_A$ | $\checkmark$ | $\checkmark$ | ? |
| $S_B$ (probability) | $n_A+1$ | $\checkmark$ | ? | $\checkmark$ |
|  | … | $\checkmark$ | ? | $\checkmark$ |
|  | $n_A+n_B$ | $\checkmark$ | ? | $\checkmark$ |

## Basic functionalities

Suppose $Y$ is the target variable, $\boldsymbol{X}$ is a matrix of
auxiliary variables, $R$ is the inclusion indicator. Then, if we are
interested in estimating the mean $\bar{\tau}_Y$ or the sum $\tau_Y$ of
the of the target variable given the observed data set
$(y_k, \boldsymbol{x}_k, R_k)$, we can approach this problem with the
possible scenarios:

- unit-level data is available for the non-probability sample $S_{A}$,
  i.e. $(y_{k}, \boldsymbol{x}_{k})$ is available for all units
  $k \in S_{A}$, and population-level data is available for
  $\boldsymbol{x}_{1}, ..., \boldsymbol{x}_{p}$, denoted as
  $\tau_{x_{1}}, \tau_{x_{2}}, ..., \tau_{x_{p}}$ and population size
  $N$ is known. We can also consider situations where population data
  are estimated (e.g. on the basis of a survey to which we do not have
  access),
- unit-level data is available for the non-probability sample $S_A$ and
  the probability sample $S_B$, i.e. $(y_k, \boldsymbol{x}_k, R_k)$ is
  determined by the data. is determined by the data: $R_k=1$ if
  $k \in S_A$ otherwise $R_k=0$, $y_k$ is observed only for sample $S_A$
  and $\boldsymbol{x}_k$ is observed in both in both $S_A$ and $S_B$,

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
  outcome = y ~ x1 + x2 + ... + xk,
  data = nonprob,
  pop_totals = c(`(Intercept)`= N,
                 x1 = tau_x1,
                 x2 = tau_x2,
                 ...,
                 xk = tau_xk),
  method_outcome = "glm",
  family_outcome = "gaussian"
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
  selection =  ~ x1 + x2 + ... + xk, 
  target = ~ y, 
  data = nonprob, 
  pop_totals = c(`(Intercept)` = N, 
                 x1 = tau_x1, 
                 x2 = tau_x2, 
                 ..., 
                 xk = tau_xk), 
  method_selection = "logit"
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
  selection =  ~ x1 + x2 + ... + xk, 
  target = ~ y, 
  data = nonprob, 
  pop_totals = c(`(Intercept)`= N, 
                 x1 = tau_x1, 
                 x2 = tau_x2, 
                 ..., 
                 xk = tau_xk), 
  method_selection = "logit", 
  control_selection = control_sel(est_method = "gee", gee_h_fun = 1)
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
  selection = ~ x1 + x2 + ... + xk, 
  outcome = y ~ x1 + x2 + …, + xk, 
  pop_totals = c(`(Intercept)` = N, 
                 x1 = tau_x1, 
                 x2 = tau_x2, 
                 ..., 
                 xk = tau_xk), 
  svydesign = prob, 
  method_outcome = "glm", 
  family_outcome = "gaussian"
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
  outcome = y ~ x1 + x2 + ... + xk, 
  data = nonprob, 
  svydesign = prob, 
  method_outcome = "glm", 
  family_outcome = "gaussian"
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
  outcome = y ~ x1 + x2 + ... + xk, 
  data = nonprob, 
  svydesign = prob, 
  method_outcome = "nn", 
  family_outcome = "gaussian", 
  control_outcome = control_outcome(k = 2)
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
  outcome = y ~ x1 + x2 + ... + xk, 
  data = nonprob, 
  svydesign = prob, 
  method_outcome = "pmm", 
  family_outcome = "gaussian"
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
  outcome = y ~ x1 + x2 + ... + xk, 
  data = nonprob, 
  svydesign = prob, 
  method_outcome = "pmm", 
  family_outcome = "gaussian", 
  control_outcome = control_out(penalty = "lasso"), 
  control_inference = control_inf(vars_selection = TRUE)
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
  selection =  ~ x1 + x2 + ... + xk, 
  target = ~ y, 
  data = nonprob, 
  svydesign = prob, 
  method_selection = "logit"
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
  selection =  ~ x1 + x2 + ... + xk, 
  target = ~ y, 
  data = nonprob, 
  svydesign = prob, 
  method_selection = "logit", 
  control_selection = control_sel(est_method = "gee", gee_h_fun = 1)
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
  selection =  ~ x1 + x2 + ... + xk, 
  target = ~ y, 
  data = nonprob, 
  svydesign = prob, 
  method_outcome = "pmm", 
  family_outcome = "gaussian", 
  control_inference = control_inf(vars_selection = TRUE)
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
  selection = ~ x1 + x2 + ... + xk, 
  outcome = y ~ x1 + x2 + ... + xk, 
  data = nonprob, 
  svydesign = prob, 
  method_outcome = "glm", 
  family_outcome = "gaussian"
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
  selection = ~ x1 + x2 + ... + xk, 
  outcome = y ~ x1 + x2 + ... + xk, 
  data = nonprob, 
  svydesign = prob,
  method_outcome = "glm", 
  family_outcome = "gaussian", 
  control_inference = control_inf(
    vars_selection = TRUE, 
    bias_correction = TRUE
  )
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
#>  - auxiliary variables source: survey
#>  - vars selection: false
#>  - variance estimator: analytic
#>  - population size fixed: false
#>  - naive (uncorrected) estimators:
#>    - variable y1: 3.1817
#>    - variable y2: 1.8087
#>  - selected estimators:
#>    - variable y1: 2.9500 (se=0.0414, ci=(2.8689, 3.0312))
#>    - variable y2: 1.5762 (se=0.0498, ci=(1.4786, 1.6739))
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
#>  - estimator type: inverse probability weighting
#>  - method: logit (mle)
#>  - auxiliary variables source: survey
#>  - vars selection: false
#>  - variance estimator: analytic
#>  - population size fixed: false
#>  - naive (uncorrected) estimators:
#>    - variable y1: 3.1817
#>    - variable y2: 1.8087
#>  - selected estimators:
#>    - variable y1: 2.9981 (se=0.0137, ci=(2.9713, 3.0249))
#>    - variable y2: 1.5906 (se=0.0137, ci=(1.5639, 1.6174))
```

## Funding

Work on this package is supported by the National Science Centre, OPUS
20 grant no. 2020/39/B/HS4/00941.

## References (selected)

<div id="refs" class="references csl-bib-body hanging-indent"
entry-spacing="0">

<div id="ref-chen2022nonparametric" class="csl-entry">

Chen, Sixia, Shu Yang, and Jae Kwang Kim. 2022. “Nonparametric Mass
Imputation for Data Integration.” *Journal of Survey Statistics and
Methodology* 10 (1): 1–24.

</div>

<div id="ref-chen2020" class="csl-entry">

Chen, Yilin, Pengfei Li, and Changbao Wu. 2020. “Doubly Robust Inference
With Nonprobability Survey Samples.” *Journal of the American
Statistical Association* 115 (532): 2011–21.
<https://doi.org/10.1080/01621459.2019.1677241>.

</div>

<div id="ref-chlebicki2025" class="csl-entry">

Chlebicki, Piotr, Łukasz Chrostowski, and Maciej Beręsewicz. 2025. “Data
Integration of Non-Probability and Probability Samples with Predictive
Mean Matching.” <https://arxiv.org/abs/2403.13750>.

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

———. 2023. “Survey: Analysis of Complex Survey Samples.”

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
