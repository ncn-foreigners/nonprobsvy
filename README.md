
<!-- README.md is generated from README.Rmd. Please edit that file -->

# `nonprobsvy`: an R package for modern statistical inference methods based on non-probability samples

<!-- badges: start -->

[![R-CMD-check](https://github.com/ncn-foreigners/nonprobsvy/workflows/R-CMD-check/badge.svg)](https://github.com/ncn-foreigners/nonprobsvy/actions)
[![Codecov test
coverage](https://codecov.io/gh/ncn-foreigners/nonprobsvy/branch/main/graph/badge.svg)](https://app.codecov.io/gh/ncn-foreigners/nonprobsvy?branch=main)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.10280114.svg)](https://doi.org/10.5281/zenodo.10280114)
[![CRAN
status](https://www.r-pkg.org/badges/version/nonprobsvy)](https://CRAN.R-project.org/package=nonprobsvy)
[![CRAN
downloads](https://cranlogs.r-pkg.org/badges/nonprobsvy)](https://cran.r-project.org/package=nonprobsvy)

<!-- badges: end -->

## Basic information

The goal of this package is to provide R users access to modern methods
for non-probability samples when auxiliary information from the
population or probability sample is available:

- inverse probability weighting estimators with possible calibration
  constraints ([Chen, Li, and Wu 2020](#ref-chen2020)),
- mass imputation estimators based in nearest neighbours ([Yang, Kim,
  and Hwang 2021](#ref-yang2021)), predictive mean matching and
  regression imputation ([Kim et al. 2021](#ref-kim2021)),
- doubly robust estimators with bias minimization Yang, Kim, and Song
  ([2020](#ref-yang2020)).

The package allows for:

- variable section in high-dimensional space using SCAD ([Yang, Kim, and
  Song 2020](#ref-yang2020)), Lasso and MCP penalty (via the `ncvreg`,
  `Rcpp`, `RcppArmadillo` packages),
- estimation of variance using analytical and bootstrap approach (see Wu
  ([2023](#ref-wu2023))),
- integration with the `survey` package when probability sample is
  available Lumley ([2023](#ref-Lumley2023)),
- different links for selection (`logit`, `probit` and `cloglog`) and
  outcome (`gaussian`, `binomial` and `poisson`) variables.

Details on use of the package be found:

- on the draft (and not proofread) version the book [Modern inference
  methods for non-probability samples with
  R](https://ncn-foreigners.github.io/nonprobsvy-book/),
- example codes that reproduce papers are available at github in the
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

| Sample                  |           | Auxiliary variables $\boldsymbol{X}$ | Target variable $Y$ | Design ($\boldsymbol{d}$) or calibrated ($\boldsymbol{w}$) weights |
|-------------------------|----------:|:------------------------------------:|:-------------------:|:------------------------------------------------------------------:|
| $S_A$ (non-probability) |         1 |             $\checkmark$             |    $\checkmark$     |                                 ?                                  |
|                         |         … |             $\checkmark$             |    $\checkmark$     |                                 ?                                  |
|                         |     $n_A$ |             $\checkmark$             |    $\checkmark$     |                                 ?                                  |
| $S_B$ (probability)     |   $n_A+1$ |             $\checkmark$             |          ?          |                            $\checkmark$                            |
|                         |         … |             $\checkmark$             |          ?          |                            $\checkmark$                            |
|                         | $n_A+n_B$ |             $\checkmark$             |          ?          |                            $\checkmark$                            |

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
  control_selection = controlSel(est_method_sel = "gee", h = 1)
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
  control_outcome = controlOutcome(k = 2)
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
  control_outcome = controlOut(penalty = "lasso"), 
  control_inference = controlInf(vars_selection = TRUE)
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
  control_selection = controlSel(est_method_sel = "gee", h = 1)
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
  control_inference = controlInf(vars_selection = TRUE)
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
  control_inference = controlInf(
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
population <- data.frame(x1,x2,y1,y2,p1,p2,base_w_srs, flag_bd1, flag_srs)
base_w_bd <- N/sum(population$flag_bd1)
```

Declare `svydesign` object with `survey` package

``` r
sample_prob <- svydesign(ids= ~1, weights = ~ base_w_srs, 
                         data = subset(population, flag_srs == 1))
```

Estimate population mean of `y1` based on doubly robust estimator using
IPW with calibration constraints.

``` r
result_dr <- nonprob(
  selection = ~ x2,
  outcome = y1 ~ x1 + x2,
  data = subset(population, flag_bd1 == 1),
  svydesign = sample_prob
)
```

Results

``` r
summary(result_dr)
#> 
#> Call:
#> nonprob(data = subset(population, flag_bd1 == 1), selection = ~x2, 
#>     outcome = y1 ~ x1 + x2, svydesign = sample_prob)
#> 
#> -------------------------
#> Estimated population mean: 2.95 with overall std.err of: 0.04195
#> And std.err for nonprobability and probability samples being respectively:
#> 0.000783 and 0.04195
#> 
#> 95% Confidence inverval for popualtion mean:
#>    lower_bound upper_bound
#> y1    2.867789     3.03224
#> 
#> 
#> Based on: Doubly-Robust method
#> For a population of estimate size: 1025063
#> Obtained on a nonprobability sample of size: 693011
#> With an auxiliary probability sample of size: 1000
#> -------------------------
#> 
#> Regression coefficients:
#> -----------------------
#> For glm regression on outcome variable:
#>             Estimate Std. Error z value P(>|z|)    
#> (Intercept) 0.996282   0.002139   465.8  <2e-16 ***
#> x1          1.001931   0.001200   835.3  <2e-16 ***
#> x2          0.999125   0.001098   910.2  <2e-16 ***
#> ---
#> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
#> 
#> -----------------------
#> For glm regression on selection variable:
#>              Estimate Std. Error z value P(>|z|)    
#> (Intercept) -0.498997   0.003702  -134.8  <2e-16 ***
#> x2           1.885629   0.005303   355.6  <2e-16 ***
#> -------------------------
#> 
#> Weights:
#>    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#>   1.000   1.071   1.313   1.479   1.798   2.647 
#> -------------------------
#> 
#> Residuals:
#>     Min.  1st Qu.   Median     Mean  3rd Qu.     Max. 
#> -0.99999  0.06603  0.23778  0.26046  0.44358  0.62222 
#> 
#> AIC: 1010622
#> BIC: 1010645
#> Log-Likelihood: -505309 on 694009 Degrees of freedom
```

Mass imputation estimator

``` r
result_mi <- nonprob(
  outcome = y1 ~ x1 + x2,
  data = subset(population, flag_bd1 == 1),
  svydesign = sample_prob
)
```

Results

``` r
summary(result_mi)
#> 
#> Call:
#> nonprob(data = subset(population, flag_bd1 == 1), outcome = y1 ~ 
#>     x1 + x2, svydesign = sample_prob)
#> 
#> -------------------------
#> Estimated population mean: 2.95 with overall std.err of: 0.04203
#> And std.err for nonprobability and probability samples being respectively:
#> 0.001227 and 0.04201
#> 
#> 95% Confidence inverval for popualtion mean:
#>    lower_bound upper_bound
#> y1    2.867433    3.032186
#> 
#> 
#> Based on: Mass Imputation method
#> For a population of estimate size: 1e+06
#> Obtained on a nonprobability sample of size: 693011
#> With an auxiliary probability sample of size: 1000
#> -------------------------
#> 
#> Regression coefficients:
#> -----------------------
#> For glm regression on outcome variable:
#>             Estimate Std. Error z value P(>|z|)    
#> (Intercept) 0.996282   0.002139   465.8  <2e-16 ***
#> x1          1.001931   0.001200   835.3  <2e-16 ***
#> x2          0.999125   0.001098   910.2  <2e-16 ***
#> ---
#> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
#> -------------------------
```

Inverse probability weighting estimator

``` r
result_ipw <- nonprob(
  selection = ~ x2,
  target = ~y1,
  data = subset(population, flag_bd1 == 1),
  svydesign = sample_prob)
```

Results

``` r
summary(result_ipw)
#> 
#> Call:
#> nonprob(data = subset(population, flag_bd1 == 1), selection = ~x2, 
#>     target = ~y1, svydesign = sample_prob)
#> 
#> -------------------------
#> Estimated population mean: 2.925 with overall std.err of: 0.05
#> And std.err for nonprobability and probability samples being respectively:
#> 0.001586 and 0.04997
#> 
#> 95% Confidence inverval for popualtion mean:
#>    lower_bound upper_bound
#> y1     2.82679    3.022776
#> 
#> 
#> Based on: Inverse probability weighted method
#> For a population of estimate size: 1025063
#> Obtained on a nonprobability sample of size: 693011
#> With an auxiliary probability sample of size: 1000
#> -------------------------
#> 
#> Regression coefficients:
#> -----------------------
#> For glm regression on selection variable:
#>              Estimate Std. Error z value P(>|z|)    
#> (Intercept) -0.498997   0.003702  -134.8  <2e-16 ***
#> x2           1.885629   0.005303   355.6  <2e-16 ***
#> ---
#> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
#> -------------------------
#> 
#> Weights:
#>    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#>   1.000   1.071   1.313   1.479   1.798   2.647 
#> -------------------------
#> 
#> Residuals:
#>     Min.  1st Qu.   Median     Mean  3rd Qu.     Max. 
#> -0.99999  0.06603  0.23778  0.26046  0.44358  0.62222 
#> 
#> AIC: 1010622
#> BIC: 1010645
#> Log-Likelihood: -505309 on 694009 Degrees of freedom
```

## Funding

Work on this package is supported by the National Science Centre, OPUS
22 grant no. 2020/39/B/HS4/00941.

## References (selected)

<div id="refs" class="references csl-bib-body hanging-indent">

<div id="ref-chen2020" class="csl-entry">

Chen, Yilin, Pengfei Li, and Changbao Wu. 2020. “Doubly Robust Inference
With Nonprobability Survey Samples.” *Journal of the American
Statistical Association* 115 (532): 2011–21.
<https://doi.org/10.1080/01621459.2019.1677241>.

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
