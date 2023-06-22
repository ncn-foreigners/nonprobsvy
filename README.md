
<!-- README.md is generated from README.Rmd. Please edit that file -->

# nonprobsvy an R package for modern inference methods based on nonprobability samples.

<!-- badges: start -->

[![R-CMD-check](https://github.com/ncn-foreigners/nonprobsvy/workflows/R-CMD-check/badge.svg)](https://github.com/ncn-foreigners/nonprobsvy/actions)
[![Codecov test
coverage](https://codecov.io/gh/ncn-foreigners/nonprobsvy/branch/main/graph/badge.svg)](https://app.codecov.io/gh/ncn-foreigners/nonprobsvy?branch=main)
<!-- badges: end -->

## Basic information

This package implements mass imputation, inverse probability weighting
and doubly robust estimators based on the following papers:

- Yang, S., Kim, J. K., & Song, R. (2020). Doubly robust inference when
  combining probability and non-probability samples with high
  dimensional data. Journal of the Royal Statistical Society. Series B,
  Statistical Methodology, 82(2), 445.
- Kim, J. K., Park, S., Chen, Y., & Wu, C. (2021). Combining
  non-probability and probability survey samples through mass
  imputation. Journal of the Royal Statistical Society Series A:
  Statistics in Society, 184(3), 941-963.
- Chen, Y., Li, P., & Wu, C. (2020). Doubly robust inference with
  nonprobability survey samples. Journal of the American Statistical
  Association, 115(532), 2011-2021.

## Installation

You can install the development version of `nonprobsvy` package from
[Github](https://github.com/ncn-foreigners/nonprobsvy) with:

``` r
remotes::install_github("ncn-foreigners/nonprobsvy")
```

## Example

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
#> nonprob(selection = ~x2, outcome = y1 ~ x1 + x2, data = subset(population, 
#>     flag_bd1 == 1), svydesign = sample_prob)
#> 
#> -------------------------
#> Estimated population mean: 2.95 with overall std.err of: 0.04195
#> And std.err for nonprobability and probability samples being respectively:
#> 0.000783 and 0.04195
#> 
#> Based on: Doubly-Robust method
#> 
#> 95% Confidence inverval for popualtion mean:
#>        lower_bound upper_bound
#> normal    2.867789     3.03224
#> 
#> For a population of estimate size: 1025111
#> Obtained on a nonprobability sample of size: 693011
#> With an auxiliary probability sample of size: 1000
#> -------------------------
#> 
#> Regression coefficients:
#> -----------------------
#> For glm regression on selection variable:
#>              Estimate Std. Error z value P(>|z|)    
#> (Intercept) -0.499105   0.003702  -134.8  <2e-16 ***
#> x2           1.885532   0.005303   355.6  <2e-16 ***
#> 
#> -----------------------
#> For glm regression on outcome variable:
#>             Estimate Std. Error z value P(>|z|)    
#> (Intercept) 0.996282   0.002139   465.8  <2e-16 ***
#> x1          1.001931   0.001200   835.3  <2e-16 ***
#> x2          0.999125   0.001098   910.2  <2e-16 ***
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
#> -0.99999  0.06604  0.23781  0.26048  0.44361  0.62225 
#> 
#> AIC: 1010626
#> BIC: 1010672
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
#> nonprob(outcome = y1 ~ x1 + x2, data = subset(population, flag_bd1 == 
#>     1), svydesign = sample_prob)
#> 
#> -------------------------
#> Estimated population mean: 2.95 with overall std.err of: 0.04203
#> And std.err for nonprobability and probability samples being respectively:
#> 0.001227 and 0.04201
#> 
#> Based on: Mass Imputation method
#> 
#> 95% Confidence inverval for popualtion mean:
#>        lower_bound upper_bound
#> normal    2.867433    3.032186
#> 
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
#> nonprob(selection = ~x2, target = ~y1, data = subset(population, 
#>     flag_bd1 == 1), svydesign = sample_prob)
#> 
#> -------------------------
#> Estimated population mean: 2.925 with overall std.err of: 0.04999
#> And std.err for nonprobability and probability samples being respectively:
#> 0.001117 and 0.04997
#> 
#> Based on: Inverse probability weighted method
#> 
#> 95% Confidence inverval for popualtion mean:
#>        lower_bound upper_bound
#> normal    2.826789    3.022736
#> 
#> For a population of estimate size: 1025111
#> Obtained on a nonprobability sample of size: 693011
#> With an auxiliary probability sample of size: 1000
#> -------------------------
#> 
#> Regression coefficients:
#> -----------------------
#> For glm regression on selection variable:
#>              Estimate Std. Error z value P(>|z|)    
#> (Intercept) -0.499105   0.003702  -134.8  <2e-16 ***
#> x2           1.885532   0.005303   355.6  <2e-16 ***
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
#> -0.99999  0.06604  0.23781  0.26048  0.44361  0.62225 
#> 
#> AIC: 1010626
#> BIC: 1010672
#> Log-Likelihood: -505309 on 694009 Degrees of freedom
```

## Funding

Work on this package is supported by the the National Science Center,
OPUS 22 grant no. 2020/39/B/HS4/00941.
