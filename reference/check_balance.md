# Checks Variable Balance Between the Probability and Non-Probability Samples

Function compares totals for auxiliary variables specified in the `x`
argument for an `object` that contains either IPW or DR estimator.

## Usage

``` r
check_balance(x, object, dig)
```

## Arguments

- x:

  formula specifying variables to check

- object:

  object of `nonprob` class

- dig:

  number of digits for rounding (default = 2)

## Value

A `list` containing totals for non-probability and probability samples
and their differences

## Examples

``` r

sample_a <- data.frame(y = c(1, 0, 1, 0, 1), x = c(0, 1, 2, 3, 4))
sample_b <- data.frame(x = c(0.5, 1.5, 2.5, 3.5), w = c(4, 4, 4, 4))
sample_b_svy <- svydesign(ids = ~1, weights = ~w, data = sample_b)

ipw_est1 <- nonprob(
  selection = ~x,
  target = ~y,
  svydesign = sample_b_svy,
  data = sample_a,
  method_selection = "logit",
  se = FALSE
)

ipw_est2 <- nonprob(
  selection = ~x,
  target = ~y,
  svydesign = sample_b_svy,
  data = sample_a,
  method_selection = "logit",
  control_selection = control_sel(est_method = "gee", gee_h_fun = 1),
  se = FALSE
)

## check the balance for the standard IPW
check_balance(~x, ipw_est1)
#> $nonprob_totals
#> (Intercept)           x 
#>          16          32 
#> 
#> $prob_totals
#> (Intercept)           x 
#>          16          32 
#> 
#> $balance
#> (Intercept)           x 
#>           0           0 
#> 

## check the balance for the calibrated IPW
check_balance(~x, ipw_est2)
#> $nonprob_totals
#> (Intercept)           x 
#>          16          32 
#> 
#> $prob_totals
#> (Intercept)           x 
#>          16          32 
#> 
#> $balance
#> (Intercept)           x 
#>           0           0 
#> 
```
