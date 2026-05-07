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

data(admin)
data(jvs)

jvs_svy <- svydesign(ids = ~ 1,  weights = ~ weight,
strata = ~ size + nace + region, data = jvs)

ipw_est1 <- nonprob(selection = ~ region + private + nace + size,
target = ~ single_shift,
svydesign = jvs_svy,
data = admin, method_selection = "logit"
)

ipw_est2 <- nonprob(
selection = ~ region + private + nace + size,
target = ~ single_shift,
svydesign = jvs_svy,
data = admin, method_selection = "logit",
control_selection = control_sel(est_method = "gee", gee_h_fun = 1))

## check the balance for the standard IPW
check_balance(~size+private, ipw_est1)
#> $nonprob_totals
#> (Intercept)       sizeM       sizeS     private 
#>    52898.13    13529.55    31175.21    48214.41 
#> 
#> $prob_totals
#> (Intercept)       sizeM       sizeS     private 
#>       51870       13758       29551       47321 
#> 
#> $balance
#> (Intercept)       sizeM       sizeS     private 
#>     1028.13     -228.45     1624.21      893.41 
#> 

## check the balance for the calibrated IPW
check_balance(~size+private, ipw_est2)
#> $nonprob_totals
#> (Intercept)       sizeM       sizeS     private 
#>       51870       13758       29551       47321 
#> 
#> $prob_totals
#> (Intercept)       sizeM       sizeS     private 
#>       51870       13758       29551       47321 
#> 
#> $balance
#> (Intercept)       sizeM       sizeS     private 
#>           0           0           0           0 
#> 

## check balance for a more complicated example
check_balance(~ I(size=="M") + I(nace == "C"), ipw_est1)
#> $nonprob_totals
#>        (Intercept) I(size == "M")TRUE I(nace == "C")TRUE 
#>          52898.131          13529.550           9113.914 
#> 
#> $prob_totals
#>        (Intercept) I(size == "M")TRUE I(nace == "C")TRUE 
#>              51870              13758               9405 
#> 
#> $balance
#>        (Intercept) I(size == "M")TRUE I(nace == "C")TRUE 
#>            1028.13            -228.45            -291.09 
#> 
```
