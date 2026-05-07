# Returns Population Size (Estimated or Fixed)

Returns population size that is assumed to be

- `fixed` – if it is based on the `pop_size` argument,

- `estimated` – if it is based on the probability survey specified in
  the `svydesign` or based on the estimated propensity scores for the
  non-probability sample.

## Usage

``` r
pop_size(object)
```

## Arguments

- object:

  object returned by the `nonprob` function.

## Value

a scalar returning the value of the population size.

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

## estimated population size based on the non-calibrated IPW (MLE)
pop_size(ipw_est1)
#> pop_size 
#>    51870 

## estimated population size based on the calibrated IPW (GEE)
pop_size(ipw_est2)
#> pop_size 
#>    51870 

```
