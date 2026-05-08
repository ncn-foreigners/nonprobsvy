# Returns Population Size (Estimated or Fixed)

Returns population size that is assumed to be

- `fixed` – if it is based on the `pop_size` argument or population
  totals,

- `estimated` – if it is based on the probability survey specified in
  the `svydesign` or, for Hajek-type IPW-MLE, on the estimated IPW total
  for the non-probability sample.

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

## estimated population size based on the non-calibrated IPW (MLE)
pop_size(ipw_est1)
#> pop_size 
#>       16 

## estimated population size based on the calibrated IPW (GEE)
pop_size(ipw_est2)
#> pop_size 
#>       16 

```
