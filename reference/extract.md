# Extracts Estimates from the Nonprob Class Object

Returns a `data.frame` of estimated mean(s) or standard error(s)

## Usage

``` r
extract(object, what)
```

## Arguments

- object:

  object of of the `nonprob` class

- what:

  what to extract: all estimates (mean(s), SE(s) and CI(s); `"all"`;
  default), estimated mean(s) (`"mean"`) or their standard error(s)
  (`"se"`)

## Value

a `data.frame` with selected information

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
extract(ipw_est1)
#>         target      mean          SE lower_bound upper_bound
#> 1 single_shift 0.7083229 0.009847656   0.6890218   0.7276239
extract(ipw_est1, "se")
#>         target          SE
#> 1 single_shift 0.009847656
```
