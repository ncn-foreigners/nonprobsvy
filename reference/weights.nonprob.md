# Extracts the Inverse Probability (Propensity Score) Weights

A generic function `weights` that returns inverse probability weights
(if present)

## Usage

``` r
# S3 method for class 'nonprob'
weights(object, ...)
```

## Arguments

- object:

  a `nonprob` class object

- ...:

  other arguments passed to methods (currently not supported)

## Value

A vector of weights or a `NULL` extracted from the `nonprob` object i.e.
element `"ipw_weights"`

## Examples

``` r

data(admin)
data(jvs)

jvs_svy <- svydesign(ids = ~ 1,  weights = ~ weight,
strata = ~ size + nace + region, data = jvs)

ipw_est1 <- nonprob(selection = ~ region + private + nace + size,
target = ~ single_shift,
svydesign = jvs_svy,
data = admin, method_selection = "logit", se = FALSE
)

summary(weights(ipw_est1))
#>    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#>   1.169   2.673   4.333   5.661   7.178  49.950 
```
