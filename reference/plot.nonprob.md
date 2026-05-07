# Plots the Estimated Mean(s) and Their Confidence Interval(s)

Simple plotting method that compares the estimated mean(s) and CI(s)
with the naive (uncorrected) estimates.

## Usage

``` r
# S3 method for class 'nonprob'
plot(x, ...)
```

## Arguments

- x:

  the `nonprob` class object

- ...:

  other arguments passed to the plot method (currently not supported)

## Examples

``` r

data(admin)
data(jvs)

jvs_svy <- svydesign(ids = ~ 1,  weights = ~ weight,
strata = ~ size + nace + region, data = jvs)

ipw_est1 <- nonprob(selection = ~ region + private + nace + size,
target = ~ single_shift,
svydesign = jvs_svy,
data = admin, method_selection = "logit")

plot(ipw_est1)

```
