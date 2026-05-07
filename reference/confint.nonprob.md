# Returns Confidence Intervals for Estimated Mean

A generic function that returns the confidence interval for the
estimated mean. If standard errors have not been estimated, the function
updates the `nonprob` object to obtain standard errors.

## Usage

``` r
# S3 method for class 'nonprob'
confint(object, parm, level = 0.95, ...)
```

## Arguments

- object:

  object of `nonprob` class.

- parm:

  names of parameters for which confidence intervals are to be computed,
  if missing all parameters will be considered.

- level:

  confidence level for intervals.

- ...:

  additional arguments

## Value

returns a `data.frame` with confidence intervals for the target
variables

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

confint(ipw_est1)
#> Calculating standard errors...
#> Error in eval(call, parent.frame()): object 'jvs_svy' not found
```
