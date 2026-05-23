# Update Method for the Nonprob Object with Changed Arguments or Parameters

The `update` method for the `nonprob` class object that allows to
re-estimate a given model with changed parameters. This is in particular
useful if a user would like to change method or estimate standard errors
if they were not estimated in the first place.

## Usage

``` r
# S3 method for class 'nonprob'
update(object, ..., evaluate = TRUE)
```

## Arguments

- object:

  the `nonprob` class object

- ...:

  arguments passed to the `nonprob` class object

- evaluate:

  If true evaluate the new call else return the call

## Value

returns a `nonprob` object

## Author

Maciej Beręsewicz

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

ipw_est1
#> A nonprob object
#>  - estimator type: IPW (Hajek, denominator: 52898)
#>  - method: logit (mle)
#>  - auxiliary variables source: survey
#>  - vars selection: false
#>  - variance estimator: analytic
#>  - population size fixed: false
#>  - naive (uncorrected) estimator: 0.6605
#>  - selected estimator: 0.7083 (se=NA, ci=(NA, NA))

update(ipw_est1, se = TRUE)
#> A nonprob object
#>  - estimator type: IPW (Hajek, denominator: 52898)
#>  - method: logit (mle)
#>  - auxiliary variables source: survey
#>  - vars selection: false
#>  - variance estimator: analytic
#>  - population size fixed: false
#>  - naive (uncorrected) estimator: 0.6605
#>  - selected estimator: 0.7083 (se=0.0098, ci=(0.6890, 0.7276))
```
