# Summary Statistics for Models of the Nonprob Class

Summarises the `nonprob` class object. The summary depends on the type
of the estimator (i.e. IPW, MI, DR)

## Usage

``` r
# S3 method for class 'nonprob'
summary(object, ...)
```

## Arguments

- object:

  object of the `nonprob` class

- ...:

  Additional optional arguments

## Value

An object of `nonprob_summary` class containing:

- `call` call

- `estimator` type of estimator

- `control` list of controls

- `ipw_weights` estimated IPW weights

- `ipw_weights_total` estimated IPW total (sum)

- `ps_scores_nonprob` estimated propensity scores for non-probability
  sample

- `ps_scores_prob` estimated propensity scores for probability sample

- `case_weights` case weights

- `output` estimated means and standard errors

- `SE` estimated standard errors of V1 and V2

- `confidence_interval` confidence intervals

- `nonprob_size` size of the non-probability sample

- `prob_size` size of the probability sample

- `pop_size` population size

- `pop_size_fixed` whether the population size is treated as fixed

- `ipw_estimator` IPW point-estimator family (`"ht"` or `"hajek"`)

- `ipw_denominator` denominator used for the IPW point estimator

- `ipw_denominator_source` source of the IPW point-estimator denominator

- `no_prob` whether probability sample was provided

- `outcome` model details

- `selection` selection details

- `estimator_method` estimator method

- `selection_formula` selection formula

- `outcome_formula` outcome formula

- `vars_selection` whether variable selection algorithm was applied

- `vars_outcome` variables of the outcome models

- `ys_rand_pred` predicted values for the random sample (if applies)

- `ys_nons_pred` predicted values for the non-probability sample

- `ys_resid` residuals for the non-probability sample

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
summary(ipw_est1)
#> A nonprob_summary object
#>  - call: nonprob(data = admin, selection = ~region + private + nace + 
#>     size, target = ~single_shift, svydesign = jvs_svy, method_selection = "logit")
#>  - estimator type: inverse probability weighting
#>  - nonprob sample size: 9344 (17.7%)
#>  - prob sample size: 6523 (12.3%)
#>  - population size: 52898 (fixed: false)
#>  - IPW point estimator: Hajek (denominator: estimated IPW weights = 52898.1311)
#>  - detailed information about models are stored in list element(s): "selection"
#> ----------------------------------------------------------------
#>  - sum of IPW weights: 52898.13 
#>  - distribution of IPW weights (nonprob sample):
#>    - min: 1.1693; mean: 5.6612; median: 4.3334; max: 49.9504
#>  - distribution of IPW probabilities (nonprob sample):
#>    - min: 0.0189; mean: 0.2894; median: 0.2525; max: 0.8552
#>  - distribution of IPW probabilities (prob sample):
#>    - min: 0.0200; mean: 0.2687; median: 0.2291; max: 0.8552
#> ----------------------------------------------------------------
```
