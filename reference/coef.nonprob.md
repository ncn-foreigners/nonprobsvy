# Returns Coefficients of the Underlying Models

Returns a `list` of coefficients for the selection and the outcome
models

## Usage

``` r
# S3 method for class 'nonprob'
coef(object, ...)
```

## Arguments

- object:

  a `nonprob` class object

- ...:

  other arguments passed to methods (currently not supported)

## Value

a `list` with two entries:

- `"coef_sel"` a matrix of coefficients for the selection equation if
  possible, else NULL

- `"coef_dr"` a matrix of coefficients for the outcome equation(s) if
  possible, else NULL

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

coef(ipw_est1)
#> $coef_sel
#>                    [,1]
#> (Intercept) -0.65277121
#> region04     0.83779549
#> region06     0.19953188
#> region08     0.10479761
#> region10    -0.15756693
#> region12    -0.60987005
#> region14    -0.84150300
#> region16     0.76385615
#> region18     1.17810496
#> region20     0.22251032
#> region22    -0.03753066
#> region24    -0.40670702
#> region26     0.20287478
#> region28     0.57862243
#> region30    -0.61021662
#> region32     0.32742204
#> private      0.05899377
#> naceD.E      0.77273741
#> naceF       -0.37783229
#> naceG       -0.33370352
#> naceH       -0.65174566
#> naceI        0.41179116
#> naceJ       -1.42636886
#> naceK.L      0.06171389
#> naceM       -0.40677813
#> naceN        0.80034162
#> naceO       -0.69354695
#> naceP        1.25095429
#> naceQ        0.30286775
#> naceR.S      0.22227967
#> sizeM       -0.36412124
#> sizeS       -1.02915884
#> 
#> $coef_out
#> NULL
#> 
```
