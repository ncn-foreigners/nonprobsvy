nonprobsvy News and Updates


# nonprobsvy 0.2.0.9001 (development)

------------------------------------------------------------------------
 
+ titles corrected 
+ new S3 method `extract` added which allows to extract results from the `nonprob` object
+ new S3 method `coef` added which allows to obtain the coefficients of underlying models (if possible)
+ fixed CRAN notes (unit tests for the IPW estimator `cloglog`) 
+ removed `sampling` package from suggested package 
+ added simple `plot` method
+ improvements in the linear algebra
+ corrected the `check_balance` error as in 
+ code cleaning

# nonprobsvy 0.2.0

------------------------------------------------------------------------

### Breaking changes

- functions `pop.size`, `controlSel`, `controlOut` and `controlInf`
  were renamed to `pop_size`, `control_sel`, `control_out` and
  `control_inf` respectively.
- function `genSimData` removed completely as it is not used anywhere
  in the package.
- argument `maxLik_method` renamed to `maxlik_method` in the
  `control_sel` function.
- `control_out` function:
  - `predictive_match` renamed to `pmm_match_type` to align with the
    PMM (Predictive Mean Matching) estimator naming convention,
    where all related parameters start with `pmm_`
- `control_sel` function:
  - argument `method` removed as it was not used
  - argument `est_method_sel` renamed to `est_method`
  - argument `h` renamed to `gee_h_fun` to make this more readable
    to the user
  - `start_type` now accepts only `zero` and `mle` (for `gee` models
    only).
- `control_inf` function:
  - `bias_inf` renamed to `vars_combine` and type changed to
    `logical`. `TRUE` if variables (its levels) should be combined
    after variable selection algorithm for the doubly robust
    approach.
  - `pi_ij` -- argument removed as it is not used.
- `nonprobsvy` class renamed to `nonprob` and all related method
  adjusted to this change
- functions `logit_model_nonprobsvy`, `probit_model_nonprobsvy` and
  `cloglog_model_nonprobsvy` removed in the favour of more readable
  `method_ps` function that specifies the propensity score model
- new option `control_inference=control_inf(vars_combine=TRUE)` which
  allows doubly robust estimator to combine variables prior estimation
  i.e. if `selection=~x1+x2` and `y~x1+x3` then the following models
  are fitted `selection=~x1+x2+x3` and `y~x1+x2+x3`. By default we set
  `control_inference=control_inf(vars_combine=FALSE)`. Note that this
  behaviour is assumed independently from variable selection.
- argument `nonprob(weights=NULL)` replaced to `nonprob(case_weights=NULL)`
  to stress that this refer to case weights not sampling or other weights
  in non-probability sample

### Features

- two additional datasets have been included: `jvs` (Job Vacancy
  Survey; a probability sample survey) and `admin` (Central Job Offers
  Database; a non-probability sample survey). The units and auxiliary
  variables have been aligned in a way that allows the data to be
  integrated using the methods implemented in this package.
- a `check_balance` function was added to check the balance in the
  totals of the variables based on the weighted weights between the
  non-probability and probability samples.
- citation file added.
- argument `na_action` with default `na.omit`
- new generic methods added:
  - `weights` -- returns IPW weights
  - `update` -- allows to update the `nonprob` class object
- new functions added and exported:
  - `method_ps` -- for modelling propensity score
  - `method_glm` -- for modelling y using `glm` function
  - `method_nn` -- for the NN method
  - `method_pmm` -- for the PMM method
  - `method_npar` -- for the non-parametric method
- new `print.nonprob`, `summary.nonprob` and `print.nonprob_summary`
  methods

``` r
> result_mi
A nonprob object
 - estimator type: mass imputation
 - method: glm (gaussian)
 - auxiliary variables source: survey
 - vars selection: false
 - variance estimator: analytic
 - population size fixed: false
 - naive (uncorrected) estimators:
   - variable y1: 3.1817
   - variable y2: 1.8087
 - selected estimators:
   - variable y1: 2.9498 (se=0.0420, ci=(2.8674, 3.0322))
   - variable y2: 1.5760 (se=0.0326, ci=(1.5122, 1.6399))
```

number of digits can be changed using `print(x, digits)` as shown below

``` r
> print(result_mi,2)
A nonprob object
 - estimator type: mass imputation
 - method: glm (gaussian)
 - auxiliary variables source: survey
 - vars selection: false
 - variance estimator: analytic
 - population size fixed: false
 - naive (uncorrected) estimators:
   - variable y1: 3.18
   - variable y2: 1.81
 - selected estimators:
   - variable y1: 2.95 (se=0.04, ci=(2.87, 3.03))
   - variable y2: 1.58 (se=0.03, ci=(1.51, 1.64))
```

``` r
> summary(result_mi) |> print(digits=2)
A nonprob_summary object
 - call: nonprob(data = subset(population, flag_bd1 == 1), outcome = y1 + 
    y2 ~ x1 + x2, svydesign = sample_prob)
 - estimator type: mass imputation
 - nonprob sample size: 693011 (69.3%)
 - prob sample size: 1000 (0.1%)
 - population size: 1000000 (fixed: false)
 - detailed information about models are stored in list element(s): "outcome"
----------------------------------------------------------------
 - distribution of outcome residuals:
   - y1: min: -4.79; mean: 0.00; median: 0.00; max: 4.54
   - y2: min: -4.96; mean: -0.00; median: -0.07; max: 12.25
 - distribution of outcome predictions (nonprob sample):
   - y1: min: -2.72; mean: 3.18; median: 3.04; max: 16.28
   - y2: min: -1.55; mean: 1.81; median: 1.58; max: 13.92
 - distribution of outcome predictions (prob sample):
   - y1: min: -0.46; mean: 2.95; median: 2.84; max: 10.31
   - y2: min: -0.58; mean: 1.58; median: 1.39; max: 7.87
----------------------------------------------------------------
```

### Bugfixes

-   basic methods and functions related to variance estimation, weights
    and probability linking methods have been rewritten in a more
    optimal and readable way.

### Other

-   more informative error messages added.
-   documentation improved.
-   switching completely to snake_case.
-   extensive cleaning of the code.
-   more unit-tests added.
-   new dependencies: `formula.tools`

### Documentation

-   annotation has been added that argument `strata` is not supported for the time being.

### Replication materials

-   to verify the quality of the software please refer to the
    replication materials available here:
    <https://github.com/ncn-foreigners/software-tutorials>

# nonprobsvy 0.1.1

------------------------------------------------------------------------

### Bugfixes

-   bug Fix occurring when estimation was based on auxiliary variable,
    which led to compression of the data from the frame to the vector.
-   bug Fix related to not passing `maxit` argument from `controlSel`
    function to internally used `nleqslv` function
-   bug Fix related to storing `vector` in `model_frame` when predicting
    `y_hat` in mass imputation `glm` model when X is based in one
    auxiliary variable only - fix provided converting it to `data.frame`
    object.

### Features

-   added information to `summary` about quality of estimation basing on
    difference between estimated and known total values of auxiliary
    variables
-   added estimation of exact standard error for k-nearest neighbor
    estimator.
-   added breaking change to `controlOut` function by switching values
    for `predictive_match` argument. From now on, the
    `predictive_match = 1` means $\hat{y}-\hat{y}$ in predictive mean
    matching imputation and `predictive_match = 2` corresponds to
    $\hat{y}-y$ matching.
-   implemented `div` option when variable selection (more in
    documentation) for doubly robust estimation.
-   added more insights to `nonprob` output such as gradient, hessian
    and jacobian derived from IPW estimation for `mle` and `gee` methods
    when `IPW` or `DR` model executed.
-   added estimated inclusion probabilities and its derivatives for
    probability and non-probability samples to `nonprob` output when
    `IPW` or `DR` model executed.
-   added `model_frame` matrix data from probability sample used for
    mass imputation to `nonprob` when `MI` or `DR` model executed.

### Unit tests

-   added unit tests for variable selection models and mi estimation
    with vector of population totals available

# nonprobsvy 0.1.0

------------------------------------------------------------------------

### Features

-   implemented population mean estimation using doubly robust, inverse
    probability weighting and mass imputation methods
-   implemented inverse probability weighting models with Maximum
    Likelihood Estimation and Generalized Estimating Equations methods
    with `logit`, `complementary log-log` and `probit` link functions.
-   implemented `generalized linear models`, `nearest neighbours` and
    `predictive mean matching` methods for Mass Imputation
-   implemented bias correction estimators for doubly-robust approach
-   implemented estimation methods when vector of population
    means/totals is available
-   implemented variables selection with `SCAD`, `LASSO` and `MCP`
    penalization equations
-   implemented `analytic` and `bootstrap` (with parallel computation -
    `doSNOW` package) variance for described estimators
-   added control parameters for models
-   added S3 methods for object of `nonprob` class such as
    -   `nobs` for samples size
    -   `pop.size` for population size estimation
    -   `residuals` for residuals of the inverse probability weighting
        model
    -   `cooks.distance` for identifying influential observations that
        have a significant impact on the parameter estimates
    -   `hatvalues` for measuring the leverage of individual
        observations
    -   `logLik` for computing the log-likelihood of the model,
    -   `AIC` (Akaike Information Criterion) for evaluating the model
        based on the trade-off between goodness of fit and complexity,
        helping in model selection
    -   `BIC` (Bayesian Information Criterion) for a similar purpose as
        AIC but with a stronger penalty for model complexity
    -   `confint` for calculating confidence intervals around parameter
        estimates
    -   `vcov` for obtaining the variance-covariance matrix of the
        parameter estimates
    -   `deviance` for assessing the goodness of fit of the model

### Unit tests

-   added unit tests for IPW estimators.

### Github repository

-   added automated `R-cmd` check

### Documentation

-   added documentation for `nonprob` function.
