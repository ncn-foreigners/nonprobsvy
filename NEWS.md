# nonprobsvy News and Updates

## nonprobsvy 0.2

------------------------------------------------------------------------

### Breaking changes

-   functions `pop.size`, `controlSel`, `controlOut` and `controlInf`
    were renamed to `pop_size`, `control_sel`, `control_out` and
    `control_inf` respectively.
-   function `genSimData` removed completely as it is not used anywhere
    in the package.
-   argument `maxLik_method` renamed to `maxlik_method` in the
    `control_sel` function. 
    
### Features

-   two additional datasets have been included: `jvs` (Job Vacancy
    Survey; a probability sample survey) and `admin` (Central Job Offers
    Database; a non-probability sample survey). The units and auxiliary
    variables have been aligned in a way that allows the data to be
    integrated using the methods implemented in this package.
-   a `nonprobsvycheck` function was added to check the balance in the
    totals of the variables based on the weighted weights between the
    non-probability and probability samples.
-   citation file added.

### Bugfixes

-   basic methods and functions related to variance estimation, weights
    and probability linking methods have been rewritten in a more
    optimal and readable way.

### Other

-   more informative error messages added.
-   documentation improved.
-   switching completely to snake_case.

### Documentation

-   annotation has been added that arguments such as `strata`, `subset`
    and `na_action` are not supported for the time being.

## nonprobsvy 0.1.1

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

## nonprobsvy 0.1.0

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
