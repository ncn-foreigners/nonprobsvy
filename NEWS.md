## version 0.1.0

------------------------------------------------------------------------

### Features

-   implemented population mean estimation using doubly robust, inverse probability weighting and mass imputation methods
-   implemented inverse probability weighting models with Maximum Likelihood Estimation and Generalized Estimating Equations methods with `logit`, `complementary log-log` and `probit` link functions.
-   implemented `generalized linear models`, `nearest neighbours` and `predictive mean matching` methods for Mass Imputation
-   implemented bias correction estimators for doubly-robust approach
-   implemented estimation methods when vector of population means/totals is available
-   implemented variables selection with `SCAD`, `LASSO` and `MCP` penalization equations
-   implemented `analytic` and `bootstrap` (with parallel computation - `doSNOW` package) variance for described estimators
-   added control parameters for models
-   added S3 methods for object of `nonprob` class such as
    -   `nobs` for samples size
    -   `pop.size` for population size estimation
    -   `residuals` for residuals of the inverse probability weighting model
    -   `cooks.distance` for identifying influential observations that have a significant impact on the parameter estimates
    -   `hatvalues` for measuring the leverage of individual observations
    -   `logLik` for computing the log-likelihood of the model,
    -   `AIC` (Akaike Information Criterion) for evaluating the model based on the trade-off between goodness of fit and complexity, helping in model selection
    -   `BIC` (Bayesian Information Criterion) for a similar purpose as AIC but with a stronger penalty for model complexity
    -   `confint` for calculating confidence intervals around parameter estimates
    -   `vcov` for obtaining the variance-covariance matrix of the parameter estimates
    -   `deviance` for assessing the goodness of fit of the model

### Unit tests

-   added unit tests for IPW estimators.

### Github repository

-   added automated `R-cmd` check

### Documentation

-   added documentation for `nonprob` function.
