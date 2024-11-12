# nonprobsvy 0.1.1

-   bugfixes
    -   bug Fix occuring when estimation was based on auxiliary variable, which led to compression of the data from the frame to the vector.
    -   bug Fix related to not passing `maxit` argument from `controlSel` function to internally used `nleqslv` function
    -   bug Fix related to storing `vector` in `model_frame` when predicting `y_hat` in mass imputation `glm` model when X is based in one auxiliary variable only - fix provided converting it to `data.frame` object.
-   features
    -   add information to `summary` about quality of estimation basing on difference between estimated and known total values of auxiliary variables
    -   add estimation of exact standard error for k-nearest neighbor estimator.
    -   add breaking change to `controlOut` function by switching values for `predictive_match` argument. From now on, the `predictive_match = 1` means $\hat{y}-\hat{y}$ in predictive mean matching imputation and `predictive_match = 2` corresponds to $\hat{y}-y$ matching.
    - implement `div` option when variable selection (more in documentation) for doubly robust estimation.
    - add more insights to `nonprob` output such as gradient, hessian and jacobian derived from IPW estimation for `mle` and `gee` methods when `IPW` or `DR` model executed.
    - add estimated inclusion probabilities and its derivatives for probability and non-probability samples to `nonprob` output when `IPW` or `DR` model executed.
    - add `model_frame` matrix data from probability sample used for mass imputation to `nonprob` when `MI` or `DR` model executed.

## nonprobsvy 0.1.0

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
