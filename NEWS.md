## version 0.0.1

------------------------------------------------------------------------

### Features

- implemented population mean estimation using doubly robust, inverse probability weighting and mass imputation methods
- implemented inverse probability weighting models with Maximum Likelihood Estimation and Generalized Estimating Equations methods with `logit`, `complementary log-log` and `probit` link functions.
- implemented `generalized linear models` and `nearest neighbours` methods for Mass Imputation
- implemented estimation methods when vector of population means/totals is available
- implemented variables selection with `SCAD`, `LASSO` and `MCP` penalization equations
- implemented `analytic` and `bootstrap` (with parallel computation) variance for described estimators
- added control parameters for models

### Unit tests

- added unit tests for IPW estimators.

### Github repository

- added automated `R-cmd` check

### Documentation

- added documentation for `nonprob` function.
