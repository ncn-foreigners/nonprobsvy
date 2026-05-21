# Inference with Non-Probability Survey Samples

`nonprob` function provides an access to the various methods for
inference based on non-probability surveys (including big data). The
function allows to estimate the population mean based on the access to a
reference probability sample (via the `survey` package), as well as
totals or means of covariates.

The package implements state-of-the-art approaches recently proposed in
the literature: Chen et al. (2020), Yang et al. (2020), Wu (2022) and
uses the [Lumley 2004](https://CRAN.R-project.org/package=survey)
`survey` package for inference (if a reference probability sample is
provided).

It provides various inverse probability weighting (e.g. with calibration
constraints), mass imputation (e.g. nearest neighbour, predictive mean
matching) and doubly robust estimators (e.g. that take into account
minimisation of the asymptotic bias of the population mean estimators).

The package uses the `survey` package functionality when a probability
sample is available.

All optional parameters are set to `NULL`. The obligatory ones include
`data` as well as one of the following three: `selection`, `outcome`, or
`target` – depending on which method has been selected. In the case of
`outcome` and `target` multiple \\y\\ variables can be specified.

## Usage

``` r
nonprob(
  data,
  selection = NULL,
  outcome = NULL,
  target = NULL,
  svydesign = NULL,
  pop_totals = NULL,
  pop_means = NULL,
  pop_size = NULL,
  method_selection = c("logit", "cloglog", "probit"),
  method_outcome = c("glm", "nn", "pmm", "npar"),
  family_outcome = c("gaussian", "binomial", "poisson"),
  subset = NULL,
  strata = NULL,
  case_weights = NULL,
  na_action = na.omit,
  control_selection = control_sel(),
  control_outcome = control_out(),
  control_inference = control_inf(),
  start_selection = NULL,
  start_outcome = NULL,
  verbose = FALSE,
  se = TRUE,
  ...
)
```

## Arguments

- data:

  a `data.frame` with dataset containing the non-probability sample

- selection:

  a `formula` (default `NULL`) for the selection (propensity) score
  model

- outcome:

  a `formula` (default `NULL`) for the outcome (target) model

- target:

  a `formula` (default `NULL`) with target variable(s). We allow
  multiple target variables (e.g. `~y1 + y2 + y3`)

- svydesign:

  an optional `svydesign2` class object containing a probability sample
  and design weights. If finite population correction should affect
  survey-side variance estimates, include the fpc in this object.

- pop_totals:

  an optional `named vector` with population totals of the covariates

- pop_means:

  an optional `named vector` with population means of the covariates;
  `pop_size` must be supplied when `pop_means` is used

- pop_size:

  an optional `double` value with population size. If omitted when a
  probability sample is supplied, the survey-weight denominator defaults
  to `sum(weights(svydesign))`. Supplying `pop_size` fixes the
  population-size denominator used by known-\\N\\ estimators such as
  IPW-MLE; it does not add or modify finite population correction in
  `svydesign`. For IPW-MLE without a fixed `pop_size`, `pop_totals`, or
  `pop_means`, the point estimator is Hajek-type, with the estimated
  IPW-weight total as its denominator. For IPW-GEE with `svydesign`, the
  point estimator uses the survey-weight denominator
  `sum(weights(svydesign))`; if a supplied `pop_size` differs from that
  denominator, a warning is issued.

- method_selection:

  a `character` (default `logit`) indicating the method for the
  propensity score link function.

- method_outcome:

  a `character` (default `glm`) indicating the method for the outcome
  model.

- family_outcome:

  a `character` (default `gaussian`) describing the error distribution
  and the link function to be used in the model. Currently supports:
  `gaussian` with the identity link, `poisson` and `binomial`.

- subset:

  an optional `vector` specifying a subset of observations to be used in
  the fitting process

- strata:

  an optional `vector` specifying strata (not yet supported, for further
  development)

- case_weights:

  an optional positive, finite `numeric vector` of prior weights to be
  used in the fitting process. It is assumed that this vector contains
  frequency or analytic weights (i.e. rows of the `data` argument are
  repeated according to the values of the `case_weights` argument), not
  probability/design weights.

- na_action:

  a function which indicates what should happen when the data contain
  `NAs` (default `na.omit` and it is the only method currently
  supported)

- control_selection:

  a `list` (default
  [`control_sel()`](https://ncn-foreigners.github.io/nonprobsvy/reference/control_sel.md)
  result) indicating parameters to be used when fitting the selection
  model for propensity scores. To change the parameters one should use
  the
  [`control_sel()`](https://ncn-foreigners.github.io/nonprobsvy/reference/control_sel.md)
  function

- control_outcome:

  a `list` (default
  [`control_out()`](https://ncn-foreigners.github.io/nonprobsvy/reference/control_out.md)
  result) indicating parameters to be used when fitting the model for
  the outcome variable. To change the parameters one should use the
  [`control_out()`](https://ncn-foreigners.github.io/nonprobsvy/reference/control_out.md)
  function

- control_inference:

  a `list` (default
  [`control_inf()`](https://ncn-foreigners.github.io/nonprobsvy/reference/control_inf.md)
  result) indicating parameters to be used for inference based on
  probability and non-probability samples. To change the parameters one
  should use the
  [`control_inf()`](https://ncn-foreigners.github.io/nonprobsvy/reference/control_inf.md)
  function

- start_selection:

  an optional `vector` with starting values for the parameters of the
  selection equation

- start_outcome:

  an optional `vector` with starting values for the parameters of the
  outcome equation

- verbose:

  a logical value (default `FALSE`) whether detailed information on the
  fitting should be presented

- se:

  Logical value (default `TRUE`) indicating whether to calculate and
  return standard error of estimated mean

- ...:

  Additional, optional arguments (not yet supported)

## Value

Returns an object of the `nonprob` class (it is actually a `list`) which
contains the following elements:  

- `call` – the call of the `nonprob` function

- `data` – a `data.frame` passed from the `nonprob` function `data`
  argument

- `X` – a `model.matrix` containing data from probability (first
  \\n\_{S\_{\text{P}}}\\ rows) and non-probability samples (next
  \\n\_{S\_{\text{P}}}\\ rows) if specified at a function call

- `y` – a `list` of vector of outcome variables if specified at a
  function call

- `R` – a `numeric vector` indicating whether a unit belongs to the
  probability (0) or non-probability (1) units in the matrix X

- `ps_scores` – a `numeric vector` of estimated propensity scores for
  probability and non-probability sample

- `case_weights` – a `vector` of case weights for non-probability sample
  based on the call

- `ipw_weights` – a `vector` of inverse probability weights for
  non-probability sample (if applicable)

- `boot_ipw_weights` – a `matrix` of bootstrap inverse probability
  weights for the non-probability sample when bootstrap variance
  estimation is used and bootstrap results are kept (if applicable)

- `control` – a `list` of control functions based on the call

- `output` – a `data.frame` with the estimated means and standard errors
  for the variables specified in the `target` or `outcome` arguments

- `SE` – a `data.frame` with standard error of the estimator of the
  population mean, divided into errors from probability and
  non-probability samples (if applicable)

- `confidence_interval` – a `data.frame` with confidence interval of
  population mean estimator

- `nonprob_size` – a scalar `numeric vector` denoting the size of
  non-probability sample

- `prob_size` – a scalar `numeric vector` denoting the size of
  probability sample

- `pop_size` – a scalar `numeric vector` estimated population size
  derived from estimated weights (non-probability sample) or known
  design weights (probability sample)

- `pop_size_fixed` – a `logical` value whether the population size was
  fixed (known) or estimated (unknown)

- `ipw_estimator` – a `character` value with the IPW point-estimator
  family (`"ht"` or `"hajek"`, if applicable)

- `ipw_denominator` – a scalar `numeric vector` with the denominator
  used for the IPW point estimator (if applicable)

- `ipw_denominator_source` – a `character` value describing the source
  of the IPW point-estimator denominator (if applicable)

- `pop_totals` – a `numeric vector` with the total values of the
  auxiliary variables derived from a probability sample or based on the
  call

- `pop_means` – a `numeric vector` with the mean values of the auxiliary
  variables derived from a probability sample or based on the call

- `outcome` – a `list` containing information about the fitting of the
  mass imputation model. Structure of the object is based on the
  `method_outcome` and `family_outcome` arguments which point to
  specific methods as defined by functions `method_*` (if specified in
  the call)

- `selection` – a `list` containing information about the fitting of the
  propensity score model. Structure of the object is based on the
  `method_selection` argument which point to specific methods as defined
  by functions `method_ps` (if specified in the call)

- `boot_sample` – a `matrix` with bootstrap estimates of the target
  variable(s) (if specified)

- `svydesign` – a `svydesign2` object (if specified)

- `ys_rand_pred` – a `list` of predicted values for the target
  variable(s) for the probability sample (for the MI and DR estimator)

- `ys_nons_pred` – a `list` of predicted values for the target
  variable(s) for the non-probability sample (for the MI and DR
  estimator)

- `ys_resid` – a `list` of residuals for the target variable(s) for the
  non-probability sample (for the MI and DR estimator)

- `estimator` – a `character vector` with information what type of
  estimator was selected (one of `c("ipw", "mi", "dr")`).

- `selection_formula` – a `formula` based on the `selection` argument
  (if specified)

- `estimator_method` – a `character vector` with information on the
  detailed method applied (for the `print` method)

## Details

Let \\y\\ be the response variable for which we want to estimate the
population mean, given by \$\$\mu\_{y} = \frac{1}{N} \sum\_{i=1}^N
y\_{i}.\$\$ For this purpose we consider data integration with the
following structure. Let \\S\_{\text{NP}}\\ be the non-probability
sample with the design matrix of covariates as \$\$
\boldsymbol{X}\_{\text{NP}} = \begin{bmatrix} x\_{11} & x\_{12} & \cdots
& x\_{1p} \cr x\_{21} & x\_{22} & \cdots & x\_{2p} \cr \vdots & \vdots &
\ddots & \vdots \cr x\_{n\_{\text{NP}1}} & x\_{n\_{\text{NP}2}} & \cdots
& x\_{n\_{\text{NP}p}} \cr \end{bmatrix}, \$\$ and vector of outcome
variable \$\$ \boldsymbol{y} = \begin{bmatrix} y\_{1} \cr y\_{2} \cr
\vdots \cr y\_{n\_{\text{NP}}} \end{bmatrix}. \$\$ On the other hand,
let \\S\_{\text{P}}\\ be the probability sample with design matrix of
covariates be \$\$ \boldsymbol{X}\_{\text{P}} = \begin{bmatrix} x\_{11}
& x\_{12} & \cdots & x\_{1p} \cr x\_{21} & x\_{22} & \cdots & x\_{2p}
\cr \vdots & \vdots & \ddots & \vdots \cr x\_{n\_{\text{P}1}} &
x\_{n\_{\text{P}2}} & \cdots & x\_{n\_{\text{P}p}}\cr \end{bmatrix}.
\$\$ Instead of a sample of units we can consider a vector of population
sums in the form of \\\tau_x = (\sum\_{i \in
\mathcal{U}}\boldsymbol{x}\_{i1}, \sum\_{i \in
\mathcal{U}}\boldsymbol{x}\_{i2}, ..., \sum\_{i \in
\mathcal{U}}\boldsymbol{x}\_{ip})\\ or means \\\frac{\tau_x}{N}\\, where
\\\mathcal{U}\\ refers to a finite population. Note that we do not
assume access to the response variable for \\S\_{\text{P}}\\. The
implemented estimators assume that outcome values are observed for the
non-probability sample \\S\_{\text{NP}}\\; outcome values observed in
the probability sample \\S\_{\text{P}}\\ are not used. Linked overlap
handling for units appearing in both samples is not currently
implemented;
[`control_sel()`](https://ncn-foreigners.github.io/nonprobsvy/reference/control_sel.md)
arguments `dependence` and `key` are placeholders for future
development.

Supported outcome types depend on the estimator family:

- IPW: numeric targets whose population mean is meaningful, including
  continuous, count, and 0/1 binary variables. No outcome model is
  fitted.

- Mass imputation with `method_outcome = "glm"`: continuous, count, or
  binary variables through `family_outcome = "gaussian"`, `"poisson"`,
  or `"binomial"`.

- Mass imputation with `method_outcome = "nn"`, `"pmm"`, or `"npar"`:
  numeric targets; categorical, ordinal, survival, and other structured
  outcomes are not supported.

- Doubly robust: GLM outcome models only; use
  `family_outcome = "gaussian"`, `"poisson"`, or `"binomial"`.

In general we make the following assumptions:

1.  The selection indicator of belonging to non-probability sample
    \\I\_{\text{NP}, i}\\ and the response variable \\y_i\\ are
    independent given the set of covariates \\\boldsymbol{x}\_i\\.

2.  All units have a non-zero propensity score, i.e., \\\pi\_{\text{NP},
    i} \> 0\\ for all i.

3.  The indicator variables \\I\_{\text{NP}, i}\\ and \\I\_{\text{NP},
    j}\\ are independent for given \\\boldsymbol{x}\_i\\ and
    \\\boldsymbol{x}\_j\\ for \\i \neq j\\.

There are three possible approaches to the problem of estimating
population mean using non-probability samples:

1.  Inverse probability weighting – the main drawback of non-probability
    sampling is the unknown selection mechanism for a unit to be
    included in the sample. This is why we talk about the so-called
    "biased sample" problem. The inverse probability approach is based
    on the assumption that a reference probability sample is available
    and therefore we can estimate the propensity score of the selection
    mechanism. With inverse probability weights \\\hat{d}\_{\text{NP},
    i} = 1 / \hat{\pi}\_{\text{NP}, i}\\, the package supports two IPW
    point-estimator families. The Horvitz-Thompson-type estimator uses
    an external denominator \\N_0\\, \$\$\hat{\mu}\_{IPW,HT} =
    \frac{1}{N_0}\sum\_{i \in S\_{\text{NP}}} w_i \hat{d}\_{\text{NP},
    i} y_i,\$\$ where \\w_i\\ are optional `case_weights`. The
    Hajek-type estimator uses the estimated IPW total as the
    denominator, \$\$\hat{\mu}\_{IPW,H} = \frac{\sum\_{i \in
    S\_{\text{NP}}} w_i \hat{d}\_{\text{NP}, i} y_i}{\sum\_{i \in
    S\_{\text{NP}}} w_i \hat{d}\_{\text{NP}, i}}.\$\$ For IPW-MLE,
    omitting a fixed population size (`pop_size`, `pop_totals`, or
    `pop_means` with `pop_size`) gives the Hajek-type estimator.
    Supplying a fixed population size or population totals gives the
    Horvitz-Thompson-type estimator. For IPW-GEE with a reference
    probability sample, the denominator is `sum(weights(svydesign))`;
    `gee_h_fun` changes the selection-model estimating equation, not the
    point-estimator family. For IPW-GEE with population totals or means,
    the denominator comes from those totals. For this purpose several
    estimation methods can be considered. The first approach is maximum
    likelihood estimation with a corrected log-likelihood function,
    which is given by the following formula \$\$
    \ell^{\*}(\boldsymbol{\theta}) = \sum\_{i \in S\_{\text{NP}}}\log
    \left\lbrace \frac{\pi(\boldsymbol{x}\_{i},
    \boldsymbol{\theta})}{1 -
    \pi(\boldsymbol{x}\_{i},\boldsymbol{\theta})}\right\rbrace +
    \sum\_{i \in S\_{\text{P}}}d\_{\text{P}, i}\log \left\lbrace 1 -
    \pi({\boldsymbol{x}\_{i},\boldsymbol{\theta})}\right\rbrace.\$\$ In
    the literature, the main approach to modelling propensity scores is
    based on the logit link function. However, we extend the propensity
    score model with the additional link functions such as cloglog and
    probit. The pseudo-score equations derived from ML methods can be
    replaced by the idea of generalised estimating equations with
    calibration constraints defined by equations. \$\$
    \mathbf{U}(\boldsymbol{\theta})=\sum\_{i \in S\_{\text{NP}}}
    \mathbf{h}\left(\mathbf{x}\_i, \boldsymbol{\theta}\right)-\sum\_{i
    \in S\_{\text{P}}} d\_{\text{P}, i} \pi\left(\mathbf{x}\_i,
    \boldsymbol{\theta}\right) \mathbf{h}\left(\mathbf{x}\_i,
    \boldsymbol{\theta}\right).\$\$ Notice that for \\
    \mathbf{h}\left(\mathbf{x}\_i, \boldsymbol{\theta}\right) =
    \frac{\boldsymbol{x}}{\pi(\boldsymbol{x}, \boldsymbol{\theta})}\\ We
    do not need a probability sample and can use a vector of population
    totals/means.

2.  Mass imputation – This method is based on a framework where imputed
    values of outcome variables are created for the entire probability
    sample. In this case, we treat the large sample as a training data
    set that is used to build an imputation model. Using the imputed
    values for the probability sample and the (known) design weights, we
    can build a population mean estimator of the form:
    \$\$\hat{\mu}\_{MI} = \frac{1}{N\_{\text{P}}}\sum\_{i \in
    S\_{\text{P}}} d\_{\text{P}, i} \hat{y}\_i.\$\$ It opens the door to
    a very flexible method for imputation models. The package uses
    generalized linear models from
    [`stats::glm()`](https://rdrr.io/r/stats/glm.html), the nearest
    neighbour algorithm using
    [`RANN::nn2()`](https://jefferislab.github.io/RANN/reference/nn2.html)
    and predictive mean matching.

3.  Doubly robust estimation – The IPW and MI estimators are sensitive
    to misspecified models for the propensity score and outcome
    variable, respectively. To this end, so-called doubly robust methods
    are presented that take these problems into account. It is a simple
    idea to combine propensity score and imputation models during
    inference, leading to the following estimator \$\$\hat{\mu}\_{DR} =
    \frac{1}{N\_{\text{NP}}}\sum\_{i \in S\_{\text{NP}}}
    \hat{d}\_{\text{NP}, i} (y_i - \hat{y}\_i) +
    \frac{1}{N\_{\text{P}}}\sum\_{i \in S\_{\text{P}}} d\_{\text{P}, i}
    \hat{y}\_i.\$\$ In addition, an approach based directly on bias
    minimisation has been implemented. The following formula \$\$
    \begin{aligned} bias(\hat{\mu}\_{DR}) = & \mathbb{E}
    (\hat{\mu}\_{DR} - \mu) \cr = & \mathbb{E} \left\lbrace \frac{1}{N}
    \sum\_{i=1}^N (\frac{I\_{\text{NP}, i}}{\pi\_{\text{NP}, i}
    (\boldsymbol{x}\_i^{\mathrm{T}} \boldsymbol{\theta})} - 1 ) (y_i -
    \operatorname{m}(\boldsymbol{x}\_i^{\mathrm{T}} \boldsymbol{\beta}))
    \right\rbrace \cr + & \mathbb{E} \left\lbrace \frac{1}{N}
    \sum\_{i=1}^N (I\_{\text{P}, i} d\_{\text{P}, i} - 1)
    \operatorname{m}( \boldsymbol{x}\_i^{\mathrm{T}} \boldsymbol{\beta})
    \right\rbrace, \end{aligned} \$\$ lead us to system of equations
    \$\$ \begin{aligned} J(\theta, \beta) = \left\lbrace
    \begin{array}{c} J_1(\theta, \beta) \cr J_2(\theta, \beta)
    \end{array}\right\rbrace = \left\lbrace \begin{array}{c}
    \sum\_{i=1}^N I\_{\text{NP}, i}\\ \left\lbrace
    \frac{1}{\pi(\boldsymbol{x}\_i, \boldsymbol{\theta})}-1
    \right\rbrace \left\lbrace y_i-m(\boldsymbol{x}\_i,
    \boldsymbol{\beta}) \right\rbrace \boldsymbol{x}\_i \cr
    \sum\_{i=1}^N \frac{I\_{\text{NP}, i}}{\pi(\boldsymbol{x}\_i,
    \boldsymbol{\theta})} \frac{\partial m(\boldsymbol{x}\_i,
    \boldsymbol{\beta})}{\partial \boldsymbol{\beta}} - \sum\_{i \in
    \mathcal{S}\_{\text{P}}} d\_{\text{P}, i} \frac{\partial
    m(\boldsymbol{x}\_i, \boldsymbol{\beta})}{\partial
    \boldsymbol{\beta}} \end{array} \right\rbrace, \end{aligned} \$\$
    where \\m\left(\boldsymbol{x}\_{i}, \boldsymbol{\beta}\right)\\ is a
    mass imputation (regression) model for the outcome variable and
    propensity scores \\\pi\_{\text{NP}, i}\\ are estimated using a
    `logit` function for the model. As with the `MLE` and `GEE`
    approaches we have extended this method to `cloglog` and `probit`
    links.

As it is not straightforward to calculate the variances of these
estimators, asymptotic equivalents of the variances derived using the
Taylor approximation have been proposed in the literature. Details can
be found [here](https://ncn-foreigners.ue.poznan.pl/nonprobsvy-book/).
In addition, the bootstrap approach can be used for variance estimation.

The doubly robust analytic (Taylor-linearisation) variance is derived
under the logistic propensity model (Chen, Li & Wu 2020, Theorem 2). For
the `probit` and `cloglog` selection links the probability-sample
(design) variance component is a conservative approximation: it tends to
over-estimate the standard error and can be numerically unstable when
fitted propensities approach 1. For doubly robust inference with those
links, `control_inf(var_method = "bootstrap")` is recommended.

The function also allows variables selection using known methods that
have been implemented to handle the integration of probability and
non-probability sampling. In the presence of high-dimensional data,
variable selection is important, because it can reduce the variability
in the estimate that results from using irrelevant variables to build
the model. Let \\\operatorname{U}\left( \boldsymbol{\theta},
\boldsymbol{\beta} \right)\\ be the joint estimating function for
\\\left( \boldsymbol{\theta}, \boldsymbol{\beta} \right)\\. We define
the penalized estimating functions as \$\$ \operatorname{U}^p
\left(\boldsymbol{\theta}, \boldsymbol{\beta}\right) =
\operatorname{U}\left(\boldsymbol{\theta}, \boldsymbol{\beta}\right) -
\left\lbrace \begin{array}{c}
q\_{\lambda\_\theta}(\|\boldsymbol{\theta}\|)
\operatorname{sgn}(\boldsymbol{\theta}) \\
q\_{\lambda\_\beta}(\|\boldsymbol{\beta}\|)
\operatorname{sgn}(\boldsymbol{\beta}) \end{array} \right\rbrace, \$\$
where \\\lambda\_{\theta}\\ and \\q\_{\lambda\_{\beta}}\\ are some
smooth functions. We let \\q\_{\lambda} \left(x\right) = \frac{\partial
p\_{\lambda}}{\partial x}\\, where \\p\_{\lambda}\\ is some penalization
function. Details of penalization functions and techniques for solving
this type of equation can be found
[here](https://ncn-foreigners.ue.poznan.pl/nonprobsvy-book/). To use the
variable selection model, set the `vars_selection` parameter in the
[`control_inf()`](https://ncn-foreigners.github.io/nonprobsvy/reference/control_inf.md)
function to `TRUE`. In addition, in the other control functions such as
[`control_sel()`](https://ncn-foreigners.github.io/nonprobsvy/reference/control_sel.md)
and
[`control_out()`](https://ncn-foreigners.github.io/nonprobsvy/reference/control_out.md)
you can set parameters for the selection of the relevant variables, such
as the number of folds during cross-validation algorithm or the lambda
value for penalizations. Details can be found in the documentation of
the control functions for `nonprob`.

## References

Kim JK, Park S, Chen Y, Wu C. Combining non-probability and probability
survey samples through mass imputation. J R Stat Soc Series A.
2021;184:941– 963.

Shu Yang, Jae Kwang Kim, Rui Song. Doubly robust inference when
combining probability and non-probability samples with high dimensional
data. J. R. Statist. Soc. B (2020)

Yilin Chen , Pengfei Li & Changbao Wu (2020) Doubly Robust Inference
With Nonprobability Survey Samples, Journal of the American Statistical
Association, 115:532, 2011-2021

Shu Yang, Jae Kwang Kim and Youngdeok Hwang Integration of data from
probability surveys and big found data for finite population inference
using mass imputation. Survey Methodology, June 2021 29 Vol. 47, No. 1,
pp. 29-58

## See also

[`stats::optim()`](https://rdrr.io/r/stats/optim.html) – For more
information on the `optim` function used in the `optim` method of
propensity score model fitting.

[`maxLik::maxLik()`](https://rdrr.io/pkg/maxLik/man/maxLik.html) – For
more information on the `maxLik` function used in `maxLik` method of
propensity score model fitting.

[`ncvreg::cv.ncvreg()`](https://pbreheny.github.io/ncvreg/reference/cv.ncvreg.html)
– For more information on the `cv.ncvreg` function used in variable
selection for the outcome model.

[`nleqslv::nleqslv()`](https://bertcarnell.github.io/nleqslv/reference/nleqslv.html)
– For more information on the `nleqslv` function used in estimation
process of the bias minimization approach.

[`stats::glm()`](https://rdrr.io/r/stats/glm.html) – For more
information about the generalised linear models used during mass
imputation process.

[`RANN::nn2()`](https://jefferislab.github.io/RANN/reference/nn2.html) –
For more information about the nearest neighbour algorithm used during
mass imputation process.

[`control_sel()`](https://ncn-foreigners.github.io/nonprobsvy/reference/control_sel.md)
– For the control parameters related to selection model.

[`control_out()`](https://ncn-foreigners.github.io/nonprobsvy/reference/control_out.md)
– For the control parameters related to outcome model.

[`control_inf()`](https://ncn-foreigners.github.io/nonprobsvy/reference/control_inf.md)
– For the control parameters related to statistical inference.

## Author

Łukasz Chrostowski, Maciej Beręsewicz, Piotr Chlebicki

## Examples

``` r
# \donttest{
# generate data based on Doubly Robust Inference With Non-probability Survey Samples (2021)
# Yilin Chen , Pengfei Li & Changbao Wu
set.seed(123)
# sizes of population and probability sample
N <- 20000 # population
n_b <- 1000 # probability
# data
z1 <- rbinom(N, 1, 0.7)
z2 <- runif(N, 0, 2)
z3 <- rexp(N, 1)
z4 <- rchisq(N, 4)

# covariates
x1 <- z1
x2 <- z2 + 0.3 * z2
x3 <- z3 + 0.2 * (z1 + z2)
x4 <- z4 + 0.1 * (z1 + z2 + z3)
epsilon <- rnorm(N)
sigma_30 <- 10.4
sigma_50 <- 5.2
sigma_80 <- 2.4

# response variables
y30 <- 2 + x1 + x2 + x3 + x4 + sigma_30 * epsilon
y50 <- 2 + x1 + x2 + x3 + x4 + sigma_50 * epsilon
y80 <- 2 + x1 + x2 + x3 + x4 + sigma_80 * epsilon

# population
sim_data <- data.frame(y30, y50, y80, x1, x2, x3, x4)
## propensity score model for non-probability sample (sum to 1000)
eta <- -4.461 + 0.1 * x1 + 0.2 * x2 + 0.1 * x3 + 0.2 * x4
rho <- plogis(eta)

# inclusion probabilities for probability sample
z_prob <- x3 + 0.2051
sim_data$p_prob <- n_b* z_prob/sum(z_prob)

# data
sim_data$flag_nonprob <- as.numeric(runif(N) < rho) ## sampling nonprob
sim_data$flag_prob <- as.numeric(runif(n_b) < sim_data$p_prob) ## sampling prob
nonprob_df <- subset(sim_data, flag_nonprob == 1) ## non-probability sample
svyprob <- svydesign(
  ids = ~1, probs = ~p_prob,
  data = subset(sim_data, flag_prob == 1),
  pps = "brewer"
) ## probability sample

## mass imputation estimator
mi_res <- nonprob(
  outcome = y30 + y50 + y80 ~ x1 + x2 + x3 + x4,
  data = nonprob_df,
  svydesign = svyprob
)
mi_res
#> A nonprob object
#>  - estimator type: mass imputation
#>  - method: glm (gaussian)
#>  - auxiliary variables source: survey
#>  - vars selection: false
#>  - variance estimator: analytic
#>  - population size fixed: false
#>  - naive (uncorrected) estimators:
#>    - variable y30: 11.8211
#>    - variable y50: 11.9242
#>    - variable y80: 11.9797
#>  - selected estimators:
#>    - variable y30: 9.4791 (se=0.3927, ci=(8.7094, 10.2488))
#>    - variable y50: 9.5458 (se=0.2241, ci=(9.1067, 9.9850))
#>    - variable y80: 9.5818 (se=0.1517, ci=(9.2845, 9.8791))
## inverse probability weighted estimator
ipw_res <- nonprob(
  selection = ~ x1 + x2 + x3 + x4,
  target = ~y30 + y50 + y80,
  data = nonprob_df,
  svydesign = svyprob
)
ipw_res
#> A nonprob object
#>  - estimator type: inverse probability weighting
#>  - method: logit (mle)
#>  - IPW point estimator: Hajek (denominator: estimated IPW weights = 21240.5738)
#>  - auxiliary variables source: survey
#>  - vars selection: false
#>  - variance estimator: analytic
#>  - population size fixed: false
#>  - naive (uncorrected) estimators:
#>    - variable y30: 11.8211
#>    - variable y50: 11.9242
#>    - variable y80: 11.9797
#>  - selected estimators:
#>    - variable y30: 9.5067 (se=0.4079, ci=(8.7072, 10.3063))
#>    - variable y50: 9.6443 (se=0.2479, ci=(9.1584, 10.1302))
#>    - variable y80: 9.7184 (se=0.1820, ci=(9.3617, 10.0751))
## doubly robust estimator
dr_res <- nonprob(
  outcome = y30 + y50 + y80 ~ x1 + x2 + x3 + x4,
  selection = ~ x1 + x2 + x3 + x4,
  data = nonprob_df,
  svydesign = svyprob
)
dr_res
#> A nonprob object
#>  - estimator type: doubly robust
#>  - method: glm (gaussian)
#>  - auxiliary variables source: survey
#>  - vars selection: false
#>  - variance estimator: analytic
#>  - population size fixed: false
#>  - naive (uncorrected) estimators:
#>    - variable y30: 11.8211
#>    - variable y50: 11.9242
#>    - variable y80: 11.9797
#>  - selected estimators:
#>    - variable y30: 9.3388 (se=0.4548, ci=(8.4473, 10.2302))
#>    - variable y50: 9.4757 (se=0.2528, ci=(8.9802, 9.9712))
#>    - variable y80: 9.5494 (se=0.1618, ci=(9.2322, 9.8666))
# }
```
