#' @title Inference with non-probability survey samples
#' @author Łukasz Chrostowski, Maciej Beręsewicz, Piotr Chlebicki
#'
#' @description \code{nonprob} function provides an access to the various methods for inference based on non-probability surveys (including big data). The function allows to estimate the population mean based on the access to a reference probability sample (via the `survey` package),  as well as totals or means of covariates.
#'
#' The package implements state-of-the-art approaches recently proposed in the literature: Chen et al. (2020),
#' Yang et al. (2020), Wu (2022) and uses the [Lumley 2004](https://CRAN.R-project.org/package=survey) `survey` package for inference (if a reference probability sample is provided).
#'
#' It provides various inverse probability weighting (e.g. with calibration constraints), mass imputation (e.g. nearest neighbour, predictive mean matching) and doubly robust estimators (e.g. that take into account minimisation of the asymptotic bias of the population mean estimators).
#'
#' The package uses the `survey` package functionality when a probability sample is available.
#'
#' All optional parameters are set to `NULL`. The obligatory ones include `data` as well as one of the following three:
#'  \code{selection}, \code{outcome}, or \code{target} -- depending on which method has been selected.
#'  In the case of \code{outcome} and \code{target} multiple \eqn{y} variables can be specified.
#'
#' @param data a `data.frame` with dataset containing the non-probability sample
#' @param selection a `formula` (default `NULL`) for the selection (propensity) score model
#' @param outcome a `formula` (default `NULL`) for the outcome (target) model
#' @param target a `formula` (default `NULL`) with target variable(s). We allow multiple target variables (e.g. `~y1 + y2 + y3`)
#' @param svydesign an optional `svydesign2` class object containing a probability sample and design weights
#' @param pop_totals an optional `named vector` with population totals of the covariates
#' @param pop_means an optional `named vector` with population means of the covariates
#' @param pop_size an optional `double` value with population size
#' @param method_selection a `character` (default `logit`) indicating the method for the propensity score link function.
#' @param method_outcome a `character` (default `glm`) indicating the method for the outcome model.
#' @param family_outcome a `character` (default `gaussian`)  describing the error distribution and the link function to be used in the model. Currently supports: `gaussian` with the identity link, `poisson` and `binomial`.
#' @param subset an optional `vector` specifying a subset of observations to be used in the fitting process
#' @param strata an optional `vector` specifying strata (not yet supported, for further development)
#' @param case_weights an optional `vector` of prior weights to be used in the fitting process.
#' It is assumed that this vector contains frequency or analytic weights (i.e. rows of the `data` argument are repeated according to the values of the `case_weights` argument), not probability/design weights.
#' @param na_action a function which indicates what should happen when the data contain `NAs` (default `na.omit` and it is the only method currently supported)
#' @param control_selection a `list` (default `control_sel()` result) indicating parameters to be used when fitting the selection model for propensity scores. To change the parameters one should use the `control_sel()` function
#' @param control_outcome a `list` (default `control_out()` result) indicating parameters to be used when fitting the model for the outcome variable. To change the parameters one should use the `control_out()` function
#' @param control_inference a `list` (default `control_inf()` result) indicating parameters to be used for inference based on probability and non-probability samples. To change the parameters one should use the `control_inf()` function
#' @param start_selection an optional `vector` with starting values for the parameters of the selection equation
#' @param start_outcome an optional `vector` with starting values for the parameters of the outcome equation
#' @param verbose a numerical value (default `TRUE`) whether detailed information on the fitting should be presented
#' @param se Logical value (default `TRUE`) indicating whether to calculate and return standard error of estimated mean.
#' @param ... Additional, optional arguments
#'
#' @details Let \eqn{y} be the response variable for which we want to estimate the population mean,
#' given by \deqn{\mu_{y} = \frac{1}{N} \sum_{i=1}^N y_{i}.} For this purpose we consider data integration
#' with the following structure. Let \eqn{S_A} be the non-probability sample with the design matrix of covariates as
#' \deqn{
#' \boldsymbol{X}_A =
#'   \begin{bmatrix}
#' x_{11} & x_{12} & \cdots & x_{1p} \cr
#' x_{21} & x_{22} & \cdots & x_{2p} \cr
#' \vdots & \vdots & \ddots & \vdots \cr
#' x_{n_{A1}} & x_{n_{A2}} & \cdots & x_{n_{Ap}} \cr
#' \end{bmatrix},
#' }
#' and vector of outcome variable
#' \deqn{
#' \boldsymbol{y} =
#'   \begin{bmatrix}
#' y_{1} \cr
#' y_{2} \cr
#' \vdots \cr
#' y_{n_{A}}
#' \end{bmatrix}.
#' }
#' On the other hand, let \eqn{S_B} be the probability sample with design matrix of covariates be
#' \deqn{
#' \boldsymbol{X}_B =
#'   \begin{bmatrix}
#' x_{11} & x_{12} & \cdots & x_{1p} \cr
#' x_{21} & x_{22} & \cdots & x_{2p} \cr
#' \vdots & \vdots & \ddots & \vdots \cr
#' x_{n_{B1}} & x_{n_{B2}} & \cdots & x_{n_{Bp}}\cr
#' \end{bmatrix}.
#' }
#' Instead of a sample of units we can consider a vector of population sums in the form of \eqn{\tau_x = (\sum_{i \in \mathcal{U}}\boldsymbol{x}_{i1}, \sum_{i \in \mathcal{U}}\boldsymbol{x}_{i2}, ..., \sum_{i \in \mathcal{U}}\boldsymbol{x}_{ip})} or means
#' \eqn{\frac{\tau_x}{N}}, where \eqn{\mathcal{U}} refers to a finite population. Note that we do not assume access to the response variable for \eqn{S_B}.
#' In general we make the following assumptions:
#' 1.  The selection indicator of belonging to non-probability sample \eqn{R_{i}} and the response variable \eqn{y_i} are independent given the set of covariates \eqn{\boldsymbol{x}_i}.
#' 2.  All units have a non-zero propensity score, i.e., \eqn{\pi_{i}^{A} > 0} for all i.
#' 3.  The indicator variables \eqn{R_{i}^{A}} and \eqn{R_{j}^{A}} are independent for given \eqn{\boldsymbol{x}_i} and \eqn{\boldsymbol{x}_j} for \eqn{i \neq j}.
#'
#' There are three possible approaches to the problem of estimating population mean using non-probability samples:
#'
#' 1. Inverse probability weighting -- the main drawback of non-probability sampling is the unknown selection mechanism for a unit to be included in the sample.
#'  This is why we talk about the so-called "biased sample" problem. The inverse probability approach is based on the assumption that a reference probability sample
#'  is available and therefore we can estimate the propensity score of the selection mechanism.
#'  The estimator has the following form:
#'  \deqn{\hat{\mu}_{IPW} = \frac{1}{N^{A}}\sum_{i \in S_{A}} \frac{y_{i}}{\hat{\pi}_{i}^{A}}.}
#'  For this purpose several estimation methods can be considered. The first approach is maximum likelihood estimation with a corrected
#'  log-likelihood function, which is given by the following formula
#'  \deqn{
#'  \ell^{*}(\boldsymbol{\theta}) = \sum_{i \in S_{A}}\log \left\lbrace \frac{\pi(\boldsymbol{x}_{i}, \boldsymbol{\theta})}{1 - \pi(\boldsymbol{x}_{i},\boldsymbol{\theta})}\right\rbrace + \sum_{i \in S_{B}}d_{i}^{B}\log \left\lbrace 1 - \pi({\boldsymbol{x}_{i},\boldsymbol{\theta})}\right\rbrace.}
#'  In the literature, the main approach to modelling propensity scores is based on the logit link function.
#'  However, we extend the propensity score model with the additional link functions such as cloglog and probit.
#'  The pseudo-score equations derived from ML methods can be replaced by the idea of generalised estimating equations
#'  with calibration constraints defined by equations.
#'  \deqn{
#'  \mathbf{U}(\boldsymbol{\theta})=\sum_{i \in S_A} \mathbf{h}\left(\mathbf{x}_i, \boldsymbol{\theta}\right)-\sum_{i \in S_B} d_i^B \pi\left(\mathbf{x}_i, \boldsymbol{\theta}\right) \mathbf{h}\left(\mathbf{x}_i, \boldsymbol{\theta}\right).}
#'  Notice that for \eqn{ \mathbf{h}\left(\mathbf{x}_i, \boldsymbol{\theta}\right) = \frac{\pi(\boldsymbol{x}, \boldsymbol{\theta})}{\boldsymbol{x}}} We do not need a probability sample and can use a vector of population totals/means.
#'
#' 2. Mass imputation -- This method is based on a framework where imputed values of outcome variables are created for the entire probability sample. In this case, we treat the large sample as a training data set that is used to build an imputation model.
#'    Using the imputed values for the probability sample and the (known) design weights,
#'    we can build a population mean estimator of the form:
#'    \deqn{\hat{\mu}_{MI} = \frac{1}{N^B}\sum_{i \in S_{B}} d_{i}^{B} \hat{y}_i.}
#'    It opens the door to a very flexible method for imputation models. The package uses generalized linear models from [stats::glm()],
#'    the nearest neighbour algorithm using [RANN::nn2()] and predictive mean matching.
#'
#' 3. Doubly robust estimation -- The IPW and MI estimators are sensitive to misspecified models for the propensity score and outcome variable, respectively.
#'    To this end, so-called doubly robust methods are presented that take these problems into account.
#'    It is a simple idea to combine propensity score and imputation models during inference, leading to the following estimator
#'    \deqn{\hat{\mu}_{DR} = \frac{1}{N^A}\sum_{i \in S_A} \hat{d}_i^A (y_i - \hat{y}_i) + \frac{1}{N^B}\sum_{i \in S_B} d_i^B \hat{y}_i.}
#'    In addition, an approach based directly on bias minimisation has been implemented. The following formula
#'    \deqn{
#'    \begin{aligned}
#'    bias(\hat{\mu}_{DR}) = & \mathbb{E} (\hat{\mu}_{DR} - \mu) \cr = & \mathbb{E} \left\lbrace \frac{1}{N} \sum_{i=1}^N (\frac{R_i^A}{\pi_i^A (\boldsymbol{x}_i^{\mathrm{T}} \boldsymbol{\theta})}
#'    - 1 ) (y_i - \operatorname{m}(\boldsymbol{x}_i^{\mathrm{T}} \boldsymbol{\beta})) \right\rbrace \cr + & \mathbb{E} \left\lbrace \frac{1}{N} \sum_{i=1}^N (R_i^B d_i^B - 1) \operatorname{m}( \boldsymbol{x}_i^{\mathrm{T}} \boldsymbol{\beta}) \right\rbrace,
#'    \end{aligned}
#'    }
#'    lead us to system of equations
#'    \deqn{
#'    \begin{aligned}
#'    J(\theta, \beta) =
#'    \left\lbrace
#'    \begin{array}{c}
#'                       J_1(\theta, \beta) \cr
#'                       J_2(\theta, \beta)
#'                       \end{array}\right\rbrace = \left\lbrace \begin{array}{c}
#'                                                \sum_{i=1}^N R_i^A\ \left\lbrace \frac{1}{\pi(\boldsymbol{x}_i, \boldsymbol{\theta})}-1 \right\rbrace \left\lbrace y_i-m(\boldsymbol{x}_i, \boldsymbol{\beta}) \right\rbrace \boldsymbol{x}_i \cr
#'                                                \sum_{i=1}^N \frac{R_i^A}{\pi(\boldsymbol{x}_i, \boldsymbol{\theta})} \frac{\partial m(\boldsymbol{x}_i, \boldsymbol{\beta})}{\partial \boldsymbol{\beta}}
#'                                                - \sum_{i \in \mathcal{S}_{\mathrm{B}}} d_i^{\mathrm{B}} \frac{\partial m(\boldsymbol{x}_i, \boldsymbol{\beta})}{\partial \boldsymbol{\beta}}
#'   \end{array} \right\rbrace,
#'   \end{aligned}
#'   }
#'   where \eqn{m\left(\boldsymbol{x}_{i}, \boldsymbol{\beta}\right)} is a mass imputation (regression) model for the outcome variable and
#'   propensity scores \eqn{\pi_i^A} are estimated using a `logit` function for the model. As with the `MLE` and `GEE` approaches we have extended
#'   this method to `cloglog` and `probit` links.
#'
#'   As it is not straightforward to calculate the variances of these estimators, asymptotic equivalents of the variances derived using the Taylor approximation have been proposed in the literature.
#'   Details can be found [here](https://ncn-foreigners.github.io/nonprobsvy-book/).
#'   In addition, the bootstrap approach can be used for variance estimation.
#'
#'   The function also allows variables selection using known methods that have been implemented to handle the integration of probability and non-probability sampling.
#'   In the presence of high-dimensional data, variable selection is important, because it can reduce the variability in the estimate that results from using irrelevant variables to build the model.
#'   Let \eqn{\operatorname{U}\left( \boldsymbol{\theta}, \boldsymbol{\beta} \right)} be the joint estimating function for \eqn{\left( \boldsymbol{\theta}, \boldsymbol{\beta} \right)}. We define the
#'   penalized estimating functions as
#'   \deqn{
#'   \operatorname{U}^p \left(\boldsymbol{\theta}, \boldsymbol{\beta}\right) =
#'   \operatorname{U}\left(\boldsymbol{\theta}, \boldsymbol{\beta}\right) -
#'   \left\lbrace
#'   \begin{array}{c}
#'   q_{\lambda_\theta}(|\boldsymbol{\theta}|) \operatorname{sgn}(\boldsymbol{\theta}) \\
#'   q_{\lambda_\beta}(|\boldsymbol{\beta}|) \operatorname{sgn}(\boldsymbol{\beta})
#'   \end{array}
#'   \right\rbrace,
#'   }
#'   where \eqn{\lambda_{\theta}} and \eqn{q_{\lambda_{\beta}}} are some smooth functions. We let \eqn{q_{\lambda} \left(x\right) = \frac{\partial p_{\lambda}}{\partial x}}, where \eqn{p_{\lambda}} is some penalization function.
#'   Details of penalization functions and techniques for solving this type of equation can be found [here](https://ncn-foreigners.github.io/nonprobsvy-book/).
#'   To use the variable selection model, set the `vars_selection` parameter in the [control_inf()] function to `TRUE`. In addition, in the other control functions such as [control_sel()] and [control_out()]
#'   you can set parameters for the selection of the relevant variables, such as the number of folds during cross-validation algorithm or the lambda value for penalizations. Details can be found
#'   in the documentation of the control functions for `nonprob`.
#'
#'
#'
#' @references
#' Kim JK, Park S, Chen Y, Wu C. Combining non-probability and
#' probability survey samples through mass imputation. J R Stat Soc Series A. 2021;184:941–
#' 963.
#'
#' Shu Yang, Jae Kwang Kim, Rui Song. Doubly robust inference when combining probability
#' and non-probability samples with high dimensional data. J. R. Statist. Soc. B (2020)
#'
#' Yilin Chen , Pengfei Li & Changbao Wu (2020) Doubly Robust Inference
#' With Nonprobability Survey Samples, Journal of the American Statistical Association, 115:532,
#' 2011-2021
#'
#' Shu Yang, Jae Kwang Kim and Youngdeok Hwang Integration of data from
#' probability surveys and big found data for finite population inference using mass imputation.
#' Survey Methodology, June 2021 29 Vol. 47, No. 1, pp. 29-58
#'
#' @return Returns an object of the `nonprob` class (it is actually a `list`) which contains the following elements:\cr

#' \itemize{
#'  \item{\code{call} -- the call of the `nonprob` function}
#'  \item{\code{data} -- a `data.frame` passed from the `nonprob` function `data` argument}
#'  \item{\code{X} -- a `model.matrix` containing data from probability (first \eqn{n_{S_B}} rows) and non-probability samples (next \eqn{n_{S_B}} rows) if specified at a function call}
#'  \item{\code{y} -- a `list` of vector of outcome variables if specified at a function call}
#'  \item{\code{R} -- a `numeric vector` indicating whether a unit belongs to the probability (0) or non-probability (1) units in the matrix X}
#'  \item{\code{ps_scores} -- a `numeric vector` of estimated propensity scores for probability and non-probability sample}
#'  \item{\code{case_weights} -- a `vector` of case weights for non-probability sample based on the call}
#'  \item{\code{ipw_weights} -- a `vector` of inverse probability weights for non-probability sample (if applicable)}
#'  \item{\code{control} -- a `list` of control functions based on the call}
#'  \item{\code{output} -- a `data.frame` with the estimated means and standard errors for the variables specified in the `target` or `outcome` arguments}
#'  \item{\code{SE} -- a `data.frame` with standard error of the estimator of the population mean, divided into errors from probability and non-probability samples (if applicable)}
#'  \item{\code{confidence_interval} -- a `data.frame` with confidence interval of population mean estimator}
#'  \item{\code{nonprob_size} -- a scalar `numeric vector` denoting the size of non-probability sample}
#'  \item{\code{prob_size} -- a scalar `numeric vector` denoting the size of probability sample}
#'  \item{\code{pop_size} -- a scalar `numeric vector` estimated population size derived from estimated weights (non-probability sample) or known design weights (probability sample)}
#'  \item{\code{pop_size_fixed} -- a `logical` value whether the population size was fixed (known) or estimated (unknown)}
#'  \item{\code{pop_totals} -- a `numeric vector` with the total values of the auxiliary variables derived from a probability sample or based on the call}
#'  \item{\code{pop_means} -- a `numeric vector` with the mean values of the auxiliary variables derived from a probability sample or based on the call}
#'  \item{\code{outcome} -- a `list` containing information about the fitting of the mass imputation model. Structure of the object is based on the `method_outcome` and `family_outcome` arguments which point to specific methods as defined by functions `method_*` (if specified in the call)}
#'  \item{\code{selection} -- a `list` containing information about the fitting of the propensity score model. Structure of the object is based on the `method_selection` argument which point to specific methods as defined by functions `method_ps` (if specified in the call)}
#'  \item{\code{boot_sample} -- a `matrix` with bootstrap estimates of the target variable(s) (if specified)}
#'  \item{\code{svydesign} -- a `svydesign2` object (if specified)}
#'  \item{\code{ys_rand_pred} -- a `list` of predicted values for the target variable(s) for the probability sample (for the MI and DR estimator)}
#'  \item{\code{ys_nons_pred} -- a `list` of predicted values for the target variable(s) for the non-probability sample (for the MI and DR estimator)}
#'  \item{\code{ys_resid} -- a `list` of residuals for the target variable(s) for the non-probability sample (for the MI and DR estimator)}
#'  \item{\code{estimator} -- a `character vector` with information what type of estimator was selected (one of `c("ipw", "mi", "dr")`).}
#'  \item{\code{selection_formula} -- a `formula` based on the `selection` argument (if specified)}
#'  \item{\code{estimator_method} -- a `character vector` with information on the detailed method applied (for the `print` method)}
#' }
#' @seealso
#' [stats::optim()] -- For more information on the \code{optim} function used in the
#' \code{optim} method of propensity score model fitting.
#'
#' [maxLik::maxLik()] -- For more information on the \code{maxLik} function used in
#' \code{maxLik} method of propensity score model fitting.
#'
#' [ncvreg::cv.ncvreg()] -- For more information on the \code{cv.ncvreg} function used in
#' variable selection for the outcome model.
#'
#' [nleqslv::nleqslv()] -- For more information on the \code{nleqslv} function used in
#' estimation process of the bias minimization approach.
#'
#' [stats::glm()] -- For more information about the generalised linear models used during mass imputation process.
#'
#' [RANN::nn2()] -- For more information about the nearest neighbour algorithm used during mass imputation process.
#'
#' [control_sel()] -- For the control parameters related to selection model.
#'
#' [control_out()] -- For the control parameters related to outcome model.
#'
#' [control_inf()] -- For the control parameters related to statistical inference.

#' @examples
#' \donttest{
#' # generate data based on Doubly Robust Inference With Non-probability Survey Samples (2021)
#' # Yilin Chen , Pengfei Li & Changbao Wu
#' library(sampling)
#' set.seed(123)
#' # sizes of population and probability sample
#' N <- 20000 # population
#' n_b <- 1000 # probability
#' # data
#' z1 <- rbinom(N, 1, 0.7)
#' z2 <- runif(N, 0, 2)
#' z3 <- rexp(N, 1)
#' z4 <- rchisq(N, 4)
#'
#' # covariates
#' x1 <- z1
#' x2 <- z2 + 0.3 * z2
#' x3 <- z3 + 0.2 * (z1 + z2)
#' x4 <- z4 + 0.1 * (z1 + z2 + z3)
#' epsilon <- rnorm(N)
#' sigma_30 <- 10.4
#' sigma_50 <- 5.2
#' sigma_80 <- 2.4
#'
#' # response variables
#' y30 <- 2 + x1 + x2 + x3 + x4 + sigma_30 * epsilon
#' y50 <- 2 + x1 + x2 + x3 + x4 + sigma_50 * epsilon
#' y80 <- 2 + x1 + x2 + x3 + x4 + sigma_80 * epsilon
#'
#' # population
#' sim_data <- data.frame(y30, y50, y80, x1, x2, x3, x4)
#' ## propensity score model for non-probability sample (sum to 1000)
#' eta <- -4.461 + 0.1 * x1 + 0.2 * x2 + 0.1 * x3 + 0.2 * x4
#' rho <- plogis(eta)
#'
#' # inclusion probabilities for probability sample
#' z_prob <- x3 + 0.2051
#' sim_data$p_prob <- inclusionprobabilities(z_prob, n = n_b)
#'
#' # data
#' sim_data$flag_nonprob <- UPpoisson(rho) ## sampling nonprob
#' sim_data$flag_prob <- UPpoisson(sim_data$p_prob) ## sampling prob
#' nonprob_df <- subset(sim_data, flag_nonprob == 1) ## non-probability sample
#' svyprob <- svydesign(
#'   ids = ~1, probs = ~p_prob,
#'   data = subset(sim_data, flag_prob == 1),
#'   pps = "brewer"
#' ) ## probability sample
#'
#' ## mass imputation estimator
#' mi_res <- nonprob(
#'   outcome = y30 + y50 + y80 ~ x1 + x2 + x3 + x4,
#'   data = nonprob_df,
#'   svydesign = svyprob
#' )
#' mi_res
#' ## inverse probability weighted estimator
#' ipw_res <- nonprob(
#'   selection = ~ x1 + x2 + x3 + x4,
#'   target = ~y30 + y50 + y80,
#'   data = nonprob_df,
#'   svydesign = svyprob
#' )
#' ipw_res
#' ## doubly robust estimator
#' dr_res <- nonprob(
#'   outcome = y30 + y50 + y80 ~ x1 + x2 + x3 + x4,
#'   selection = ~ x1 + x2 + x3 + x4,
#'   data = nonprob_df,
#'   svydesign = svyprob
#' )
#' dr_res
#' }
#' @name nonprob
NULL
