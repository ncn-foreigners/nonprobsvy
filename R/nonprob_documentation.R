#' @import mathjaxr
NULL
#' @title Inference with non-probability survey samples
#' @author Łukasz Chrostowski, Maciej Beręsewicz
#'
#' \loadmathjax
#' @description \code{nonprob} fits a model for inference based on non-probability surveys (including big data) using various methods.
#' The function allows you to estimate the population mean with access to a reference probability sample, as well as sums and means of covariates.
#'
#' The package implements state-of-the-art approaches recently proposed in the literature: Chen et al. (2020),
#' Yang et al. (2020), Wu (2022) and uses the [Lumley 2004](https://CRAN.R-project.org/package=survey) `survey` package for inference.
#'
#' It provides propensity score weighting (e.g. with calibration constraints), mass imputation (e.g. nearest neighbour) and
#' doubly robust estimators that take into account minimisation of the asymptotic bias of the population mean estimators or
#' variable selection.
#' The package uses the `survey` package functionality when a probability sample is available.
#'
#'
#' @param data a `data.frame` with data from the non-probability sample.
#' @param selection a `formula`, the selection (propensity) equation.
#' @param outcome a `formula`, the outcome equation.
#' @param target a `formula` with target variables.
#' @param svydesign an optional `svydesign` object (from the survey package) containing a probability sample and design weights.
#' @param pop_totals an optional `named vector` with population totals of the covariates.
#' @param pop_means an optional `named vector` with population means of the covariates.
#' @param pop_size an optional `double` value with population size.
#' @param method_selection a `character` indicating the method for propensity scores estimation.
#' @param method_outcome a `character` indicating the method for response variable estimation.
#' @param family_outcome a `character` string describing the error distribution and the link function to be used in the model, set to `gaussian` by default. Currently supports: gaussian with identity link, poisson and binomial.
#' @param subset an optional `vector` specifying a subset of observations to be used in the fitting process - not yet supported.
#' @param strata an optional `vector` specifying strata - not yet supported.
#' @param weights an optional `vector` of prior weights to be used in the fitting process. Should be NULL or a numeric vector. It is assumed that this vector contains frequency or analytic weights.
#' @param na_action a function which indicates what should happen when the data contain `NAs` - not yet supported.
#' @param control_selection a `list` indicating parameters to be used when fitting the selection model for propensity scores.
#' @param control_outcome a `list` indicating parameters to be used when fitting the model for the outcome variable.
#' @param control_inference a `list` indicating parameters to be used for inference based on probability and non-probability samples, contains parameters such as the estimation method or the variance method.
#' @param start_selection an optional `vector` with starting values for the parameters of the selection equation.
#' @param start_outcome an optional `vector` with starting values for the parameters of the outcome equation.
#' @param verbose verbose, numeric.
#' @param x a logical value indicating whether to return model matrix of covariates as a part of the output.
#' @param y a logical value indicating whether to return vector of the outcome variable as a part of the output.
#' @param se Logical value indicating whether to calculate and return standard error of estimated mean.
#' @param ... Additional, optional arguments.
#'
#' @details Let \mjseqn{y} be the response variable for which we want to estimate the population mean,
#' given by \mjsdeqn{\mu_{y} = \frac{1}{N} \sum_{i=1}^N y_{i}.} For this purpose we consider data integration
#' with the following structure. Let \mjseqn{S_A} be the non-probability sample with the design matrix of covariates as
#' \mjsdeqn{
#' \boldsymbol{X}_A =
#'   \begin{bmatrix}
#' x_{11} & x_{12} & \cdots & x_{1p} \cr
#' x_{21} & x_{22} & \cdots & x_{2p} \cr
#' \vdots & \vdots & \ddots & \vdots \cr
#' x_{n_{A}1} & x_{n_{A}2} & \cdots & x_{n_{A}p} \cr
#' \end{bmatrix}
#' }
#' and vector of outcome variable
#' \mjsdeqn{
#' \boldsymbol{y} =
#'   \begin{bmatrix}
#' y_{1} \cr
#' y_{2} \cr
#' \vdots \cr
#' y_{n_{A}}.
#' \end{bmatrix}
#' }
#' On the other hand, let \mjseqn{S_B} be the probability sample with design matrix of covariates be
#' \mjsdeqn{
#' \boldsymbol{X}_B =
#'   \begin{bmatrix}
#' x_{11} & x_{12} & \cdots & x_{1p} \cr
#' x_{21} & x_{22} & \cdots & x_{2p} \cr
#' \vdots & \vdots & \ddots & \vdots \cr
#' x_{n_{B}1} & x_{n_{B}2} & \cdots & x_{n_{B}p}. \cr
#' \end{bmatrix}
#' }
#' Instead of a sample of units we can consider a vector of population sums in the form of \mjseqn{\tau_x = (\sum_{i \in \mathcal{U}}\boldsymbol{x}_{i1}, \sum_{i \in \mathcal{U}}\boldsymbol{x}_{i2}, ..., \sum_{i \in \mathcal{U}}\boldsymbol{x}_{ip})} or means
#' \mjseqn{\frac{\tau_x}{N}}, where \mjseqn{\mathcal{U}} refers to a finite population. Note that we do not assume access to the response variable for \mjseqn{S_B}.
#' In general we make the following assumptions:
#' 1.  The selection indicator of belonging to non-probability sample \mjseqn{R_{i}} and the response variable \mjseqn{y_i} are independent given the set of covariates \mjseqn{\boldsymbol{x}_i}.
#' 2.  All units have a non-zero propensity score, i.e., \mjseqn{\pi_{i}^{A} > 0} for all i.
#' 3.  The indicator variables \mjseqn{R_{i}^{A}} and \mjseqn{R_{j}^{A}} are independent for given \mjseqn{\boldsymbol{x}_i} and \mjseqn{\boldsymbol{x}_j} for \mjseqn{i \neq j}.
#'
#' There are three possible approaches to the problem of estimating population mean using non-probability samples:
#'
#' 1. Inverse probability weighting - The main drawback of non-probability sampling is the unknown selection mechanism for a unit to be included in the sample.
#'  This is why we talk about the so-called "biased sample" problem. The inverse probability approach is based on the assumption that a reference probability sample
#'  is available and therefore we can estimate the propensity score of the selection mechanism.
#'  The estimator has the following form:
#'  \mjsdeqn{\hat{\mu}_{IPW} = \frac{1}{N^{A}}\sum_{i \in S_{A}} \frac{y_{i}}{\hat{\pi}_{i}^{A}}.}
#'  For this purpose several estimation methods can be considered. The first approach is maximum likelihood estimation with a corrected
#'  log-likelihood function, which is given by the following formula
#'  \mjsdeqn{
#'  \ell^{*}(\boldsymbol{\theta}) = \sum_{i \in S_{A}}\log \left\lbrace \frac{\pi(\boldsymbol{x}_{i}, \boldsymbol{\theta})}{1 - \pi(\boldsymbol{x}_{i},\boldsymbol{\theta})}\right\rbrace + \sum_{i \in S_{B}}d_{i}^{B}\log \left\lbrace 1 - \pi({\boldsymbol{x}_{i},\boldsymbol{\theta})}\right\rbrace.}
#'  In the literature, the main approach to modelling propensity scores is based on the logit link function.
#'  However, we extend the propensity score model with the additional link functions such as cloglog and probit.
#'  The pseudo-score equations derived from ML methods can be replaced by the idea of generalised estimating equations
#'  with calibration constraints defined by equations.
#'  \mjsdeqn{
#'  \mathbf{U}(\boldsymbol{\theta})=\sum_{i \in S_A} \mathbf{h}\left(\mathbf{x}_i, \boldsymbol{\theta}\right)-\sum_{i \in S_B} d_i^B \pi\left(\mathbf{x}_i, \boldsymbol{\theta}\right) \mathbf{h}\left(\mathbf{x}_i, \boldsymbol{\theta}\right).}
#'  Notice that for \mjseqn{ \mathbf{h}\left(\mathbf{x}_i, \boldsymbol{\theta}\right) = \frac{\pi(\boldsymbol{x}, \boldsymbol{\theta})}{\boldsymbol{x}}} We do not need a probability sample and can use a vector of population totals/means.
#'
#' 2. Mass imputation -- This method is based on a framework where imputed values of outcome variables are created for the entire probability sample. In this case, we treat the large sample as a training data set that is used to build an imputation model.
#'    Using the imputed values for the probability sample and the (known) design weights,
#'    we can build a population mean estimator of the form:
#'    \mjsdeqn{\hat{\mu}_{MI} = \frac{1}{N^B}\sum_{i \in S_{B}} d_{i}^{B} \hat{y}_i.}
#'    It opens the door to a very flexible method for imputation models. The package uses generalized linear models from [stats::glm()],
#'    the nearest neighbour algorithm using [RANN::nn2()] and predictive mean matching.
#'
#' 3. Doubly robust estimation -- The IPW and MI estimators are sensitive to misspecified models for the propensity score and outcome variable, respectively.
#'    To this end, so-called doubly robust methods are presented that take these problems into account.
#'    It is a simple idea to combine propensity score and imputation models during inference, leading to the following estimator
#'    \mjsdeqn{\hat{\mu}_{DR} = \frac{1}{N^A}\sum_{i \in S_A} \hat{d}_i^A (y_i - \hat{y}_i) + \frac{1}{N^B}\sum_{i \in S_B} d_i^B \hat{y}_i.}
#'    In addition, an approach based directly on bias minimisation has been implemented. The following formula
#'    \mjsdeqn{
#'    \begin{aligned}
#'    bias(\hat{\mu}_{DR}) = & \mathbb{E} (\hat{\mu}_{DR} - \mu) \cr = & \mathbb{E} \left\lbrace \frac{1}{N} \sum_{i=1}^N (\frac{R_i^A}{\pi_i^A (\boldsymbol{x}_i^{\mathrm{T}} \boldsymbol{\theta})}
#'    - 1 ) (y_i - \operatorname{m}(\boldsymbol{x}_i^{\mathrm{T}} \boldsymbol{\beta})) \right\rbrace \cr + & \mathbb{E} \left\lbrace \frac{1}{N} \sum_{i=1}^N (R_i^B d_i^B - 1) \operatorname{m}( \boldsymbol{x}_i^{\mathrm{T}} \boldsymbol{\beta}) \right\rbrace,
#'    \end{aligned}
#'    }
#'    lead us to system of equations
#'    \mjsdeqn{
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
#'   where \mjseqn{m\left(\boldsymbol{x}_{i}, \boldsymbol{\beta}\right)} is a mass imputation (regression) model for the outcome variable and
#'   propensity scores \mjseqn{\pi_i^A} are estimated using a `logit` function for the model. As with the `MLE` and `GEE` approaches we have extended
#'   this method to `cloglog` and `probit` links.
#'
#'   As it is not straightforward to calculate the variances of these estimators, asymptotic equivalents of the variances derived using the Taylor approximation have been proposed in the literature.
#'   Details can be found [here](https://ncn-foreigners.github.io/nonprobsvy-book/intro.html).
#'   In addition, a bootstrap approach can be used for variance estimation.
#'
#'   The function also allows variables selection using known methods that have been implemented to handle the integration of probability and non-probability sampling.
#'   In the presence of high-dimensional data, variable selection is important, because it can reduce the variability in the estimate that results from using irrelevant variables to build the model.
#'   Let \mjseqn{\operatorname{U}\left( \boldsymbol{\theta}, \boldsymbol{\beta} \right)} be the joint estimating function for \mjseqn{\left( \boldsymbol{\theta}, \boldsymbol{\beta} \right)}. We define the
#'   penalized estimating functions as
#'   \mjsdeqn{
#'   \operatorname{U}^p \left(\boldsymbol{\theta}, \boldsymbol{\beta}\right) = \operatorname{U}\left(\boldsymbol{\theta}, \boldsymbol{\beta}\right) - \left\lbrace \begin{array}{c}
#'   q_{\lambda_\theta}(|\boldsymbol{\theta}|) \operatorname{sgn}(\boldsymbol{\theta}) \\
#'   q_{\lambda_\beta}(|\\boldsymbol{\beta}|) \operatorname{sgn}(\boldsymbol{\beta})
#'   \end{array} \right\rbrace,
#'   }
#'   where \mjseqn{\lambda_{\theta}} and \mjseqn{q_{\lambda_{\beta}}} are some smooth functions. We let \mjseqn{q_{\lambda} \left(x\right) = \frac{\partial p_{\lambda}}{\partial x}}, where \mjseqn{p_{\lambda}} is some penalization function.
#'   Details of penalization functions and techniques for solving this type of equation can be found [here](https://ncn-foreigners.github.io/nonprobsvy-book/variableselection.html).
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
#' @return Returns an object of class \code{c("nonprobsvy", "nonprobsvy_dr")} in case of doubly robust estimator,
#' \code{c("nonprobsvy", "nonprobsvy_mi")} in case of mass imputation estimator and
#' \code{c("nonprobsvy", "nonprobsvy_ipw")} in case of inverse probability weighting estimator
#' with type \code{list} containing:\cr
#' \itemize{
#'  \item{\code{X} -- model matrix containing data from probability and non-probability samples if specified at a function call.}
#'  \item{\code{y}} -- list of vector of outcome variables if specified at a function call.
#'  \item{\code{R}} -- vector indicating the probablistic (0) or non-probablistic (1) units in the matrix X.
#'  \item{\code{prob} -- vector of estimated propensity scores for non-probability sample.}
#'  \item{\code{weights} -- vector of estimated weights for non-probability sample.}
#'  \item{\code{control} -- list of control functions.}
#'  \item{\code{output} -- output of the model with information on the estimated population mean and standard errors.}
#'  \item{\code{SE} -- standard error of the estimator of the population mean, divided into errors from probability and non-probability samples.}
#'  \item{\code{confidence_interval} -- confidence interval of population mean estimator.}
#'  \item{\code{nonprob_size} -- size of non-probability sample.}
#'  \item{\code{prob_size} -- size of probability sample.}
#'  \item{\code{pop_size} -- estimated population size derived from estimated weights (non-probability sample) or known design weights (probability sample).}
#'  \item{\code{pop_totals} -- the total values of the auxiliary variables derived from a probability sample or vector of total/mean values.}
#'  \item{\code{outcome} -- list containing information about the fitting of the mass imputation model, in the case of regression model the object containing the list returned by
#'  [stats::glm()], in the case of the nearest neighbour imputation the object containing list returned by [RANN::nn2()]. If `bias_correction` in [control_inf()] is set to `TRUE`, the estimation is based on
#'  the joint estimating equations for the `selection` and `outcome` model and therefore, the list is different from the one returned by the [stats::glm()] function and contains elements such as
#'  \itemize{
#'  \item{\code{coefficients} -- estimated coefficients of the regression model.}
#'  \item{\code{std_err} -- standard errors of the estimated coefficients.}
#'  \item{\code{residuals} -- The response residuals.}
#'  \item{\code{variance_covariance} -- The variance-covariance matrix of the coefficient estimates.}
#'  \item{\code{df_residual} -- The degrees of freedom for residuals.}
#'  \item{\code{family} -- specifies the error distribution and link function to be used in the model.}
#'  \item{\code{fitted.values} -- The predicted values of the response variable based on the fitted model.}
#'  \item{\code{linear.predictors} -- The linear fit on link scale.}
#'  \item{\code{X} -- The design matrix.}
#'  \item{\code{method} -- set on `glm`, since the regression method.}
#'  \item{\code{model_frame} -- Matrix of data from probability sample used for mass imputation.}
#'  }
#'  }
#'  In addition, if the variable selection model for the outcome variable is fitting, the list includes the
#'  \itemize{
#'  \item{\code{cve} -- the error for each value of `lambda`, averaged across the cross-validation folds.}
#'  }
#'  \item{\code{selection} -- list containing information about fitting of propensity score model, such as
#'  \itemize{
#'  \item{\code{coefficients} -- a named vector of coefficients.}
#'  \item{\code{std_err} -- standard errors of the estimated model coefficients.}
#'  \item{\code{residuals} -- the response residuals.}
#'  \item{\code{variance} -- the root mean square error.}
#'  \item{\code{fitted_values} -- the fitted mean values, obtained by transforming the linear predictors by the inverse of the link function.}
#'  \item{\code{link} -- the `link` object used.}
#'  \item{\code{linear_predictors} -- the linear fit on link scale.}
#'  \item{\code{aic} --	A version of Akaike's An Information Criterion, minus twice the maximized log-likelihood plus twice the number of parameters.}
#'  \item{\code{weights} -- vector of estimated weights for non-probability sample.}
#'  \item{\code{prior.weights} -- the weights initially supplied, a vector of 1s if none were.}
#'  \item{\code{est_totals} -- the estimated total values of auxiliary variables derived from a non-probability sample}.
#'  \item{\code{formula} -- the formula supplied.}
#'  \item{\code{df_residual} -- the residual degrees of freedom.}
#'  \item{\code{log_likelihood} -- value of log-likelihood function if `mle` method, in the other case `NA`.}
#'  \item{\code{cve} -- the error for each value of the `lambda`, averaged across the cross-validation folds for the variable selection model
#'  when the propensity score model is fitting. Returned only if selection of variables for the model is used.}
#'  \item{\code{method_selection} -- Link function, e.g. `logit`, `cloglog` or `probit`.}
#'  \item{\code{hessian} -- Hessian Gradient of the log-likelihood function from `mle` method}.
#'  \item{\code{gradient} -- Gradient of the log-likelihood function from `mle` method.}
#'  \item{\code{method} -- An estimation method for selection model, e.g. `mle` or `gee`.}
#'  \item{\code{prob_der} -- Derivative of the inclusion probability function for units in a non--probability sample.}
#'  \item{\code{prob_rand} -- Inclusion probabilities for unit from a probabiliy sample from `svydesign` object.}
#'  \item{\code{prob_rand_est} -- Inclusion probabilites to a non--probabiliy sample for unit from probability sample.}
#'  \item{\code{prob_rand_est_der} -- Derivative of the inclusion probabilites to a non--probabiliy sample for unit from probability sample.}
#'   }
#'  }
#'  \item{\code{stat} -- matrix of the estimated population means in each bootstrap iteration.
#'                       Returned only if a bootstrap method is used to estimate the variance and \code{keep_boot} in
#'                       [control_inf()] is set on `TRUE`.}
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
#' MI_res <- nonprob(
#'   outcome = y80 ~ x1 + x2 + x3 + x4,
#'   data = nonprob_df,
#'   svydesign = svyprob
#' )
#' summary(MI_res)
#' ## inverse probability weighted estimator
#' IPW_res <- nonprob(
#'   selection = ~ x1 + x2 + x3 + x4,
#'   target = ~y80,
#'   data = nonprob_df,
#'   svydesign = svyprob
#' )
#' summary(IPW_res)
#' ## doubly robust estimator
#' DR_res <- nonprob(
#'   outcome = y80 ~ x1 + x2 + x3 + x4,
#'   selection = ~ x1 + x2 + x3 + x4,
#'   data = nonprob_df,
#'   svydesign = svyprob
#' )
#' summary(DR_res)
#' }
#' @name nonprob
NULL
