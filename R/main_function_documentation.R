#' @import mathjaxr
NULL
#' @title Inference with the non-probability survey samples
#' @author Łukasz Chrostowski, Maciej Beręsewicz
#'
#' \loadmathjax
#' @description \code{nonprob} fits model for inference based on non-probability surveys (including big data) using various methods.
#' The function allows you to estimate the population mean having access to a reference probability sample as well as total/mean values of covariates.
#' In the package implemented state-of-the-art approaches recently proposed in the literature: Chen et al. (2020), Yang et al. (2020), Wu (2022) and use `survey` package [Lumley 2004](https://cran.r-project.org/web/packages/survey/index.html) for inference.
#' Provided propensity score weighting (e.g. with calibration constraints), mass imputation (e.g. nearest neighbour) and doubly robust estimators that take into account minimization of the asymptotic bias of the population mean estimators,
#' variable selection or overlap between random and non-random sample.
#' The package uses `survey` package functionalities when a probability sample is available.
#'
#'
#' @param data `data.frame` with data from the nonprobability sample.
#' @param selection `formula`, the selection (propensity) equation.
#' @param outcome `formula`, the outcome equation.
#' @param target `formula` with target variables.
#' @param svydesign an optional `svydesign` object (from the survey package) containing probability sample and design weights.
#' @param pop_totals an optional `named vector` with population totals of the covariates.
#' @param pop_means an optional `named vector` with population means of the covariates.
#' @param pop_size an optional `double` with population size.
#' @param overlap logical value - `TRUE` if samples overlap.
#' @param method_selection a `character` with method for propensity scores estimation
#' @param method_outcome a `character` with method for response variable estimation
#' @param family_selection a `character` string describing the error distribution and link function to be used in the model. Default is "binomial". Currently only binomial with logit link is supported.
#' @param family_outcome a `character` string describing the error distribution and link function to be used in the model. Default is "gaussian". Currently supports: gaussian with identity link, poisson and binomial.
#' @param subset an optional `vector` specifying a subset of observations to be used in the fitting process.
#' @param strata an optional `vector` specifying strata.
#' @param weights an optional `vector` of prior weights to be used in the fitting process. Should be NULL or a numeric vector. It is assumed that this vector contains frequency or analytic weights
#' @param na_action a function which indicates what should happen when the data contain `NAs`.
#' @param control_selection a list indicating parameters to use in fitting selection model for propensity scores
#' @param control_outcome a list indicating parameters to use in fitting model for outcome variable
#' @param control_inference a list indicating parameters to use in inference based on probability and non-probability samples, contains parameters such as estimation method or variance method
#' @param start an optional `list` with starting values for the parameters of the selection and outcome equation
#' @param verbose verbose, numeric
#' @param x Logical value indicating whether to return model matrix of covariates as a part of output.
#' @param y Logical value indicating whether to return vector of outcome variable as a part of output.
#' @param ... Additional, optional arguments.
#'
#' @details Let \mjseqn{y} be the response variable for which we want to estimate population mean
#' given by \mjsdeqn{\mu_{y} = \frac{1}{N}\sum_{i=1}^N y_i.} For this purposes we consider data integration
#' with the following structure. Let \mjseqn{S_A} be the non-probability sample with design matrix of covariates as
#' \mjsdeqn{
#' \begin{equation}
#' \boldsymbol{X}_A =
#'   \begin{bmatrix}
#' x_{11} & x_{12} & \cdots & x_{1p} \cr
#' x_{21} & x_{22} & \cdots & x_{2p} \cr
#' \vdots & \vdots & \ddots & \vdots \cr
#' x_{n_{A}1} & x_{n_{A}2} & \cdots & x_{n_{A}p} \cr
#' \end{bmatrix}
#' \end{equation}
#' }
#' and vector of outcome variable
#' \mjsdeqn{
#' \begin{equation}
#' \boldsymbol{y} =
#'   \begin{bmatrix}
#' y_{1} \cr
#' y_{2} \cr
#' \vdots \cr
#' y_{n_{A}}.
#' \end{bmatrix}
#' \end{equation}
#' }
#' On the other hand let \mjseqn{S_B} be the probability sample with design matrix of covariates as
#' \mjsdeqn{
#' \begin{equation}
#' \boldsymbol{X}_B =
#'   \begin{bmatrix}
#' x_{11} & x_{12} & \cdots & x_{1p} \cr
#' x_{21} & x_{22} & \cdots & x_{2p} \cr
#' \vdots & \vdots & \ddots & \vdots \cr
#' x_{n_{B}1} & x_{n_{B}2} & \cdots & x_{n_{B}p}. \cr
#' \end{bmatrix}
#' \end{equation}
#' }
#' Instead of sample of units we can consider vector of population totals in the form of \mjseqn{\tau_x = (\sum_{i \in \mathcal{U}}\boldsymbol{x}_{i1}, \sum_{i \in \mathcal{U}}\boldsymbol{x}_{i2}, ..., \sum_{i \in \mathcal{U}}\boldsymbol{x}_{ip})} or means
#' \mjseqn{\frac{\tau_x}{N}}, where \mjseqn{\mathcal{U}} regards to some finite population. Notice that we do not assume access to the response variable for \mjseqn{S_B}.
#' Generally we provide following assumptions:
#' 1.  The selection indicator of belonging to non-probability sample \mjseqn{R_i} and the response variable \mjseqn{y_i} are independent given the set of covariates \mjseqn{\boldsymbol{x}_i}.
#' 2.  All units have a non-zero propensity score, that is, \mjseqn{\pi_i^A > 0} for all i.
#' 3.  The indicator variables \mjseqn{R_i^A} and \mjseqn{R_j^A} are independent for given \mjseqn{\boldsymbol{x}_i} and \mjseqn{\boldsymbol{x}_j} for \mjseqn{i \neq j}.
#'
#' There are three possible approaches to problem of population mean estimation using non-probability samples:
#'
#' 1. Inverse probability weighting -- The biggest drawback of the non-probability sampling is unknown selection mechanism for a unit to be included in the sample.
#'  This is why we talk about so called “biased sample” problem.
#'  Inverse probability approach is based on assumption that reference probability sample
#'  is available and therefore we can estimate propensity score of selection mechanism.
#'  Estimator has following form:
#'  \mjsdeqn{\mu_{IPW} = \frac{1}{N^A}\sum_{i \in S_A}\frac{y_i}{\hat{\pi}_{i}^A}.}
#'  For this purpose with consider multiple ways of estimation. The first approach is Maximum Likelihood Estimation with corrected
#'  log-likelihood function which is given by the following formula
#'  \mjsdeqn{
#'  \begin{equation}
#'  \ell^{*}(\boldsymbol{\theta}) = \sum_{i \in S_{A}}\log \left\lbrace \frac{\pi(\boldsymbol{x}_{i}, \boldsymbol{\theta})}{1 - \pi(\boldsymbol{x}_{i},\boldsymbol{\theta})}\right\rbrace + \sum_{i \in S_{B}}d_{i}^{B}\log \left\lbrace 1 - \pi({\boldsymbol{x}_{i},\boldsymbol{\theta})}\right\rbrace.
#'  \end{equation}}
#'  In the literature main approach is based on `logit` link function when it comes to modelling propensity scores \mjseqn{\pi_i^A},
#'  however we expand propensity score model with the additional link functions, such as `cloglog` and `probit`.
#'  The pseudo score equations derived from ML methods can be replaced by idea of generalized estimating equations with calibration constraints defined by equations
#'  \mjsdeqn{
#'  \begin{equation}
#'  \mathbf{U}(\boldsymbol{\theta})=\sum_{i \in S_A} \mathbf{h}\left(\mathbf{x}_i, \boldsymbol{\theta}\right)-\sum_{i \in S_B} d_i^B \pi\left(\mathbf{x}_i, \boldsymbol{\theta}\right) \mathbf{h}\left(\mathbf{x}_i, \boldsymbol{\theta}\right).
#'  \end{equation}}
#'  Notice that for \mjseqn{ \mathbf{h}\left(\mathbf{x}_i, \boldsymbol{\theta}\right) = \frac{\pi(\boldsymbol{x}, \boldsymbol{\theta})}{\boldsymbol{x}}} we do not require probability
#'  sample and can use vector of population totals/means.
#'
#' 2. Mass imputation -- This method relies on framework,
#'    where imputed values of the outcome variables are created for whole probability sample.
#'    In this case we treat big-data sample as a training dataset, which is used to build an imputation model. Using imputed values
#'    for probability sample and design (known) weights, we can build population mean estimator of form:
#'    \mjsdeqn{\mu_{MI} = \frac{1}{N^B}\sum_{i \in S_B} d_i^B \hat{y}_i.}
#'    It opens the the door to very flexible method for imputation model. In the package used generalized linear models from [stats::glm()]
#'    and nearest neighbour algorithm using [RANN::nn2()].
#'
#' 3. Doubly robust estimation -- The IPW and MI estimators are sensible on misspecified models for propensity score and outcome variable respectively.
#'    For this purpose so called doubly-robust methods, which take into account these problems, are presented.
#'    It is quite simple idea of combination propensity score and imputation models during inference which lead to the following estimator
#'    \mjsdeqn{\mu_{DR} = \frac{1}{N^A} \sum_{i \in S_A} \hat{d}_i^A (y_i - \hat{y}_i) + \frac{1}{N^B}\sum_{i \in S_B} d_i^B \hat{y}_i.}
#'    In addition, an approach based directly on bias minimisation has been implemented. Following formula
#'    \mjsdeqn{
#'    \begin{aligned}
#'    bias(\hat{\mu}_{DR}) = & \mathbb{E} (\hat{\mu}_{DR} - \mu) \cr = & \mathbb{E} \left\lbrace \frac{1}{N} \sum_{i=1}^N (\frac{R_i^A}{\pi_i^A (\boldsymbol{x}_i^{\mathrm{T}} \boldsymbol{\theta})}
#'    - 1 ) (y_i - \operatorname{m}(\boldsymbol{x}_i^{\mathrm{T}} \boldsymbol{\beta})) \right\rbrace \cr + & \mathbb{E} \left\lbrace \frac{1}{N} \sum_{i=1}^N (R_i^B d_i^B - 1) \operatorname{m}( \boldsymbol{x}_i^{\mathrm{T}} \boldsymbol{\beta}) \right\rbrace,
#'    \end{aligned}
#'    }
#'    lead us to system of equations
#'    \mjsdeqn{
#'    \begin{equation}
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
#'   \end{equation}
#'   }
#'   where \mjseqn{m\left(\boldsymbol{x}_i, \boldsymbol{\beta}\right)} is mass imputation (regression) model for outcome variable and
#'   propensity scores \mjseqn{\pi_i^A} are estimated using `logit` function for the model. As in the `MLE` and `GEE` approaches we have expanded
#'   this method on `cloglog` and `probit` links.
#'
#'   Since it is not straightforward thing to calculate variances of these estimators, in the literature proposed are asymptotic equivalents
#'   of variances that are derives using Taylor approximation. Details cen be found [here](https://ncn-foreigners.github.io/nonprobsvy-book/intro.html).
#'   Moreover one can use bootstrap approach for variance estimation.
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
#'  \item{\code{y}} -- Vector of outcome variable if specified at a function call.
#'  \item{\code{prob} -- vector of estimated propensity scores for non-probability sample.}
#'  \item{\code{weights} -- vector of estimated weights for non-probability sample.}
#'  \item{\code{control} -- list of control functions.}
#'  \item{\code{output} -- output of the model containing information about estimated population mean and standard errors.}
#'  \item{\code{SE} -- standard error of population mean estimator divided into errors coming from probability and non-probability samples.}
#'  \item{\code{confidence_interval} -- confidence interval of population mean estimator}
#'  \item{\code{coeff_selection} -- `data.frame` of estimated coefficients of the propensity score (inverse probability weighting) model and their standard errors.}
#'  \item{\code{coeff_outcome} -- `data.frame` of estimated coefficients of the outcome (mass imputation) model and their standard errors.}
#'  \item{\code{nonprob_size} -- size of non-probability sample}
#'  \item{\code{prob_size} -- size of probability sample}
#'  \item{\code{pop_size} -- estimated population size derived from estimated weights (non-probability sample) or known design weights (probability sample)}
#'  \item{\code{outcome} -- list containing information about fitting of mass imputation model, in case of regression model, object containing list returned by the function
#'  [stats::glm()], in case of nearest neigbour imputation object containing list returned by [RANN::nn2()].}
#'  \item{\code{selection} -- list containing information about fitting of propensity score model, such as
#'  \itemize{
#'  \item{\code{coefficients} -- a named vector of coefficients}
#'  \item{\code{std_err} -- standard errors of the estimated model coefficients}
#'  \item{\code{residuals} -- the working residuals}
#'  \item{\code{fitted_values} -- the fitted mean values, obtained by transforming the linear predictors by the inverse of the link function.}
#'  \item{\code{link} -- the `link` object used.}
#'  \item{\code{linear_predictors} -- the linear fit on link scale.}
#'  \item{\code{aic} --	A version of Akaike's An Information Criterion, minus twice the maximized log-likelihood plus twice the number of parameters.}
#'  \item{\code{weights} -- vector of estimated weights for non-probability sample.}
#'  \item{\code{prior.weights} -- the weights initially supplied, a vector of 1s if none were.}
#'  \item{\code{formula} -- the formula supplied.}
#'  \item{\code{df_residual} -- the residual degrees of freedom.}
#'  \item{\code{log_likelihood} -- value of log-likelihood function if `mle` method, in the other case `NULL`.}
#'   }
#'  }
#' }
#' @seealso
#' [stats::optim()] -- For more information on \code{optim} function used in
#' \code{optim} method of fitting propensity score model.
#'
#' [maxLik::maxLik()] -- For more information on \code{maxLik} function used in
#' \code{maxLik} method of fitting propensity score model.
#'
#' [ncvreg::cv.ncvreg()] -- For more information on \code{cv.ncvreg} function used in
#' variable selection model for outcome model.
#'
#' [nleqslv::nleqslv()] -- For more information on \code{nleqslv} function used in
#' estimation process of bias minimization approach.
#'
#' [stats::glm()] -- For more information about generalised linear models used during mass imputation process.
#'
#' [RANN::nn2()] -- For more information about nearest neighbor algorithm used during mass imputation process.
#'
#' [controlSel()] -- For control parameters related to selection model.
#'
#' [controlOut()] -- For control parameters related to outcome model.
#'
#' [controlInf()] -- For control parameters related to statistical inference.

#' @examples
#' \donttest{
#' # generate data based on Doubly Robust Inference With Nonprobability Survey Samples (2021)
#' # Yilin Chen , Pengfei Li & Changbao Wu
#' library(sampling)
#' set.seed(123)
#' # sizes of population and probability sample
#' N <-  20000 # population
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
#' y30 <- 2 + x1  + x2 + x3 + x4 + sigma_30*epsilon
#' y50 <- 2 + x1  + x2 + x3 + x4 + sigma_50*epsilon
#' y80 <- 2 + x1  + x2 + x3 + x4 + sigma_80*epsilon
#'
#' # population
#' sim_data <- data.frame(y30, y50, y80, x1, x2, x3, x4)
#' ## propensity score model for nonprobability sample (sum to 1000)
#' eta <- -4.461 + 0.1*x1 + 0.2*x2 + 0.1*x3 + 0.2*x4
#' rho <- plogis(eta)
#'
#' # inclusion probabilities for probability sample
#' z_prob <- x3+0.2051
#' sim_data$p_prob <- inclusionprobabilities(z_prob, n = n_b)
#'
#' # data
#' sim_data$flag_nonprob <- UPpoisson(rho) ## sampling nonprob
#' sim_data$flag_prob <- UPpoisson(sim_data$p_prob) ## sampling prob
#' nonprob_df <- subset(sim_data, flag_nonprob == 1) ## nonprobability sample
#' svyprob <- svydesign(ids=~1, probs = ~ p_prob,
#'                      data = subset(sim_data, flag_prob == 1),
#'                      pps = "brewer") ## probability sample
#'
#' ## mass imputation estimator
#' MI_res <- nonprob(outcome = y80 ~ x1 + x2 + x3 + x4,
#'                   data = nonprob_df,
#'                   svydesign = svyprob)
#' summary(MI_res)
#' ## inverse probability weighted estimator
#' IPW_res <- nonprob(selection = ~ x1 + x2 + x3 + x4,
#'                    target = ~ y80,
#'                    data = nonprob_df,
#'                    svydesign = svyprob)
#' summary(IPW_res)
#' ## doubly robust estimator
#' DR_res <- nonprob(outcome = y80 ~ x1 + x2 + x3 + x4,
#'                   selection = ~ x1 + x2 + x3 + x4,
#'                   data = nonprob_df,
#'                   svydesign = svyprob)
#' summary(DR_res)
#' }
#' @name nonprob
NULL
