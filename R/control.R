#' @title Control parameters for selection model
#' @author Łukasz Chrostowski, Maciej Beręsewicz
#' @description \code{controlSel} constructs a list with all necessary control parameters
#' for selection model.
#'
#' \loadmathjax
#'
#' @param method estimation method.
#' @param epsilon Tolerance for fitting algorithms by default \code{1e-6}.
#' @param maxit Maximum number of iterations.
#' @param trace logical value. If `TRUE` trace steps of the fitting algorithms. Default is `FALSE`
#' @param optimizer - optimization function for maximum likelihood estimation.
#' @param optim_method maximisation method that will be passed to [stats::optim()] function. Default is `BFGS`.
#' @param maxLik_method maximisation method that will be passed to [maxLik::maxLik()] function. Default is `NR`.
#' @param dependence logical value - `TRUE` if samples are dependent.
#' @param key binary key variable
#' @param est_method_sel Method of estimation for propensity score model.
#' @param h Smooth function for the generalized estimating equations methods taking the following values
#' \itemize{
#'   \item if \code{1} then \mjseqn{\mathbf{h}\left(\mathbf{x}, \boldsymbol{\theta}\right) =
#'   \frac{\pi(\mathbf{x}, \boldsymbol{\theta})}{\mathbf{x}}}
#'   \item if \code{2} then \mjseqn{ \mathbf{h}\left(\mathbf{x}, \boldsymbol{\theta}\right) = \mathbf{x}}
#'   }
#' @param penalty The penanlization function used during variables selection.
#' @param a_SCAD The tuning parameter of the SCAD penalty for selection model. Default is 3.7.
#' @param a_MCP The tuning parameter of the MCP penalty for selection model. Default is 3.
#' @param lambda A user-specified \mjseqn{\lambda} value during variable selection model fitting.
#' @param lambda_min The smallest value for lambda, as a fraction of `lambda.max`. Default is .001.
#' @param nlambda The number of `lambda` values. Default is 50.
#' @param nfolds The number of folds for cross validation. Default is 10.
#' @param print_level this argument determines the level of printing which is done during the optimization (for propensity score model) process.
#' @param start_type - Type of method for start points for model fitting taking the following values
#' \itemize{
#' \item if \code{glm} then start taken from the glm function called on samples.
#' \item if \code{naive} then start consists of a vector which has the value of an estimated parameter for one-dimensional data (on intercept) and 0 for the rest.
#' }
#'
#' @return List with selected parameters.
#'
#' @seealso
#'
#' [nonprob()] -- for fitting procedure with non-probability samples.
#'
#' @export

controlSel <- function(method = "glm.fit", # perhaps another control function for model with variables selection
                       epsilon = 1e-4,
                       maxit = 500,
                       trace = FALSE,
                       optimizer = c("maxLik", "optim"),
                       maxLik_method = "NR",
                       optim_method = "BFGS",
                       dependence = FALSE,
                       key = NULL,
                       est_method_sel = c("mle", "gee"),
                       h = c(1, 2),
                       penalty = c("SCAD", "lasso", "MCP"),
                       a_SCAD = 3.7,
                       a_MCP = 3,
                       lambda = -1,
                       lambda_min = .001,
                       nlambda = 50,
                       nfolds = 10,
                       print_level = 0,
                       start_type = c("glm", "naive")
                       ) {

  list(epsilon = epsilon,
       maxit = maxit,
       trace = trace,
       optimizer = if(missing(optimizer)) "optim" else optimizer,
       maxLik_method = maxLik_method,
       optim_method = optim_method,
       dependence = dependence,
       key = key,
       est_method_sel = if(missing(est_method_sel)) "mle" else est_method_sel,
       h = if(missing(h)) 1 else h,
       penalty = if(missing(penalty)) "SCAD" else penalty,
       a_SCAD = a_SCAD,
       a_MCP = a_MCP,
       lambda_min = lambda_min,
       nlambda = nlambda,
       nfolds = nfolds,
       lambda = lambda,
       print_level = print_level,
       start_type = if(missing(start_type)) "naive" else start_type
      )

}

#' @title Control parameters for outcome model
#' @description \code{controlOut} constructs a list with all necessary control parameters
#' for outcome model.
#' @param epsilon Tolerance for fitting algorithms. Default is \code{1e-6}.
#' @param maxit Maximum number of iterations.
#' @param trace logical value. If `TRUE` trace steps of the fitting algorithms. Default is `FALSE`.
#' @param k The k parameter in the [RANN::nn2()] function. Default is 5.
#' @param penalty penalty algorithm for variable selection. Default is `SCAD`
#' @param a_SCAD The tuning parameter of the SCAD penalty for outcome model. Default is 3.7.
#' @param a_MCP The tuning parameter of the MCP penalty for outcome model. Default is 3.
#' @param lambda_min The smallest value for lambda, as a fraction of lambda.max. Default is .001.
#' @param nlambda The number of lambda values. Default is 100.
#' @param nfolds The number of folds during cross-validation for variables selection model.
#' @param treetype Type of tree for nearest neighbour imputation passed to [RANN::nn2()] function.
#' @param searchtype Type of search for nearest neighbour imputation passed to [RANN::nn2()] function.
#' @param predictive_match Indicates how to select 'closest' unit from
#' nonprobability sample for each unit in probability sample. Either \code{1}
#' (default) or \code{2} where \code{1} is matching by minimising distance
#' between \mjseqn{\hat{y}_{i}} for \mjseqn{i \in S_{A}} and
#' \mjseqn{y_{j}} for \mjseqn{j \in S_{B}} and \code{2} is matching by
#' minimising distance between \mjseqn{\hat{y}_{i}} for \mjseqn{i \in S_{A}}
#' and \mjseqn{\hat{y}_{i}} for \mjseqn{i \in S_{A}}.
#'
#' @return List with selected parameters.
#'
#' @seealso
#'
#' [nonprob()] -- for fitting procedure with non-probability samples.
#'
#'
#' @export

controlOut <- function(epsilon = 1e-6,
                       maxit = 100,
                       trace = FALSE,
                       k = 5,
                       penalty = c("SCAD", "lasso", "MCP"),
                       a_SCAD = 3.7,
                       a_MCP = 3,
                       lambda_min = .001,
                       nlambda = 100,
                       nfolds = 10,
                       treetype = "kd",
                       searchtype = "standard",
                       predictive_match = 1:2
                       ) {

  if (missing(predictive_match))
    predictive_match <- 1

  list(epsilon = epsilon,
       maxit = maxit,
       trace = trace,
       k = k,
       penalty = if(missing(penalty)) "SCAD" else penalty,
       a_SCAD = a_SCAD,
       a_MCP = a_MCP,
       lambda_min = lambda_min,
       nlambda = nlambda,
       nfolds = nfolds,
       treetype = treetype,
       searchtype = searchtype,
       predictive_match = predictive_match)

}


#' @title Control parameters for inference
#' @description \code{controlInf} constructs a list with all necessary control parameters
#' for statistical inference.
#' @param vars_selection If `TRUE`, then variables selection model is used.
#' @param var_method variance method.
#' @param rep_type replication type for weights in the bootstrap method for variance estimation passed to [survey::as.svrepdesign()].
#'  Default is `subbootstrap`.
#' @param bias_inf inference method in the bias minimization.
#' \itemize{
#'   \item if \code{union} then final model is fitting on union of selected variables for selection and outcome models
#'   \item if \code{div} then final model is fitting separately on division of selected variables into relevant ones for
#'   selection and outcome model.
#'   }
#' @param bias_correction if `TRUE`, then bias minimization estimation used during fitting the model.
#' @param num_boot number of iteration for bootstrap algorithms.
#' @param alpha Significance level, Default is 0.05.
#' @param cores Number of cores in parallel computing.
#'
#'
#' @return List with selected parameters.
#'
#' @seealso
#'
#' [nonprob()] -- for fitting procedure with non-probability samples.
#'
#' @export

controlInf <- function(vars_selection = FALSE,
                       var_method = c("analytic",
                                      "bootstrap"),
                       rep_type = c("auto", "JK1", "JKn", "BRR", "bootstrap",
                                    "subbootstrap","mrbbootstrap","Fay"),
                       bias_inf = c("union", "div"),
                       num_boot = 500,
                       bias_correction = FALSE,
                       alpha = 0.05,
                       cores = 1) {

  list(vars_selection = if(missing(vars_selection)) FALSE else vars_selection,
       var_method = if(missing(var_method)) "analytic" else var_method,
       rep_type = if(missing(rep_type)) "subbootstrap" else rep_type,
       bias_inf = if(missing(bias_inf)) "union" else bias_inf,
       bias_correction = bias_correction,
       num_boot = num_boot,
       alpha = alpha,
       cores = cores)

}
