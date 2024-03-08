#' @title Control parameters for selection model
#' @author Łukasz Chrostowski, Maciej Beręsewicz
#' \loadmathjax
#'
#' @description \code{controlSel} constructs a list with all necessary control parameters
#' for selection model.
#'
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
#' \item if \code{zero} then start is a vector of zeros.
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
                       start_type = c("glm", "naive", "zero")) {


  list(
    epsilon = epsilon,
    maxit = maxit,
    trace = trace,
    optimizer = if (missing(optimizer)) "optim" else optimizer,
    maxLik_method = maxLik_method,
    optim_method = optim_method,
    dependence = dependence,
    key = key,
    est_method_sel = if (missing(est_method_sel)) "mle" else est_method_sel,
    h = if (missing(h)) 1 else h,
    penalty = if (missing(penalty)) "SCAD" else penalty,
    a_SCAD = a_SCAD,
    a_MCP = a_MCP,
    lambda_min = lambda_min,
    nlambda = nlambda,
    nfolds = nfolds,
    lambda = lambda,
    print_level = print_level,
    start_type = if (missing(start_type)) "naive" else start_type
  )
}




