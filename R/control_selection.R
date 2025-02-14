#' @title Control parameters for the selection model
#'
#' @description \code{control_sel} constructs a list with all necessary control parameters
#' for selection model.
#'
#' @param est_method Method of estimation for propensity score model.
#' @param gee_h_fun Smooth function for the generalized estimating equations (GEE) method taking the following values
#' \itemize{
#'   \item if \code{1} then \mjseqn{\mathbf{h}\left(\mathbf{x}, \boldsymbol{\theta}\right) =
#'   \frac{\pi(\mathbf{x}, \boldsymbol{\theta})}{\mathbf{x}}}
#'   \item if \code{2} then \mjseqn{ \mathbf{h}\left(\mathbf{x}, \boldsymbol{\theta}\right) = \mathbf{x}}
#'   }
#' @param optimizer  optimization function for maximum likelihood estimation.
#' @param optim_method maximisation method that will be passed to [stats::optim()] function. Default is `BFGS`.
#' @param maxlik_method maximisation method that will be passed to [maxLik::maxLik()] function. Default is `NR`.
#' @param epsilon Tolerance for fitting algorithms by default \code{1e-6}.
#' @param maxit Maximum number of iterations.
#' @param trace logical value. If `TRUE` trace steps of the fitting algorithms. Default is `FALSE`
#' @param penalty The penalization function used during variables selection.
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
#' @param nleqslv_method The method that will be passed to [nleqslv::nleqslv()] function.
#' @param nleqslv_global The global strategy that will be passed to [nleqslv::nleqslv()] function.
#' @param nleqslv_xscalm The type of x scaling that will be passed to [nleqslv::nleqslv()] function.
#' @param dependence logical value (default `TRUE`) informing whether samples overlap (NOT YET IMPLEMENTED, FOR FUTURE DEVELOPMENT).
#' @param key binary key variable allowing to identify the overlap (NOT YET IMPLEMENTED, FOR FUTURE DEVELOPMENT).
#'
#' @return List with selected parameters.
#'
#' @seealso
#'
#' [nonprob()] -- for fitting procedure with non-probability samples.
#'
#' @export

control_sel <- function(est_method = c("mle", "gee"),
                        gee_h_fun = 1,
                        optimizer = c("maxLik", "optim"),
                        maxlik_method = "NR",
                        optim_method = "BFGS",
                        epsilon = 1e-4,
                        maxit = 500,
                        trace = FALSE,
                        penalty = c("SCAD", "lasso", "MCP"),
                        a_SCAD = 3.7,
                        a_MCP = 3,
                        lambda = -1,
                        lambda_min = .001,
                        nlambda = 50,
                        nfolds = 10,
                        print_level = 0,
                        start_type = c("glm", "naive", "zero"),
                        nleqslv_method = c("Broyden", "Newton"),
                        nleqslv_global = c("dbldog", "pwldog", "cline",
                                           "qline", "gline", "hook", "none"),
                        nleqslv_xscalm = c("fixed", "auto"),
                        dependence = FALSE,
                        key = NULL) {

  # Input validation
  optimizer <- match.arg(optimizer)
  est_method <- match.arg(est_method)
  penalty <- match.arg(penalty)
  start_type <- match.arg(start_type)
  nleqslv_method <- match.arg(nleqslv_method)
  nleqslv_global <- match.arg(nleqslv_global)
  nleqslv_xscalm <- match.arg(nleqslv_xscalm)

  if (!is.numeric(epsilon) || epsilon <= 0)
    stop("'epsilon' must be a positive number")

  if (!is.numeric(maxit) || maxit < 1 || maxit %% 1 != 0)
    stop("'maxit' must be a positive integer")

  if (!is.logical(trace))
    stop("'trace' must be logical")

  if (!is.character(maxlik_method))
    stop("'maxlik_method' must be a character string")

  if (!is.character(optim_method))
    stop("'optim_method' must be a character string")

  # not checked as not implemented

  # if (!is.logical(dependence))
  #   stop("'dependence' must be logical")
  #
  # if (!is.null(key) && !is.numeric(key))
  #   stop("'key' must be NULL or numeric")

  if (!is.numeric(gee_h_fun) || !gee_h_fun %in% c(1, 2))
    stop("'gee_h_fun' must be either 1 or 2")

  if (!is.numeric(a_SCAD) || a_SCAD <= 2)
    stop("'a_SCAD' must be greater than 2")

  if (!is.numeric(a_MCP) || a_MCP <= 1)
    stop("'a_MCP' must be greater than 1")

  if (!is.numeric(lambda))
    stop("'lambda' must be numeric")

  if (!is.numeric(lambda_min) || lambda_min <= 0 || lambda_min >= 1)
    stop("'lambda_min' must be between 0 and 1")

  if (!is.numeric(nlambda) || nlambda < 1 || nlambda %% 1 != 0)
    stop("'nlambda' must be a positive integer")

  if (!is.numeric(nfolds) || nfolds < 2 || nfolds %% 1 != 0)
    stop("'nfolds' must be an integer >= 2")

  if (!is.numeric(print_level) || print_level %% 1 != 0)
    stop("'print_level' must be an integer")

  list(
    est_method = est_method,
    gee_h_fun = gee_h_fun,
    optimizer = optimizer,
    maxlik_method = maxlik_method,
    optim_method = optim_method,
    epsilon = epsilon,
    maxit = maxit,
    trace = trace,
    penalty = penalty,
    a_SCAD = a_SCAD,
    a_MCP = a_MCP,
    lambda = lambda,
    lambda_min = lambda_min,
    nlambda = nlambda,
    nfolds = nfolds,
    print_level = print_level,
    start_type = start_type,
    nleqslv_method = nleqslv_method,
    nleqslv_global = nleqslv_global,
    nleqslv_xscalm = nleqslv_xscalm,
    dependence = dependence,
    key = key
  )
}
