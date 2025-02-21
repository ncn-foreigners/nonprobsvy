#' @title Control parameters for outcome model
#'
#' @importFrom stats loess.control
#' @description \code{control_out} constructs a list with all necessary control parameters
#' for outcome model.
#'
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
#' @param treetype Type of tree for nearest neighbour imputation (for the NN and PMM estimator) passed to [RANN::nn2()] function.
#' @param searchtype Type of search for nearest neighbour imputation (for the NN and PMM estimator) passed to [RANN::nn2()] function.
#' @param pmm_match_type (Only for the PMM Estimator)
#' Indicates how to select 'closest' unit from nonprobability sample for each
#' unit in probability sample. Either \code{1} (default) or \code{2} where
#' \code{2} is matching by minimizing distance between \mjseqn{\hat{y}_{i}} for
#' \mjseqn{i \in S_{A}} and \mjseqn{y_{j}} for \mjseqn{j \in S_{B}} and \code{1}
#' is matching by minimizing distance between \mjseqn{\hat{y}_{i}} for
#' \mjseqn{i \in S_{A}} and \mjseqn{\hat{y}_{i}} for \mjseqn{i \in S_{A}}.
#' @param pmm_weights (Only for the PMM Estimator)
#' Indicate how to weight \code{k} nearest neighbours in \mjseqn{S_{B}} to
#' create imputed value for units in \mjseqn{S_{A}}. The default value
#' \code{"none"} indicates that mean of \code{k} nearest \mjseqn{y}'s from
#' \mjseqn{S_{B}} should be used whereas \code{"prop_dist"} results in
#' weighted mean of these \code{k} values where weights are inversely
#' proportional to distance between matched values.
#' @param pmm_k_choice (Only for the PMM Estimator) Character value indicating how \code{k} hyper-parameter
#' should be chosen, by default \code{"none"} meaning \code{k} provided in
#' \code{control_outcome} argument will be used. For now the only other
#' option \code{"min_var"} means that \code{k}  will be chosen by minimizing
#' estimated variance of estimator for mean. Parameter \code{k} provided in
#' this control list will be chosen as starting point.
#' @param pmm_reg_engine (Only for the PMM Estimator) whether to use parametric (`"glm"`)
#' or non-parametric (`"loess"`) regression model for the outcome. The default is `"glm"`.
#' @param npar_loess control parameters for the [stats::loess] via the [stats::loess.control] function.
#'
#' @return List with selected parameters.
#'
#' @seealso
#'
#' [nonprob()] -- for fitting procedure with non-probability samples.
#'
#'
#' @export

control_out <- function(epsilon = 1e-4,
                        maxit = 100,
                        trace = FALSE,
                        k = 1,
                        penalty = c("SCAD", "lasso", "MCP"),
                        a_SCAD = 3.7,
                        a_MCP = 3,
                        lambda_min = .001,
                        nlambda = 100,
                        nfolds = 10,
                        treetype = c("kd", "rp", "ball"),
                        searchtype = c("standard", "priority"),
                        pmm_match_type = 1,
                        pmm_weights = c("none", "dist"),
                        pmm_k_choice = c("none", "min_var"),
                        pmm_reg_engine = c("glm", "loess"),
                        npar_loess = stats::loess.control(surface = "interpolate", trace.hat = "approximate")) {

  # Input validation
  penalty <- match.arg(penalty)
  treetype <- match.arg(treetype)
  searchtype <- match.arg(searchtype)
  pmm_weights <- match.arg(pmm_weights)
  pmm_k_choice <- match.arg(pmm_k_choice)
  pmm_reg_engine <- match.arg(pmm_reg_engine)

  if (!is.numeric(epsilon) || epsilon <= 0)
    stop("'epsilon' must be a positive number")

  if (!is.numeric(maxit) || maxit < 1 || maxit %% 1 != 0)
    stop("'maxit' must be a positive integer")

  if (!is.logical(trace))
    stop("'trace' must be logical")

  if (!is.numeric(k) || k < 1 || k %% 1 != 0)
    stop("'k' must be a positive integer")

  if (!is.numeric(a_SCAD) || a_SCAD <= 2)
    stop("'a_SCAD' must be greater than 2")

  if (!is.numeric(a_MCP) || a_MCP <= 1)
    stop("'a_MCP' must be greater than 1")

  if (!is.numeric(lambda_min) || lambda_min <= 0 || lambda_min >= 1)
    stop("'lambda_min' must be between 0 and 1")

  if (!is.numeric(nlambda) || nlambda < 1 || nlambda %% 1 != 0)
    stop("'nlambda' must be a positive integer")

  if (!is.numeric(nfolds) || nfolds < 2 || nfolds %% 1 != 0)
    stop("'nfolds' must be an integer >= 2")

  if (!pmm_match_type %in% c(1, 2))
    stop("'pmm_match_type' must be either 1 or 2")

  list(
    epsilon = epsilon,
    maxit = maxit,
    trace = trace,
    k = k,
    penalty = penalty,
    a_SCAD = a_SCAD,
    a_MCP = a_MCP,
    lambda_min = lambda_min,
    nlambda = nlambda,
    nfolds = nfolds,
    treetype = treetype,
    searchtype = searchtype,
    pmm_match_type = pmm_match_type,
    pmm_weights = pmm_weights,
    pmm_k_choice = pmm_k_choice,
    pmm_reg_engine = pmm_reg_engine,
    npar_loess = npar_loess
  )
}
