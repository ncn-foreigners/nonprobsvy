#' @title Control parameters for outcome model
#'
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
#' @param treetype Type of tree for nearest neighbour imputation passed to [RANN::nn2()] function.
#' @param searchtype Type of search for nearest neighbour imputation passed to [RANN::nn2()] function.
#' @param predictive_match (Only for predictive mean matching)
#' Indicates how to select 'closest' unit from nonprobability sample for each
#' unit in probability sample. Either \code{1} (default) or \code{2} where
#' \code{2} is matching by minimizing distance between \mjseqn{\hat{y}_{i}} for
#' \mjseqn{i \in S_{A}} and \mjseqn{y_{j}} for \mjseqn{j \in S_{B}} and \code{1}
#' is matching by minimizing distance between \mjseqn{\hat{y}_{i}} for
#' \mjseqn{i \in S_{A}} and \mjseqn{\hat{y}_{i}} for \mjseqn{i \in S_{A}}.
#' @param pmm_weights (Only for predictive mean matching)
#' Indicate how to weight \code{k} nearest neighbours in \mjseqn{S_{B}} to
#' create imputed value for units in \mjseqn{S_{A}}. The default value
#' \code{"none"} indicates that mean of \code{k} nearest \mjseqn{y}'s from
#' \mjseqn{S_{B}} should be used whereas \code{"prop_dist"} results in
#' weighted mean of these \code{k} values where weights are inversely
#' proportional to distance between matched values.
#' @param pmm_k_choice Character value indicating how \code{k} hyper-parameter
#' should be chosen, by default \code{"none"} meaning \code{k} provided in
#' \code{control_outcome} argument will be used. For now the only other
#' option \code{"min_var"} means that \code{k}  will be chosen by minimizing
#' estimated variance of estimator for mean. Parameter \code{k} provided in
#' this control list will be chosen as starting point.
#' @param pmm_reg_engine TODO
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
                       treetype = "kd",
                       searchtype = "standard",
                       predictive_match = 1:2,
                       pmm_weights = c("none", "prop_dist"),
                       pmm_k_choice = c("none", "min_var"),
                       pmm_reg_engine = c("glm", "loess")) {
  if (missing(predictive_match)) {
    predictive_match <- 1
  }

  if (missing(pmm_weights)) {
    pmm_weights <- "none"
  }

  list(
    epsilon = epsilon,
    maxit = maxit,
    trace = trace,
    k = k,
    penalty = if (missing(penalty)) "SCAD" else penalty,
    a_SCAD = a_SCAD,
    a_MCP = a_MCP,
    lambda_min = lambda_min,
    nlambda = nlambda,
    nfolds = nfolds,
    treetype = treetype,
    searchtype = searchtype,
    predictive_match = predictive_match,
    pmm_weights = pmm_weights,
    pmm_k_choice = if (missing(pmm_k_choice)) "none" else pmm_k_choice,
    pmm_reg_engine = if (missing(pmm_reg_engine)) "glm" else pmm_reg_engine
  )
}
