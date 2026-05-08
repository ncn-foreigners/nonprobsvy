#' @title Control Parameters for Inference
#'
#' @description \code{control_inf} constructs a `list` with all necessary control parameters
#' for statistical inference.

#' @param var_method the variance method (default `"analytic"`).
#' @param rep_type the replication type for weights in the bootstrap method for variance estimation passed to [survey::as.svrepdesign()].
#'  Default is `"subbootstrap"`.
#' @param vars_selection logical scalar (default `FALSE`); if `TRUE`, then the variables selection model is used.
#' @param vars_combine logical scalar indicating whether variables should be combined after variable selection for doubly robust estimators (default `FALSE`)
#' @param bias_correction logical scalar (default `FALSE`); if `TRUE`, then the bias minimization estimation used during model fitting.
#' @param num_boot the number of iteration for bootstrap algorithms.
#' @param alpha significance level (default 0.05).
#' @param cores the number of cores in parallel computing (default 1).
#' @param keep_boot a logical value indicating whether statistics from bootstrap should be kept (default `TRUE`)
#' @param nn_exact_se a logical value indicating whether to compute the exact
#' standard error estimate for \code{nn} or \code{pmm} estimator. The variance estimator for
#' estimation based on \code{nn} or \code{pmm} can be decomposed into three parts, with the
#' third computed using covariance between imputed values for units in
#' the probability sample using predictive matches from the non-probability sample.
#' In most situations this term is negligible and is very computationally
#' expensive so by default it is set to \code{FALSE}, but the recommended option is to
#' set this value to \code{TRUE} before submitting the final results.
#'
#'
#' @return A `list` with selected parameters.
#'
#' @seealso
#'
#' [nonprob()] -- for fitting procedure with non-probability samples.
#'
#' @export

control_inf <- function(var_method = c("analytic", "bootstrap"),
                        rep_type = c("subbootstrap", "auto", "JK1", "JKn", "BRR", "bootstrap",
                                     "mrbbootstrap", "Fay"),
                        vars_selection = FALSE,
                        vars_combine = FALSE,
                        bias_correction = FALSE,
                        num_boot = 500,
                        alpha = 0.05,
                        cores = 1,
                        keep_boot = TRUE,
                        nn_exact_se = FALSE) {

  # Input validation
  var_method <- match.arg(var_method)
  rep_type <- match.arg(rep_type)

  is_logical_scalar <- function(x) is.logical(x) && length(x) == 1 && !is.na(x)
  is_numeric_scalar <- function(x) is.numeric(x) && length(x) == 1 && !is.na(x) && is.finite(x)

  if (!is_logical_scalar(vars_selection)) {
    stop("'vars_selection' must be a logical scalar")
  }

  if (!is_logical_scalar(vars_combine)) {
    stop("'vars_combine' must be a logical scalar")
  }

  if (!is_logical_scalar(bias_correction)) {
    stop("'bias_correction' must be a logical scalar")
  }

  if (!is_logical_scalar(keep_boot)) {
    stop("'keep_boot' must be a logical scalar")
  }

  if (!is_logical_scalar(nn_exact_se)) {
    stop("'nn_exact_se' must be a logical scalar")
  }

  if (!is_numeric_scalar(num_boot) || num_boot < 1 || num_boot %% 1 != 0) {
    stop("'num_boot' must be a positive integer")
  }

  if (!is_numeric_scalar(alpha) || alpha <= 0 || alpha >= 1) {
    stop("'alpha' must be between 0 and 1")
  }

  if (!is_numeric_scalar(cores) || cores < 1 || cores %% 1 != 0) {
    stop("'cores' must be a positive integer")
  }

  list(
    var_method = var_method,
    rep_type = rep_type,
    vars_selection = vars_selection,
    vars_combine = vars_combine,
    bias_correction = bias_correction,
    num_boot = num_boot,
    alpha = alpha,
    cores = cores,
    keep_boot = keep_boot,
    nn_exact_se = nn_exact_se
  )
}
