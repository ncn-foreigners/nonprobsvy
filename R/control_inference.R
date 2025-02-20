#' @title Control parameters for inference
#'
#' @description \code{control_inf} constructs a list with all necessary control parameters
#' for statistical inference.
#'
#' @param vars_selection If `TRUE`, then the variables selection model is used.
#' @param var_method the variance method.
#' @param rep_type the replication type for weights in the bootstrap method for variance estimation passed to [survey::as.svrepdesign()].
#'  Default is `subbootstrap`.
#' @param bias_inf the inference method in the bias minimization.
#' \itemize{
#'   \item if \code{union}, then the final model is fitted on the union of selected variables for selection and outcome models
#'   \item if \code{div}, then the final model is fitted separately on division of selected variables into relevant ones for
#'   selection and outcome model.
#'   }
#' @param bias_correction if `TRUE`, then the bias minimization estimation used during model fitting.
#' @param num_boot the number of iteration for bootstrap algorithms.
#' @param alpha significance level, 0.05 by defult.
#' @param cores the number of cores in parallel computing.
#' @param keep_boot a logical value indicating whether statistics from bootstrap should be kept.
#' By default set to \code{TRUE}
#' @param nn_exact_se a logical value indicating whether to compute the exact
#' standard error estimate for \code{nn} or \code{pmm} estimator. The variance estimator for
#' estimation based on \code{nn} or \code{pmm} can be decomposed into three parts, with the
#' third computed using covariance between imputed values for units in
#' the probability sample using predictive matches from the non-probability sample.
#' In most situations this term is negligible and is very computationally
#' expensive so by default it is set to \code{FALSE}, but the recommended option is to
#' set this value to \code{TRUE} before submitting the final results.
#' @param pi_ij TODO, either a matrix or a \code{ppsmat} class object.
#'
#'
#' @return List with selected parameters.
#'
#' @seealso
#'
#' [nonprob()] -- for fitting procedure with non-probability samples.
#'
#' @export

control_inf <- function(vars_selection = FALSE,
                       var_method = c(
                         "analytic",
                         "bootstrap"
                       ),
                       rep_type = c(
                         "auto", "JK1", "JKn", "BRR", "bootstrap",
                         "subbootstrap", "mrbbootstrap", "Fay"
                       ),
                       bias_inf = c("union", "div"),
                       num_boot = 500,
                       bias_correction = FALSE,
                       alpha = 0.05,
                       cores = 1,
                       keep_boot,
                       nn_exact_se = FALSE,
                       pi_ij) {
  list(
    vars_selection = if (missing(vars_selection)) FALSE else vars_selection,
    var_method = if (missing(var_method)) "analytic" else var_method,
    rep_type = if (missing(rep_type)) "subbootstrap" else rep_type,
    bias_inf = if (missing(bias_inf)) "union" else bias_inf,
    bias_correction = bias_correction,
    num_boot = num_boot,
    alpha = alpha,
    cores = cores,
    keep_boot = if (missing(keep_boot)) {
      TRUE
    } else {
      if (!is.logical(keep_boot)) {
        stop("keep_boot argument for controlInf must be logical")
      } else {
        keep_boot
      }
    },
    nn_exact_se = if (!is.logical(nn_exact_se) & length(nn_exact_se) == 1) {
      stop("Argument nn_exact_se must be a logical scalar")
    } else {
      nn_exact_se
    },
    pi_ij = if (missing(pi_ij)) NULL else pi_ij
  )
}
