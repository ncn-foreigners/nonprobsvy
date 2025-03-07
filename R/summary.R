#' @title Summary statistics for model of the `nonprob` class.
#'
#' @description
#' Summarises the `nonprob` class object. The summary depends on the type of
#' the estimator (i.e. IPW, MI, DR)
#'
#'
#' @param object object of the `nonprob` class
#' @param ... Additional optional arguments
#'
#' @return An object of \code{nonprob_summary} class containing:
#' \itemize{
#' \item \code{call}
#' }
#'
#'
#' @method summary nonprob
#' @exportS3Method
summary.nonprob <- function(object, ...) {

  if (object$estimator != "ipw") {
    summary_ys_rand_pred <- lapply(object$ys_rand_pred, summary)
    summary_ys_nons_pred <- lapply(object$ys_nons_pred, summary)
    summary_ys_resid <- lapply(object$ys_resid, summary)
    names(summary_ys_rand_pred) <- names(summary_ys_nons_pred) <- names(summary_ys_resid) <- names(object$y)
  } else {
    summary_ys_rand_pred <- object$ys_rand_pred
    summary_ys_nons_pred <- object$ys_nons_pred
    summary_ys_resid <- object$ys_resid
  }


  res <- structure(
    list(call = object$call,
         estimator = object$estimator,
         control = object$control,
         ipw_weights = if (object$estimator == "mi") NULL else summary(object$ipw_weights),
         ps_scores = if (object$estimator == "mi") NULL else summary(object$ps_scores),
         case_weights = summary(object$case_weights),
         output = object$output,
         SE = object$SE,
         confidence_interval = object$confidence_interval,
         nonprob_size = object$nonprob_size,
         prob_size = object$prob_size,
         pop_size = object$pop_size,
         pop_size_fixed = object$pop_size_fixed,
         outcome = object$outcome,
         selection = object$selection,
         estimator_method = object$estimator_method,
         selection_formula = object$selection_formula,
         outcome_formula = object$outcome_formula,
         outcome = object$outcome,
         selection = object$selection,
         vars_selection = names(object$selection$coefficients),
         vars_outcome = lapply(object$outcome, function(x) names(x$coefficients)),
         ys_rand_pred = summary_ys_rand_pred,
         ys_nons_pred = summary_ys_nons_pred,
         ys_resid = summary_ys_resid
         ),
    class = c("nonprob_summary")
  )
  res
}





