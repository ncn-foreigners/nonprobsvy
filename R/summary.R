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
#' \item {\code{call} call}
#' \item {\code{estimator} type of estimator}
#' \item {\code{control} list of controls}
#' \item {\code{ipw_weights} estimated IPW weights}
#' \item {\code{ipw_weights_total} estimated IPW total (sum)}
#' \item {\code{ps_scores_nonprob} estimated propensity scores for non-probability sample}
#' \item {\code{ps_scores_prob} estimated propensity scores for probability sample}
#' \item {\code{case_weights} case weights}
#' \item {\code{output} estimated means and standard errors}
#' \item {\code{SE} estimated standard errors of V1 and V2}
#' \item {\code{confidence_interval} confidence intervals}
#' \item {\code{nonprob_size} size of the non-probability sample}
#' \item {\code{prob_size} size of the probability sample}
#' \item {\code{pop_size} population size}
#' \item {\code{pop_size_fixed} whether the population size is treated as fixed}
#' \item {\code{no_prob} whether probability sample was provided}
#' \item {\code{outcome} model details}
#' \item {\code{selection} selection details}
#' \item {\code{estimator_method} estimator method}
#' \item {\code{selection_formula} selection formula}
#' \item {\code{outcome_formula} outcome formula}
#' \item {\code{vars_selection} whether variable selection algorithm was applied}
#' \item {\code{vars_outcome} variables of the outcome models}
#' \item {\code{ys_rand_pred} predicted values for the random sample (if applies)}
#' \item {\code{ys_nons_pred} predicted values for the non-probability sample}
#' \item {\code{ys_resid} residuals for the non-probability sample}
#' }
#'
#' @examples
#'
#' data(admin)
#' data(jvs)
#'
#' jvs_svy <- svydesign(ids = ~ 1,  weights = ~ weight,
#' strata = ~ size + nace + region, data = jvs)
#'
#' ipw_est1 <- nonprob(selection = ~ region + private + nace + size,
#' target = ~ single_shift,
#' svydesign = jvs_svy,
#' data = admin, method_selection = "logit"
#' )
#' summary(ipw_est1)
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
         ipw_weights_total = if (object$estimator == "mi") NULL else sum(object$ipw_weights),
         ps_scores_nonprob = if (object$estimator == "mi") NULL else summary(object$ps_scores[object$R == 1]),
         ps_scores_prob = if (object$estimator == "mi" | object$prob_size == 0) NULL else summary(object$ps_scores[object$R == 0]),
         case_weights = summary(object$case_weights),
         output = object$output,
         SE = object$SE,
         confidence_interval = object$confidence_interval,
         nonprob_size = object$nonprob_size,
         prob_size = object$prob_size,
         pop_size = object$pop_size,
         pop_size_fixed = object$pop_size_fixed,
         no_prob = is.null(object$svydesign),
         outcome = object$outcome,
         selection = object$selection,
         estimator_method = object$estimator_method,
         selection_formula = object$selection_formula,
         outcome_formula = object$outcome_formula,
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





