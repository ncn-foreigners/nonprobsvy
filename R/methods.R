
# my methods --------------------------------------------------------------

#' @method extract nonprob
#' @exportS3Method extract nonprob
extract.nonprob <- function(object, what=c("all", "mean", "se")) {

  what_selected <- match.arg(what)

  ext <- switch(what_selected,
                "all" = cbind(target = rownames(object$output),
                              object$output,
                              object$confidence_interval),
                "mean" = data.frame(target = rownames(object$output),
                                    mean=object$output$mean),
                "se" = data.frame(target = rownames(object$output),
                                  SE=object$output$SE))

  rownames(ext) <- NULL

  return(ext)
}
#' @title Extracts Estimates from the Nonprob Class Object
#' @description Returns a \code{data.frame} of estimated mean(s) or standard error(s)
#' @param object object of of the \code{nonprob} class
#' @param what what to extract: all estimates (mean(s), SE(s) and CI(s); \code{"all"}; default), estimated mean(s) (\code{"mean"}) or their standard error(s) (\code{"se"})
#' @return a \code{data.frame} with selected information
#' @examples
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
#' extract(ipw_est1)
#' extract(ipw_est1, "se")
#' @export
extract <- function(object, what) {
  UseMethod("extract")
}

#' @method pop_size nonprob
#' @exportS3Method
pop_size.nonprob <- function(object) {
  object$pop_size
}

#' @title Returns Population Size (Estimated or Fixed)
#' @description Returns population size that is assumed to be
#'\itemize{
#'  \item{\code{fixed} -- if it is based on the `pop_size` argument,}
#'  \item{\code{estimated} -- if it is based on the probability survey specified in the `svydesign` or based on the estimated propensity scores for the non-probability sample.}
#'}
#' @param object object returned by the `nonprob` function.
#' @return a scalar returning the value of the population size.
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
#'
#' ipw_est2 <- nonprob(
#' selection = ~ region + private + nace + size,
#' target = ~ single_shift,
#' svydesign = jvs_svy,
#' data = admin, method_selection = "logit",
#' control_selection = control_sel(est_method = "gee", gee_h_fun = 1))
#'
#' ## estimated population size based on the non-calibrated IPW (MLE)
#' pop_size(ipw_est1)
#'
#' ## estimated population size based on the calibrated IPW (GEE)
#' pop_size(ipw_est2)
#'
#'
#' @export
pop_size <- function(object) {
  UseMethod("pop_size")
}


# base R methods ----------------------------------------------------------

#' @title Returns the Number of Rows in Samples
#' @description
#' Returns information on the number of rows of the probability sample (if provided)
#' and non-probability sample.
#' @param object a \code{nonprob} class object
#' @param ... other arguments passed to methods (currently not supported)
#'
#' @return a named \code{vector} with row numbers
#' @method nobs nonprob
#' @importFrom stats nobs
#' @exportS3Method
nobs.nonprob <- function(object,
                         ...) {
  c("prob" = object$prob_size, "nonprob" = object$nonprob_size)
}

#' @title Extracts the Inverse Probability (Propensity Score) Weights
#' @description A generic function `weights` that returns inverse probability weights (if present)
#'
#' @param object a \code{nonprob} class object
#' @param ... other arguments passed to methods (currently not supported)
#'
#' @returns A vector of weights or a `NULL` extracted from the `nonprob` object i.e. element `"ipw_weights"`
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
#' data = admin, method_selection = "logit", se = FALSE
#' )
#'
#' summary(weights(ipw_est1))
#'
#' @method weights nonprob
#' @importFrom stats weights
#' @exportS3Method
weights.nonprob <- function(object,
                            ...) {
  object$ipw_weights
}

#' @title Update Method for the Nonprob Object with Changed Arguments or Parameters
#' @author Maciej Beręsewicz
#'
#' @description
#'
#' The `update` method for the `nonprob` class object that allows to re-estimate
#' a given model with changed parameters. This is in particular useful if a user
#' would like to change method or estimate standard errors if they were not
#' estimated in the first place.
#'
#' @param object the `nonprob` class object
#' @param ... arguments passed to the `nonprob` class object
#' @param evaluate If true evaluate the new call else return the call
#'
#' @method update nonprob
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
#' data = admin, method_selection = "logit", se = FALSE
#' )
#'
#' ipw_est1
#'
#' update(ipw_est1, se = TRUE)
#'
#' @return returns a `nonprob` object
#' @importFrom stats update
#' @exportS3Method
update.nonprob <- function(object, ..., evaluate=TRUE) {
  call <- object$call

  # Handle additional arguments
  extras <- match.call(expand.dots = FALSE)$...
  if (length(extras)) {
    existing <- !is.na(match(names(extras), names(call)))
    for (a in names(extras)[existing]) call[[a]] <- extras[[a]]
    if (any(!existing)) {
      call <- c(as.list(call), extras[!existing])
      call <- as.call(call)
    }
  }

  # Either evaluate the call and return a new model,
  # or return the call itself
  if (evaluate) {
    eval(call, parent.frame())
  } else {
    call
  }
}

#' @title Returns Confidence Intervals for Estimated Mean
#'
#' @description A generic function that returns the confidence interval
#' for the estimated mean. If standard errors have not been estimated, the function
#' updates the `nonprob` object to obtain standard errors.
#'
#' @param object object of `nonprob` class.
#' @param parm names of parameters for which confidence intervals are to be
#' computed, if missing all parameters will be considered.
#' @param level confidence level for intervals.
#' @param ... additional arguments
#'
#' @examples
#' data(admin)
#' data(jvs)
#'
#' jvs_svy <- svydesign(ids = ~ 1,  weights = ~ weight,
#' strata = ~ size + nace + region, data = jvs)
#'
#' ipw_est1 <- nonprob(selection = ~ region + private + nace + size,
#' target = ~ single_shift,
#' svydesign = jvs_svy,
#' data = admin, method_selection = "logit", se = FALSE
#' )
#'
#' confint(ipw_est1)
#'
#' @method confint nonprob
#' @return returns a `data.frame` with confidence intervals for the target variables
#' @importFrom stats confint
#' @importFrom stats quantile
#' @exportS3Method
confint.nonprob <- function(object,
                            parm,
                            level = 0.95,
                            ...) {

  call <- object$call

  if ("se" %in% names(call)) {
    if (!eval(call$se))  {
      message("Calculating standard errors...")
      object <- update(object, se = T)
    }
  }

  if (missing(parm)) parm <- rownames(object$output)

  if (level == 0.95) {
    CIs <- object$confidence_interval
    CIs$target <- rownames(CIs)
    rownames(CIs) <- NULL
  } else {
    if (is.null(object$boot_sample)) {
      CIs <- object$output
      z <- stats::qnorm(1 - (1-level) / 2)
      # confidence interval based on the normal approximation
      CIs$lower_bound <- CIs$mean - z * CIs$SE
      CIs$upper_bound <- CIs$mean + z * CIs$SE
      CIs$target <- rownames(CIs)
      rownames(CIs) <- NULL
    } else {
      CIs <- object$output
      alpha <- 1-level
      SE_q <- apply(object$boot_sample, 2, stats::quantile, probs = c(alpha/2, 1-alpha/2))
      CIs$lower_bound <- SE_q[1,]
      CIs$upper_bound <- SE_q[2,]
      CIs$target <- rownames(CIs)
      rownames(CIs) <- NULL
    }
  }
  return(CIs[CIs$target %in% parm, c("target", "lower_bound", "upper_bound")])
}


#' @title Returns Coefficients of the Underlying Models
#' @description
#' Returns a \code{list} of coefficients for the selection and the outcome models
#'
#' @param object a \code{nonprob} class object
#' @param ... other arguments passed to methods (currently not supported)
#'
#' @return a \code{list} with two entries:
#' \itemize{
#' \item{\code{"coef_sel"} a matrix of coefficients for the selection equation if possible, else NULL}
#' \item{\code{"coef_dr"} a matrix of coefficients for the outcome equation(s) if possible, else NULL}
#' }
#'
#' @method coef nonprob
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
#' data = admin, method_selection = "logit", se = FALSE
#' )
#'
#' coef(ipw_est1)
#'
#' @importFrom stats coef
#' @exportS3Method
coef.nonprob <- function(object,
                         ...) {

  if (is.null(object$selection)) {
    coef_sel <- NULL
  } else if (isTRUE(object$control$control_inference$bias_correction)) {
    coef_sel <- sapply(object$selection,"[[", "x")
  } else {
    coef_sel <- as.matrix(coef(object$selection))
  }

  if (is.null(object$outcome) || grepl("^(nn|npar)", object$estimator_method)) {
    coef_out <- NULL
  } else if (isTRUE(object$control$control_inference$bias_correction)) {
    coef_out <- sapply(object$outcome, "[[", "x")
  } else {
    coef_out <- as.matrix(sapply(object$outcome, coef))
  }

  return(list(coef_sel=coef_sel, coef_out = coef_out))
}
