#' nonprob
#
#' nonprob: Function for inference based on nonprobability surveys using various methods
#
#' @param selection - `formula`, the selection (propensity) equation.
#' @param outcome - `formula`, the outcome equation.
#' @param target - `formula` with target variables.
#' @param data - an optional `data.frame` with data from the nonprobability sample.
#' @param svydesign - an optional `svydesign` object (from the survey package) containing probability sample.
#' @param pop_totals - an optional `named vector` with population totals.
#' @param pop_means - an optional `named vector` with population means.
#' @param pop_size - an optional `double` with population size.
#' @param method_selection - a `character` with method for propensity scores estimation
#' @param method_outcome - a `character` with method for response variable estimation
#' @param family_selection - a `character` string describing the error distribution and link function to be used in the model. Default is "binomial". Currently only binomial with logit link is supported.
#' @param family_outcome - a `character` string describing the error distribution and link function to be used in the model. Default is "gaussian". Currently supports: gaussian with identity link, poisson and binomial.
#' @param subset - an optional `vector` specifying a subset of observations to be used in the fitting process.
#' @param strata - an optional `vector` specifying strata.
#' @param weights - an optional `vector` of ‘prior weights’ to be used in the fitting process. Should be NULL or a numeric vector. It is assumed that this vector contains frequency or analytic weights
#' @param na_action - a function which indicates what should happen when the data contain `NAs`.
#' @param control_selection a list indicating parameters to use in fitting selection model for propensity scores
#' @param control_outcome a list indicating parameters to use in fitting model for outcome variable
#' @param control_inference a list indicating parameters to use in inference based on probablity and nonprobability samples, contains parameters such as estimation method or variance method
#' @param start - an optional `list` with starting values for the parameters of the selection and outcome equation
#' @param verbose - verbose, numeric
#' @param contrasts a
#' @param model a
#' @param x a
#' @param y a
#' @param ... a
#'
#' @export


nonprob <- function(selection = NULL,
                    outcome = NULL,
                    target = NULL,
                    data = NULL,
                    svydesign = NULL,
                    pop_totals = NULL,
                    pop_means = NULL,
                    pop_size = NULL,
                    method_selection = c("logit", "cloglog", "probit"),
                    method_outcome = c("glm.fit", "nn"),
                    family_selection = "binomial",
                    family_outcome = c("gaussian", "binomial", "poisson"),
                    subset,
                    strata,
                    weights = NULL,
                    na_action,
                    control_selection = controlSel(),
                    control_outcome = controlOut(),
                    control_inference = controlInf(est_method = "likelihood"),
                    start = NULL,
                    verbose = 0L,
                    contrasts = NULL,
                    model = TRUE,
                    x = TRUE,
                    y = TRUE,
                    ...) {
  est_method <- control_inference$est_method

  if (!is.data.frame(data)) {
    data <- data.frame(data)
  }

  if (is.null(weights)) weights <- rep.int(1, nrow(data))

  if(missing(method_selection)) method_selection <- "logit"
  if(missing(family_outcome)) family_outcome <- "gaussian"

  if(!(method_selection %in% c("logit", "cloglog", "probit"))) stop("Invalid method for selection formula.")
  if(!(family_outcome %in% c("gaussian", "binomial", "poisson"))) stop("Invalid family for outcome formula.")

  ## basic checkers
  if (is.null(selection) & is.null(outcome)) {
    stop("Please provide selection or outcome formula.")
  }
  if (inherits(selection, "formula") && inherits(target, "formula") && (is.null(outcome) || inherits(outcome, "formula") == FALSE)) {
    model_used <- "P"
  }

  if (inherits(outcome, "formula") && (is.null(selection) || inherits(selection, "formula") == FALSE)) {
    model_used <- "M"
  }

  if (inherits(selection, "formula") && inherits(outcome, "formula")) {

    ifelse(est_method == "likelihood", model_used <- "DR", model_used <- "DRsel")
  }

  ## validate data

  ## model estimates
  model_estimates <- switch(model_used,
    P = nonprobIPW(selection,
                   target,
                   data,
                   svydesign,
                   pop_totals,
                   pop_means,
                   pop_size,
                   method_selection,
                   family_selection,
                   subset,
                   strata,
                   weights,
                   na_action,
                   control_selection,
                   control_inference,
                   start,
                   verbose,
                   contrasts,
                   model,
                   x,
                   y,
                   ...),
    M = nonprobMI(outcome,
                  data,
                  svydesign,
                  method_outcome,
                  family_outcome,
                  subset,
                  strata,
                  weights,
                  na_action,
                  control_outcome,
                  control_inference,
                  start,
                  verbose,
                  contrasts,
                  model,
                  x,
                  y,
                  ...),
    DR = nonprobDR(selection,
                   outcome,
                   data,
                   svydesign,
                   pop_totals,
                   pop_means,
                   pop_size,
                   method_selection,
                   method_outcome,
                   family_selection,
                   family_outcome,
                   subset,
                   strata,
                   weights,
                   na_action,
                   control_selection,
                   control_outcome,
                   control_inference,
                   start,
                   verbose,
                   contrasts,
                   model,
                   x,
                   y,
                   ...),
    DRsel = nonprobSel(selection,
                        outcome,
                        data,
                        svydesign,
                        pop_totals,
                        pop_means,
                        pop_size,
                        method_selection,
                        method_outcome,
                        family_selection,
                        family_outcome,
                        subset,
                        strata,
                        weights,
                        na_action,
                        control_selection,
                        control_outcome,
                        control_inference,
                        start,
                        verbose,
                        contrasts,
                        model,
                        x,
                        y,
                        ...)
  )

  model_estimates
}


