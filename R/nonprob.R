#' nonprob
#
#' nonprob: Function for inference based on nonprobability surveys using various methods
#
#' @param selection - `formula`, the selection (propensity) equation.
#' @param outcome - `formula`, the outcome equation.
#' @param data - an optional `data.frame` with data from the nonprobability sample.
#' @param svydesign - an optional `svydesign` object (from the survey package) containing probability sample.
#' @param pop.totals - an optional `named vector` with population totals.
#' @param pop.means - an optional `named vector` with population means.
#' @param pop.size - an optional `double` with population size.
#' @param method.selection - a `character` with method for propensity scores estimation
#' @param method.outcome - a `character` with method for response variable estimation
#' @param family.selection - a `character` string describing the error distribution and link function to be used in the model. Default is "binomial". Currently only binomial with logit link is supported.
#' @param family.outcome - a `character` string describing the error distribution and link function to be used in the model. Default is "gaussian". Currently supports: gaussian with identity link, poisson and binomial.
#' @param subset - an optional `vector` specifying a subset of observations to be used in the fitting process.
#' @param weights - an optional `vector` of ‘prior weights’ to be used in the fitting process. Should be NULL or a numeric vector. It is assumed that this vector contains frequency or analytic weights
#' @param na.action - a function which indicates what should happen when the data contain `NAs`.
#' @param control.selection a
#' @param control.outcome a
#' @param control.inference a
#' @param start - an optional `list` with starting values for the parameters of the selection and outcome equation
#' @param verbose - verbose, numeric
#' @param contrasts a
#' @param model a
#' @param x a
#' @param y a
#' @param ... a
#' @export


nonprob <- function(selection = NULL,
                    outcome = NULL,
                    data = NULL,
                    svydesign = NULL,
                    pop.totals = NULL,
                    pop.means = NULL,
                    pop.size = NULL,
                    method.selection = c("logit", "cloglog", "probit"),
                    method.outcome = c("glm.fit", "nn"),
                    family.selection = "binomial",
                    family.outcome = "gaussian",
                    subset,
                    weights = NULL,
                    na.action,
                    control.selection = controlSel(),
                    control.outcome = controlOut(),
                    control.inference = controlInf(),
                    start = NULL,
                    verbose = 0L,
                    contrasts = NULL,
                    model = TRUE,
                    x = TRUE,
                    y = TRUE,
                    ...) {

  ##

  est.method <- control.inference$est.method

  # if (missing(method.selection)) method.selection <- "logit"

  if (!is.data.frame(data)) {
    data <- data.frame(data)
  }

  if (is.null(weights)) weights <- rep.int(1, nrow(data))

  ## basic checkers
  if (is.null(selection) & is.null(outcome)) {
    stop("Please provide selection or outcome formula.")
  }

  if (inherits(selection, "formula") & (is.null(outcome) | inherits(outcome, "formula") == FALSE)) {
    ifelse(is.null(pop.size), model_used <- "P1", "P2")
  }

  if (inherits(outcome, "formula") & (is.null(selection) | inherits(selection, "formula") == FALSE)) {
    model_used <- "M"
  }

  if (inherits(selection, "formula") & inherits(outcome, "formula")) {

    ifelse(est.method == "likelihood", model_used <- "DR", model_used <- "DRsel")
  }

  ## validate data

  ## model estimates
  model_estimates <- switch(model_used,
    P1 = nonprobIPW(selection,
                   data,
                   svydesign,
                   pop.totals,
                   pop.means,
                   pop.size,
                   method.selection,
                   family.selection = "binomial",
                   subset,
                   weights,
                   na.action,
                   control.selection = controlSel(),
                   control.inference = controlInf(),
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
                  method.outcome,
                  family.outcome = "gaussian",
                  subset,
                  weights,
                  na.action,
                  control.outcome = controlOut(),
                  control.inference = controlInf(),
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
                   pop.totals,
                   pop.means,
                   pop.size,
                   method.selection,
                   method.outcome,
                   family.selection = "binomial",
                   family.outcome = "gaussian",
                   subset,
                   weights,
                   na.action,
                   control.selection,
                   control.outcome,
                   control.inference,
                   start,
                   verbose,
                   contrasts,
                   model,
                   x,
                   y,
                   ...),
    DRsel <- nonprobSel(selection,
                        outcome,
                        data,
                        svydesign,
                        pop.totals,
                        pop.means,
                        pop.size,
                        method.selection,
                        method.outcome,
                        family.selection = "binomial",
                        family.outcome = "gaussian",
                        subset,
                        weights,
                        na.action,
                        control.selection,
                        control.outcome,
                        control.inference,
                        start,
                        verbose,
                        contrasts,
                        model,
                        x,
                        y,
                        ...)
  )

  return(model_estimates)

}
