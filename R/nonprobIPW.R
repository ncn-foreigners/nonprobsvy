#' nonprobIPW
#
#' nonprobIPW: Function for inference based on nonprobability big data sample and estimated propensity scores.
#
#' @param selection - `formula`, the selection (propensity) equation.
#' @param outcome - `formula`, the outcome equation.
#' @param data - an optional `data.frame` with data from the nonprobability sample.
#' @param svydesign - an optional `svydesign` object (from the survey package) containing probability sample.
#' @param family.outcome - a `character` string describing the error distribution and link function to be used in the model. Default is "gaussian". Currently supports: gaussian with identity link, poisson and binomial.
#' @param subset - an optional `vector` specifying a subset of observations to be used in the fitting process.
#' @param weights - an optional `vector` of ‘prior weights’ to be used in the fitting process. Should be NULL or a numeric vector. It is assumed that this vector contains frequency or analytic weights
#' @param na.action
#' @param control.outcome
#' @param control.inference
#' @param start
#' @param verbose
#' @param contrasts
#' @param model
#' @param x
#' @param y
#' @param ...
#'
#' @export



nonprobIPW <- function(selection,
                       data,
                       svydesign,
                       pop.totals,
                       pop.means,
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
                       ...){







}
