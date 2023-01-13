#' nonprobDR
#
#' nonprobDR: Function for inference based on nonprobability surveys, nonprobabilty big data sample and estimated propensoty scores.
#
#' @param selection - `formula`, the selection (propensity) equation.
#' @param outcome - `formula`, the outcome equation.
#' @param data - an optional `data.frame` with data from the nonprobability sample.
#' @param svydesign - an optional `svydesign` object (from the survey package) containing probability sample.
#' @param family.selection - a `character` string describing the error distribution and link function to be used in the model. Default is "binomial". Currently only binomial with logit link is supported.
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


nonprobDR <- function(selection,
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
                      ...){









}
