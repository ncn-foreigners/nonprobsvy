#' @import mathjaxr
NULL
#' @title Inference with the non-probability survey samples.
#' @author Łukasz Chrostowski, Maciej Beręsewicz
#'
#' @description \code{nonprob} fits model for inference based on non-probability surveys using various methods.
#'
#'
#' @param selection `formula`, the selection (propensity) equation.
#' @param outcome `formula`, the outcome equation.
#' @param target `formula` with target variables.
#' @param data an optional `data.frame` with data from the nonprobability sample.
#' @param svydesign an optional `svydesign` object (from the survey package) containing probability sample.
#' @param pop_totals an optional `named vector` with population totals.
#' @param pop_means an optional `named vector` with population means.
#' @param pop_size an optional `double` with population size.
#' @param method_selection a `character` with method for propensity scores estimation
#' @param method_outcome a `character` with method for response variable estimation
#' @param family_selection a `character` string describing the error distribution and link function to be used in the model. Default is "binomial". Currently only binomial with logit link is supported.
#' @param family_outcome a `character` string describing the error distribution and link function to be used in the model. Default is "gaussian". Currently supports: gaussian with identity link, poisson and binomial.
#' @param subset an optional `vector` specifying a subset of observations to be used in the fitting process.
#' @param strata an optional `vector` specifying strata.
#' @param weights an optional `vector` of ‘prior weights’ to be used in the fitting process. Should be NULL or a numeric vector. It is assumed that this vector contains frequency or analytic weights
#' @param na_action a function which indicates what should happen when the data contain `NAs`.
#' @param control_selection a list indicating parameters to use in fitting selection model for propensity scores
#' @param control_outcome a list indicating parameters to use in fitting model for outcome variable
#' @param control_inference a list indicating parameters to use in inference based on probability and non-probability samples, contains parameters such as estimation method or variance method
#' @param start an optional `list` with starting values for the parameters of the selection and outcome equation
#' @param verbose verbose, numeric
#' @param contrasts a
#' @param model a
#' @param x a
#' @param y a
#' @param ... Additional, optional arguments.
#'
#' @references
#'
#' Kim JK, Park S, Chen Y, Wu C. Combining non-probability and
#' probability survey samples through mass imputation. J R Stat Soc Series A. 2021;184:941–
#' 963.
#'
#' Shu Yang, Jae Kwang Kim, Rui Song. Doubly robust inference when combining probability
#' and non-probability samples with high dimensional data. J. R. Statist. Soc. B (2020)
#'
#' Yilin Chen , Pengfei Li & Changbao Wu (2020) Doubly Robust Inference
#' With Nonprobability Survey Samples, Journal of the American Statistical Association, 115:532,
#' 2011-2021
#'
#' Shu Yang, Jae Kwang Kim and Youngdeok Hwang Integration of data from
#' probability surveys and big found data for finite population inference using mass imputation.
#' Survey Methodology, June 2021 29 Vol. 47, No. 1, pp. 29-58
#'
#' @name main_doc
NULL
