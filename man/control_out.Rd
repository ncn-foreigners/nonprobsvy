% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/control_outcome.R
\name{control_out}
\alias{control_out}
\title{Control parameters for outcome model}
\usage{
control_out(
  epsilon = 1e-08,
  maxit = 100,
  trace = FALSE,
  k = 5,
  penalty = c("SCAD", "lasso", "MCP"),
  a_SCAD = 3.7,
  a_MCP = 3,
  lambda_min = 0.001,
  nlambda = 100,
  nfolds = 10,
  treetype = c("kd", "rp", "ball"),
  searchtype = c("standard", "priority"),
  pmm_match_type = 1,
  pmm_weights = c("none", "dist"),
  pmm_k_choice = c("none", "min_var"),
  pmm_reg_engine = c("glm", "loess"),
  npar_loess = stats::loess.control(surface = "direct", trace.hat = "approximate")
)
}
\arguments{
\item{epsilon}{Tolerance for fitting algorithms. Default is \code{1e-6}.}

\item{maxit}{Maximum number of iterations.}

\item{trace}{logical value. If \code{TRUE} trace steps of the fitting algorithms. Default is \code{FALSE}.}

\item{k}{The k parameter in the \code{\link[RANN:nn2]{RANN::nn2()}} function. Default is 5.}

\item{penalty}{penalty algorithm for variable selection. Default is \code{SCAD}}

\item{a_SCAD}{The tuning parameter of the SCAD penalty for outcome model. Default is 3.7.}

\item{a_MCP}{The tuning parameter of the MCP penalty for outcome model. Default is 3.}

\item{lambda_min}{The smallest value for lambda, as a fraction of lambda.max. Default is .001.}

\item{nlambda}{The number of lambda values. Default is 100.}

\item{nfolds}{The number of folds during cross-validation for variables selection model.}

\item{treetype}{Type of tree for nearest neighbour imputation (for the NN and PMM estimator) passed to \code{\link[RANN:nn2]{RANN::nn2()}} function.}

\item{searchtype}{Type of search for nearest neighbour imputation (for the NN and PMM estimator) passed to \code{\link[RANN:nn2]{RANN::nn2()}} function.}

\item{pmm_match_type}{(Only for the PMM Estimator)
Indicates how to select 'closest' unit from non-probability sample for each
unit in probability sample. Either \code{1} (default) or \code{2} where
\code{2} is matching by minimizing distance between \eqn{\hat{y}_{i}} for
\eqn{i \in S_{A}} and \eqn{y_{j}} for \eqn{j \in S_{B}} and \code{1}
is matching by minimizing distance between \eqn{\hat{y}_{i}} for
\eqn{i \in S_{A}} and \eqn{\hat{y}_{i}} for \eqn{i \in S_{A}}.}

\item{pmm_weights}{(Only for the PMM Estimator)
Indicate how to weight \code{k} nearest neighbours in \eqn{S_{B}} to
create imputed value for units in \eqn{S_{A}}. The default value
\code{"none"} indicates that mean of \code{k} nearest \eqn{y}'s from
\eqn{S_{B}} should be used whereas \code{"prop_dist"} results in
weighted mean of these \code{k} values where weights are inversely
proportional to distance between matched values.}

\item{pmm_k_choice}{(Only for the PMM Estimator) Character value indicating how \code{k} hyper-parameter
should be chosen, by default \code{"none"} meaning \code{k} provided in
\code{control_outcome} argument will be used. For now the only other
option \code{"min_var"} means that \code{k}  will be chosen by minimizing
estimated variance of estimator for mean. Parameter \code{k} provided in
this control list will be chosen as starting point.}

\item{pmm_reg_engine}{(Only for the PMM Estimator) whether to use parametric (\code{"glm"})
or non-parametric (\code{"loess"}) regression model for the outcome. The default is \code{"glm"}.}

\item{npar_loess}{control parameters for the \link[stats:loess]{stats::loess} via the \link[stats:loess.control]{stats::loess.control} function.}
}
\value{
List with selected parameters.
}
\description{
\code{control_out} constructs a list with all necessary control parameters
for outcome model.
}
\seealso{
\code{\link[=nonprob]{nonprob()}} -- for fitting procedure with non-probability samples.
}
