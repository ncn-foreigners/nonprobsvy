% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/control_selection.R
\name{controlSel}
\alias{controlSel}
\title{Control parameters for selection model}
\usage{
controlSel(
  method = "glm.fit",
  epsilon = 1e-04,
  maxit = 500,
  trace = FALSE,
  optimizer = c("maxLik", "optim"),
  maxLik_method = "NR",
  optim_method = "BFGS",
  dependence = FALSE,
  key = NULL,
  est_method_sel = c("mle", "gee"),
  h = c(1, 2),
  penalty = c("SCAD", "lasso", "MCP"),
  a_SCAD = 3.7,
  a_MCP = 3,
  lambda = -1,
  lambda_min = 0.001,
  nlambda = 50,
  nfolds = 10,
  print_level = 0,
  start_type = c("glm", "naive", "zero")
)
}
\arguments{
\item{method}{estimation method.}

\item{epsilon}{Tolerance for fitting algorithms by default \code{1e-6}.}

\item{maxit}{Maximum number of iterations.}

\item{trace}{logical value. If \code{TRUE} trace steps of the fitting algorithms. Default is \code{FALSE}}

\item{optimizer}{\itemize{
\item optimization function for maximum likelihood estimation.
}}

\item{maxLik_method}{maximisation method that will be passed to \code{\link[maxLik:maxLik]{maxLik::maxLik()}} function. Default is \code{NR}.}

\item{optim_method}{maximisation method that will be passed to \code{\link[stats:optim]{stats::optim()}} function. Default is \code{BFGS}.}

\item{dependence}{logical value - \code{TRUE} if samples are dependent.}

\item{key}{binary key variable}

\item{est_method_sel}{Method of estimation for propensity score model.}

\item{h}{Smooth function for the generalized estimating equations methods taking the following values
\itemize{
\item if \code{1} then \mjseqn{\mathbf{h}\left(\mathbf{x}, \boldsymbol{\theta}\right) =
  \frac{\pi(\mathbf{x}, \boldsymbol{\theta})}{\mathbf{x}}}
\item if \code{2} then \mjseqn{ \mathbf{h}\left(\mathbf{x}, \boldsymbol{\theta}\right) = \mathbf{x}}
}}

\item{penalty}{The penanlization function used during variables selection.}

\item{a_SCAD}{The tuning parameter of the SCAD penalty for selection model. Default is 3.7.}

\item{a_MCP}{The tuning parameter of the MCP penalty for selection model. Default is 3.}

\item{lambda}{A user-specified \mjseqn{\lambda} value during variable selection model fitting.}

\item{lambda_min}{The smallest value for lambda, as a fraction of \code{lambda.max}. Default is .001.}

\item{nlambda}{The number of \code{lambda} values. Default is 50.}

\item{nfolds}{The number of folds for cross validation. Default is 10.}

\item{print_level}{this argument determines the level of printing which is done during the optimization (for propensity score model) process.}

\item{start_type}{\itemize{
\item Type of method for start points for model fitting taking the following values
\itemize{
\item if \code{glm} then start taken from the glm function called on samples.
\item if \code{naive} then start consists of a vector which has the value of an estimated parameter for one-dimensional data (on intercept) and 0 for the rest.
\item if \code{zero} then start is a vector of zeros.
}
}}
}
\value{
List with selected parameters.
}
\description{
\code{controlSel} constructs a list with all necessary control parameters
for selection model.
}
\seealso{
\code{\link[=nonprob]{nonprob()}} -- for fitting procedure with non-probability samples.
}
\author{
Łukasz Chrostowski, Maciej Beręsewicz
\loadmathjax
}
