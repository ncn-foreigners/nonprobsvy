# family

#' @export
poisson_nonprobsvy <- function(link = "log") {
  mu <- function(eta) exp(eta)
  variance <- function(mu, y = NULL) mu
  mu_der <- function(mu) mu
  residuals <- function(mu, y) as.vector(y - mu)

  structure(list(mu = mu,
                 variance = variance,
                 mu_der = mu_der,
                 residuals = residuals),
            class = "family")
}

#' @export
gaussian_nonprobsvy <- function(link = "identity") {
  mu <- function(eta) eta
  variance <- function(mu, y)  mean((y - mu)^2)
  mu_der <- function(mu) 1
  residuals <- function(mu, y) as.vector(y - mu)

  structure(list(mu = mu,
                 variance = variance,
                 mu_der = mu_der,
                 residuals = residuals),
             class = "family")
}

#' @export
binomial_nonprobsvy <- function(link = "logit") {
  link <- get_method(link)
  mu <- function(eta) link$make_link_inv(eta)
  variance <- function(mu, y = NULL) mu*(1-mu)
  mu_der <- function(mu) mu * (1 - mu)
  residuals <- function(mu, y) as.vector(y - mu)

  structure(list(mu = mu,
                 variance = variance,
                 mu_der = mu_der,
                 residuals = residuals),
            class = "family")
}
