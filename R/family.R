# family
## functions for the glm
## for extending existing functions

poisson_nonprobsvy <- function(link = "log") {
  x <- stats::poisson(link = link)
  # x$mu_der <- function(mu) mu # first derivative
  x$mu.eta2 <- function(mu) 1 # second derivative
  x$residuals <- function(mu, y) as.vector(y - mu)

  class(x) <- c(class(x), "nonprobsvy_family")
  x
}

gaussian_nonprobsvy <- function(link = "identity") {
  x <- stats::gaussian(link = link)
  # x$mu_der <- function(mu) 1
  x$residuals <- function(mu, y) as.vector(y - mu)

  class(x) <- c(class(x), "nonprobsvy_family")
  x
}

binomial_nonprobsvy <- function(link = "logit") {
  x <- stats::binomial(link = link)
  x$mu.eta2 <- function(mu) 1 - 2 * mu # second derivative
  x$residuals <- function(mu, y) as.vector(y - mu)

  class(x) <- c(class(x), "nonprobsvy_family")
  x
}
