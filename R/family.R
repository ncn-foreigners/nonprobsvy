# family

poisson_nonprobsvy <- function(link = "log") {
  mu <- function(eta) exp(eta)
  variance <- function(mu, y = NULL) mu
  mu_der <- function(mu) mu # first derivative
  mu_der2 <- function(mu) 1 # second derivative
  residuals <- function(mu, y) as.vector(y - mu)

  structure(list(mu = mu,
                 variance = variance,
                 mu_der = mu_der,
                 residuals = residuals),
            class = "family")
}

gaussian_nonprobsvy <- function(link = "identity") {
  mu <- function(eta) eta
  variance <- function(mu, y = NULL) mean((y - mu)^2) #rep.int(1, length(mu))
  mu_der <- function(mu) 1 # first derivative
  mu_der2 <- function(mu) 0 # second derivative
  residuals <- function(mu, y) as.vector(y - mu)

  structure(list(mu = mu,
                 variance = variance,
                 mu_der = mu_der,
                 residuals = residuals),
             class = "family")
}

binomial_nonprobsvy <- function(link = "logit") {
  link <- get_method(link)
  mu <- function(eta) link$make_link_inv(eta)
  variance <- function(mu, y = NULL) mu*(1-mu)
  mu_der <- function(mu) mu * (1 - mu) # first derivative
  mu_der2 <- function(mu) 1 - 2 * mu # second derivative
  residuals <- function(mu, y) as.vector(y - mu)

  structure(list(mu = mu,
                 variance = variance,
                 mu_der = mu_der,
                 residuals = residuals),
            class = "family")
}
