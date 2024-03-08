u_beta_mi <- function(par,
                      R,
                      X,
                      y,
                      weights,
                      family_nonprobsvy) { # TODO

  if (is.character(family_nonprobsvy)) {
    family_nonprobsvy <- paste(family_nonprobsvy, "_nonprobsvy", sep = "")
    family_nonprobsvy <- get(family_nonprobsvy, mode = "function", envir = parent.frame())
    family_nonprobsvy <- family_nonprobsvy()
  }

  p <- ncol(X)
  beta <- par
  eta <- X %*% beta
  mu <- family_nonprobsvy$mu.eta(eta)
  # mu_der <- family_nonprobsvy$mu.eta2(mu)

  n <- length(R)
  R_rand <- 1 - R
  loc_nons <- which(R == 1)
  loc_rand <- which(R == 0)
  y_mean <- mean(y[loc_nons])

  # UTB <- apply(X * y - X * R_rand * weights * as.vector(mu), 2, sum)
  # UTB <- apply(X[loc_rand,] * (weights[loc_rand] * as.vector(mu[loc_rand]) - mean(y)), 2, sum)
  # UTB <- apply(X * (weights * as.vector(mu) - mean(y)), 2, sum)
  UTB <- apply(X * (R_rand * weights * as.vector(eta) - y_mean), 2, sum)
  UTB
}
