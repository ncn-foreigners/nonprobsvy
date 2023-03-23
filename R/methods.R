# Internal functions, no need for documenting them


pearson_nonprobsvy <- function(X_nons,
                               X_rand,
                               ps_nons,
                               est_ps_rand) {

  Rnons <- c(rep(1, nrow(X_nons)), rep(0, nrow(X_rand)))
  ps <- c(ps_nons, est_ps_rand)

  res <- (Rnons - ps)/sqrt(ps*(1-ps))
  res
}


deviance_nonprobsvy <- function(X_nons,
                                X_rand,
                                ps_nons,
                                est_ps_rand) {

  Rnons <- c(rep(1, nrow(X_nons)), rep(0, nrow(X_rand)))
  s <- c(rep(1, nrow(X_nons)), rep(-1, nrow(X_rand)))
  ps <- c(ps_nons, est_ps_rand)

  res <- s * sqrt(-2 * (Rnons * log(ps) + (1 - Rnons) * log(1 - ps)))
  res
}

# example using
#pearson_residuals <- pearson_nonprobsvy(X_nons, X_rand, , est_ps_rand) # pearson residuals for propensity score model
#deviance_residuals <- deviance_nonprobsvy(X_nons, X_rand, ps_nons, est_ps_rand) # deviance residuals for propensity score model

cooks_distance_nonprobsvy <- function() {

}


hatvalues_nonprobsvy <- function() {


}


AIC_nonprobsvy <- function() {


}

BIC_nonprobsvy <- function() {



}
