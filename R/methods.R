#' Summary statistics for model
#'
#'
#' pearson.nonprobsvy
#'
#' pearson.nonprobsvy: Function for pearson residuals computing for propensity score model

pearson.nonprobsvy <- function(X_nons,
                               X_rand,
                               ps_nons,
                               est_ps_rand){

  Rnons <- c(rep(1, nrow(X_nons)), rep(0, nrow(X_rand)))
  ps <- c(ps_nons, est_ps_rand)

  res <- (Rnons - ps)/sqrt(ps*(1-ps))
  return(res)
}


#' deviance.nonprobsvy
#'
#' deviance.nonprobsvy: Function for deviance residuals computing for propensity score model

deviance.nonprobsvy <- function(X_nons,
                                X_rand,
                                ps_nons,
                                est_ps_rand){

  Rnons <- c(rep(1, nrow(X_nons)), rep(0, nrow(X_rand)))
  s <- c(rep(1, nrow(X_nons)), rep(-1, nrow(X_rand)))
  ps <- c(ps_nons, est_ps_rand)

  res <- s * sqrt(-2 * (Rnons * log(ps) + (1 - Rnons) * log(1 - ps)))

  return(res)
}

cooks.distance.nonprobsvy <- function(){

}


hatvalues.nonprobsvy <- function(){


}


AIC.nonprobsvy <- function(){


}

BIC.nonprobsvy <- function(){



}
