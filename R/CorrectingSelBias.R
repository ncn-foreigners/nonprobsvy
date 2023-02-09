#' overlap
#'
#' overlap: Function for correcting selection bias considering the possible overlapping and dependency between the nonprobability and probability sample
#'
#' @importFrom stats glm.fit
#' @importFrom betareg betareg.fit



overlap <- function(X_nons,
                    X_rand,
                    d_rand,
                    dependent,
                    method.selection){

  betamodel <- betareg.fit(x = X_rand,
                           y = 1/d_rand,
                           link = method.selection)

  method <- method.selection
  if (is.character(method)) {
    method <- get(method, mode = "function", envir = parent.frame())
  }
  if (is.function(method)) {
    method <- method()
  }

  inv_link <- method$linkInv
  prob <- inv_link(X_nons %*% betamodel$coefficients$mean)
  d_rnons <- 1/prob

  d <- c(d_rnons, d_rand) #drnons is a propensity score for the probability sample for units in the nonprobability sample

  Rnons <- rep(1, nrow(X_nons))
  Rrand <- rep(0, nrow(X_rand))
  R <- c(Rnons, Rrand)
  loc.nons <- which(R == 1)
  loc.rand <- which(R == 0)
  X <- rbind(X_nons, X_rand)
  XY <- cbind(X, R)

  modelO <- stats::glm.fit(x = X,
                           y = R,
                           family = binomial(link = "logit"))


  O <- exp(X %*% modelO$coefficients)

  if(dependent){

    modelL <- stats::glm.fit(x = X_rand,
                             y = Rrand,
                             family = binomial(link = "logit")) # glm.fit: algorithm did not converge

    L <- 1/(1 + exp(X %*% modelL$coefficients))

    ps <- O * L/(d - 1)
    weights <- 1/ps[loc.nons]

  } else {
    ps <- O/(O + d - 1)
    weights <- 1/ps[loc.nons]
  }

  return(weights)
}



