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

  R_nons <- rep(1, nrow(X_nons))
  R_rand <- rep(0, nrow(X_rand))
  R <- c(R_nons, R_rand)
  loc_nons <- which(R == 1)
  loc_rand <- which(R == 0)
  X <- rbind(X_nons, X_rand)
  XY <- cbind(X, R)

  O <- O_hat_model(x = X,
                   X = X,
                   R = R)

  if(dependent){

    L <- L_hat_model(x = X_rand,
                     X = X,
                     r = R_rand)

    ps <- O * L/(d - 1)
    weights <- 1/ps[loc_nons]

  } else {
    ps <- O/(O + d - 1)
    weights <- 1/ps[loc_nons]
  }

  return(weights)
}


bootstrap <- function(X_rand,
                      X_nons,
                      y,
                      weights, #for probability and nonprobability sample
                      O_hat,
                      L_hat,
                      d_rand,
                      d,
                      num_boot,
                      dependency,
                      N){


  n_nons <- nrow(X_nons)
  n_rand <- nrow(X_rand)
  mu_hats <- vector(mode = "numeric", length = num_boot)
  R_nons <- rep(1, n_nons)
  R_rand <- rep(0, n_rand)
  R <- c(R_nons, R_rand)
  X <- rbind(X_nons, X_rand)
  R_r <- 1-R


  weights <- weights*N/sum(weights)
  prob_floor <- weights - floor(weights)
  weights_strap <- weights
  weights_idx <- sample.int(nrow(weights), prob = prob_floor)
  weights_strap[weights_idx] <- floor(weights_strap[weights_idx])
  weights_strap[!weights_idx] <- ceiling(weights_strap[weights_idx])

  pseudo_pop <- rep(X[1:nrow(X), ], weights) # pseudo statistics - to consider
  pseudo_R <- rep(R_r, weights)
  pseudo_L <- rep(L_hat, weights)
  pseudo_O <- rep(O_hat, weights)
  pseudo_d <- rep(d, weights)

  k <- 1

  if(!dependency){


    while(k <= num_boot){

      strap_rand <- sample(x = pseudo_pop, n = n_rand, replace = TRUE, prob = 1/d_rand)
      strap_nons <- sample(x = pseudo_pop, n = n_nons, replace = TRUE, prob = 1/weights)

      X_rand_strap <- pseudo_pop[strap_rand]
      X_nons_strap <- pseudo_pop[strap_nons]
      X <- rbind(X_nons_strap, X_rand_strap)

      #Removing overlapping units - to do

      O_strap <- O_hat_model(x = X,
                             X = X,
                             R = R)

      ps <- O/(O + d - 1)
      weights <- 1/ps[loc_nons]
      N <- sum(weights)

      mu_hat_boot <- mu_hatIPW(y = y,
                              weights = weights,
                              N = N)

      mu_hats[k] <- mu_hat_boot

      k <- k + 1

    }

  } else {

    #L <- L_hat_model(x = pseudo_pop,
                 #    X = X,
                  #   R = pseudo_R)

    ps1 <- vector(mode = "numeric", length = nrow(pseudo_R))

    for (i in 1:nrow(pseudo_R)){

      ifelse(pseudo_R[i] == 1, ps1[i] <- 1 - pseudo_L[i], ps1[i] <- 1 - pseudo_d[i] + d[i]*pseudo_O[i]*pseudo_L[i] + pseudo_L[i]*(pseudo_d[i] - 1))

    }


    #weights for strap_nons to compute

    while (k <= num_boot){

      strap_rand <- sample(x = pseudo_pop, n = n_rand, replace = TRUE, prob = 1/d_rand)
      strap_nons <- sample(x = pseudo_pop, n = n_nons, replace = TRUE, prob = ps1)

      X_rand_strap <- pseudo_pop[strap_rand]
      X_nons_strap <- pseudo_pop[strap_nons]
      X <- rbind(X_nons_strap, X_rand_strap)


      O_strap <- O_hat_model(x = X,
                             X = X,
                             R = R)


      L <- L_hat_model(x = X_rand_strap,
                       X = X,
                       R = R_rand)

      ps <- O * L/(d - 1)
      weights <- 1/ps[loc_nons]

      mu_hat_boot <- mu_hatIPW(y = y,
                               weights = weights,
                               N = N)

      mu_hats[k] <- mu_hat_boot

      k <- k + 1

    }

  }

  var <- 1/(num_boot-1)*sum((mu_hats - mean(mu_hats))^2)

  se <- sqrt(var)

}





L_hat_model <- function(x,
                        X,
                        R){

  modelL <- stats::glm.fit(x = x,
                           y = R,
                           family = binomial(link = "logit")) # glm.fit: algorithm did not converge

  L <- 1/(1 + exp(X %*% modelL$coefficients))

  L

}



O_hat_model <- function(x,
                        X,
                        R){

  modelO <- stats::glm.fit(x = x,
                           y = R,
                           family = binomial(link = "logit"))


  O <- exp(X %*% modelO$coefficients)

  O

}


