#' Internal functions

bootMI <- function(X_rand,
                   X_nons,
                   weights,
                   y,
                   family_outcome,
                   num_boot,
                   d,
                   mu_hat,
                   ...
                   ){


  mu_hats <- vector(mode = "numeric", length = num_boot)
  n_nons <- nrow(X_nons)
  n_rand <- nrow(X_rand)

  k <- 1

  while (k <= num_boot) {

    strap_rand <- sample.int(replace = FALSE, n = n_rand, prob = 1/d) # to change -> replicate weights based on sampling design
    strap <- sample.int(replace = TRUE, n = n_nons)
    weights_strap <- weights[strap]
    Xnons_strap <- X_nons[strap,]
    y_strap <- y[strap]
    d_strap <- d[strap_rand]
    Nstrap <- sum(d_strap)


    model_strap <- nonprobMI_fit(x = Xnons_strap,
                                 y = y_strap,
                                 weights = weights_strap,
                                 family_outcome = family_outcome)

    beta <- model_strap$coefficients

    ystrap_nons <- as.numeric(Xnons_strap %*% beta)

    ystrap_rand <- as.numeric(X_rand %*% beta)

    mu_hat_boot <- mu_hatMI(ystrap_rand, d_strap, Nstrap)
    mu_hats[k] <- mu_hat_boot

    k <- k + 1

  }

  boot_var <- 1/num_boot * sum((mu_hats - mean(mu_hats))^2) #mean may be replaced by mu_hatMI

  return(boot_var)

}

bootIPW <- function(X_rand,
                    X_nons,
                    weights,
                    y,
                    family_outcome,
                    num_boot,
                    d,
                    mu_hat,
                    dependency,
                    ...){





}

bootDR <- function(X_rand,
                   X_nons,
                   weights,
                   y,
                   family_outcome,
                   num_boot,
                   d,
                   mu_hat,
                   dependency,
                   ...){



}
