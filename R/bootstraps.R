#' Internal functions
#' @importFrom survey as.svrepdesign

bootMI <- function(X_rand,
                   X_nons,
                   weights,
                   y,
                   family_outcome,
                   num_boot,
                   weights_rand,
                   mu_hat,
                   svydesign,
                   rep_type = NULL,
                   ...
                   ){


  mu_hats <- vector(mode = "numeric", length = num_boot)
  n_nons <- nrow(X_nons)
  n_rand <- nrow(X_rand)
  N <- sum(weights_rand)
  rep_weights <- survey::as.svrepdesign(svydesign, type = rep_type, replicates = num_boot)$repweights$weights

  k <- 1

  while (k <= num_boot) { # to analyse


    strap <- sample.int(replace = TRUE, n = n_nons)
    weights_strap <- weights[strap]
    X_nons_strap <- X_nons[strap,]
    y_strap <- y[strap]

    #using svy package
    strap_rand_svy <- which(rep_weights[,k] != 0)
    weights_rand_strap_svy <- rep_weights[,k] * weights_rand


    model_strap <- nonprobMI_fit(x = X_nons_strap,
                                 y = y_strap,
                                 weights = weights_strap,
                                 family_outcome = family_outcome)

    beta <- model_strap$coefficients

    ystrap_rand <- as.numeric(X_rand %*% beta)

    mu_hat_boot <- mu_hatMI(ystrap_rand, weights_rand_strap_svy, N)
    mu_hats[k] <- mu_hat_boot

    k <- k + 1

  }

  boot_var <- 1/num_boot * sum((mu_hats - mu_hat)^2) # mean may be replaced by mu_hatMI

  boot_var

}

bootIPW <- function(X_rand,
                    X_nons,
                    weights,
                    y,
                    family_outcome,
                    num_boot,
                    weights_rand,
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
                   weights_rand,
                   mu_hat,
                   dependency,
                   ...){



}
