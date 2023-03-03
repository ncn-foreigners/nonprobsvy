cv_nonprobsvy <- function(X_rand, X_nons, R, weights_X, method_selection, K = 10) {

  loc_nons <- which(R == 1)
  loc_rand <- which(R == 0)

  R_nons <- R[loc_nons]
  R_rand <- R[loc_rand]

  weights_nons <- weights_X[loc_nons]
  weights_rand <- weights_X[loc_rand]

  k <- 1

  loss_theta_av <- vector(mode = "numeric", length = 50)
  lambdas <- seq(from = 0.02, to = 4, by = 0.08)

  X_nons <- X_nons[sample(nrow(X_nons)), ]
  X_rand <- X_rand[sample(nrow(X_rand)), ]

  folds_nons <- cut(seq(1,nrow(X_nons)), breaks=K, labels=FALSE) #split nonprobabability sample into K parts
  folds_rand <- cut(seq(1,nrow(X_rand)), breaks=K, labels=FALSE) #split probabability sample into K parts

  # pair K subsets randomly
  sample_nons <- sample(1:K, K, replace = FALSE)
  sample_rand <- sample(1:K, K, replace = FALSE)

  for (lambda in lambdas) {

    loss_theta_vec <- vector(mode = "numeric", length = K)

    for(i in 1:K){

      # train data for X_nons
      idx_nons <- which(folds_nons==sample_nons[i], arr.ind=TRUE)
      X_nons_train <- X_nons[-idx_nons, ]
      R_nons_train <- R_nons[-idx_nons]
      weights_nons_train <- weights_nons[-idx_nons]

      # test data for X_nons
      X_nons_test <- X_nons[idx_nons, ]
      R_nons_test <- R_nons[idx_nons]
      weights_nons_test <- weights_nons[idx_nons]

      # train data for X_rand
      idx_rand <- which(folds_rand==sample_rand[i], arr.ind=TRUE)
      X_rand_train <- X_rand[-idx_rand, ]
      R_rand_train <- R_rand[-idx_rand]
      weights_rand_train <- weights_rand[-idx_rand]

      # test data for X_rand
      X_rand_test <- X_rand[idx_rand, ]
      R_rand_test <- R_rand[idx_rand]
      weights_rand_test <- weights_rand[idx_rand]


      X_train <- rbind(X_rand_train[, -1], X_nons_train[, -1])
      X_test <- rbind(X_rand_test[, -1], X_nons_test[, -1])

      R_train <- c(R_rand_train, R_nons_train)
      R_test <- c(R_rand_test, R_nons_test)

      weights_X_train <- c(weights_rand_train, weights_nons_train)
      weights_X_test <- c(weights_rand_test, weights_nons_test)

      p <- ncol(X_test)

      # initial values for set of parameters
      init_theta <- rep(0, p+1)

      # variables selection using score equation for theta

      par0 <- c(init_theta)
      LAMBDA <- Matrix::Matrix(matrix(0, p+1, p+1), sparse = TRUE)
      it <- 0

      for(jj in 1:100) { # to fix
        it <- it + 1

        u_theta0 <- u_theta(par = par0, R = R_train, X = X_train,
                            weights = weights_X_train,
                            method_selection = method_selection)

        u_theta0_der <- u_theta_der(par = par0, R = R_train, X = X_train,
                                    weights = weights_X_train,
                                    method_selection = method_selection)

        diag(LAMBDA) <- abs(q_lambda(par0, lambda))/(1e-6 + abs(par0))
        par <- par0 + solve(u_theta0_der + LAMBDA, sparse = TRUE) %*% (u_theta0 - LAMBDA %*% par0) # perhaps 'solve' function instead of 'ginv'
        # equation (13) in the article
        if (sum(abs(par - par0)) < 1e-6) break;
        if (sum(abs(par - par0)) > 1000) break;

        par0 <- par

      }

      par <- as.vector(par)
      theta_est <- par

      loss_theta_vec[i] <- loss_theta(par = theta_est, R = R_test, X = X_test,
                                      weights = weights_X_test, method_selection = method_selection)


    }

    loss_theta_av[k] <- mean(loss_theta_vec)

    k <- k + 1

  }

  i <- which.min(loss_theta_av)

  print(lambdas[i])

  lambdas[i]

}
