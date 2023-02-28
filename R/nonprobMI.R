#' nonprobMI
#
#' nonprobMI: Function for inference based on nonprobability surveys using mass imputation
#
#' @param outcome - `formula`, the outcome equation.
#' @param data - an optional `data.frame` with data from the nonprobability sample.
#' @param svydesign - an optional `svydesign` object (from the survey package) containing probability sample.
#' @param family_outcome - a `character` string describing the error distribution and link function to be used in the model. Default is "gaussian". Currently supports: gaussian with identity link, poisson and binomial.
#' @param method_outcome - a `character` with method for response variable estimation
#' @param subset - an optional `vector` specifying a subset of observations to be used in the fitting process.
#' @param strata - an optional `vector` specifying strata.
#' @param weights - an optional `vector` of ‘prior weights’ to be used in the fitting process. Should be NULL or a numeric vector. It is assumed that this vector contains frequency or analytic weights
#' @param na_action a function which indicates what should happen when the data contain `NAs`.
#' @param control_outcome a list indicating parameters to use in fitting model for outcome variable
#' @param control_inference a list indicating parameters to use in inference based on probablity and nonprobability samples, contains parameters such as estimation method or variance method
#' @param start - an optional `list` with starting values for the parameters of the selection and outcome equation
#' @param verbose - verbose, numeric
#' @param contrasts a
#' @param model a
#' @param x a
#' @param y a
#' @param ... a
#'
#' @importFrom stats glm.fit
#' @importFrom stats model.frame
#' @importFrom stats model.matrix
#' @importFrom stats update
#' @importFrom stats qnorm
#' @export

nonprobMI <- function(outcome,
                      data,
                      svydesign,
                      method_outcome,
                      family_outcome = "gaussian",
                      subset,
                      strata,
                      weights,
                      na_action,
                      control_outcome,
                      control_inference = controlInf(var_method = "analytic"),
                      start,
                      verbose,
                      contrasts,
                      model,
                      x,
                      y,
                      ...) {

  weights <- rep.int(1, nrow(data)) # to remove

  XY_nons <- model.frame(outcome, data)
  y_name <- colnames(XY_nons)[1]
  X_nons <- model.matrix(XY_nons, data) # matrix of the  nonprobability sample
  nons_names <- colnames(X_nons[,-1])
  if (all(nons_names %in% colnames(svydesign$variables))) {

    X_rand <- as.matrix(cbind(1, svydesign$variables[,nons_names])) #matrix of probability sample with intercept

  } else {

    stop("variable names in data and svydesign do not match")

  }
  y_nons <- XY_nons[,1]

  R_nons <- rep(1, nrow(X_nons))
  R_rand <- rep(0, nrow(X_rand))
  R <- c(R_nons, R_rand)

  loc_nons <- which(R == 1)
  loc_rand <- which(R == 0)

  n_nons <- nrow(X_nons)
  n_rand <- nrow(X_rand)
  X <- rbind(X_nons, X_rand)

  ps_rand <- svydesign$prob
  weights_rand <- 1/ps_rand
  N_est_rand <- sum(weights_rand)


  ## estimation

  if (control_outcome$method == "glm") {

    model_nons <- nonprobMI_fit(x = X_nons,
                                y = y_nons,
                                weights = weights,
                                family_outcome = family_outcome)

    model_nons_coefs <- as.matrix(model_nons$coefficients)


    y_rand_pred <-  as.numeric(X_rand %*% model_nons_coefs) # y_hat for probability sample

    y_nons_pred <- as.numeric(X_nons %*% model_nons_coefs)

  } else if (control_outcome$method == "nn") {

    model_rand <- nonprobMI_nn(data = X_nons, query = X_rand,
                               k = 5, treetype = "kd", searchtype = "standard")
    model_nons <- nonprobMI_nn(data = X_nons, query = X_nons,
                               k = 5, treetype = "kd", searchtype = "standard")

    y_rand_pred <- vector(mode = "numeric", length = n_rand)
    y_nons_pred <- vector(mode = "numeric", length = n_nons)


    for (i in 1:n_rand) {

      idx <- model_rand$nn.idx[i,]
      y_rand_pred[i] <- mean(y_nons[idx])

    }

    for (i in 1:n_nons) {

      idx <- model_nons$nn.idx[i,]
      y_nons_pred[i] <- mean(y_nons[idx])
    }


  }

  # updating probability sample by adding y_hat variable
  svydesign <- stats::update(svydesign,
                             y_hat_MI = y_rand_pred)

  nonprobMI_inference <- function(...) {

    mu_hat <- mu_hatMI(y = y_rand_pred,
                       weights = weights_rand,
                       N = N_est_rand)

    # design based variance estimation based on approximations of the second-order inclusion probabilities

    if (control_inference$var_method == "analytic") {

      svydesign_mean <- survey::svymean(~y_hat_MI, svydesign) # SE for y_hat_MI is equal to prob SE
      var_prob <- as.vector(attr(svydesign_mean, "var")) #probability component

      if (control_outcome$method == "nn") {

        method_selection <- "logit"
        method <- method_selection

        if (is.character(method)) {
          method <- get(method, mode = "function", envir = parent.frame())
        }
        if (is.function(method)) {
          method <- method()
        }

        ps_method <- method$make_propen_score # function for propensity score estimation
        loglike <- method$make_log_like
        gradient <- method$make_gradient
        hessian <- method$make_hessian

        optim_method <- "NR"

        # initial values for propensity score estimation
        start <- start_fit(X,
                           R,
                           weights,
                           weights_rand,
                           method_selection)

        log_like <- loglike(X_nons, X_rand, weights_rand)
        gradient <- gradient(X_nons, X_rand, weights_rand)
        hessian <- hessian(X_nons, X_rand, weights_rand)

        maxLik_nons_obj <- ps_method(X_nons, log_like, gradient, hessian, start, optim_method)

        ps_nons <- maxLik_nons_obj$ps


        sigma_hat <- switch(family_outcome,
                            "gaussian" = mean((y_nons - y_nons_pred)^2),
                            "binomial" = NULL,
                            "poisson" = NULL)

        N_est_nons <- sum(1/ps_nons)

        var_nonprob <- 1/N_est_nons^2 * sum((1 - ps_nons)/ps_nons * sigma_hat)


      } else if (control_outcome$method == "glm") {


        mx <- 1/N_est_rand * colSums(weights_rand * X_rand)
        mh <- 0
        for (i in 1:n_nons) { # matrix product instead of a loop in a near future

          xx <- t(X_nons[i,]) %*% X_nons[i,]
          mh <- mh + xx
        }

        c <- 1/(1/n_nons * mh) %*% mx
        e <- y_nons - y_nons_pred

        # nonprobability component
        var_nonprob <- 1/n_nons^2 * t(as.matrix(e^2)) %*% (X_nons %*% t(c))^2 #ok, but may be too large
        #var_nonprob <- 1/n_nons^2 * (t(as.matrix(e)) %*% (X_nons %*% t(as.matrix(c))))^2


        var_nonprob <- as.vector(var_nonprob)

      }

      se_nonprob <- sqrt(var_nonprob)
      se_prob <- sqrt(var_prob)

      # variance
      var <- var_nonprob + var_prob

    } else if (control_inference$var_method == "bootstrap") {


      # bootstrap variance
      var <- bootMI(X_rand,
                    X_nons,
                    weights,
                    y_nons,
                    family_outcome,
                    1000,
                    weights_rand,
                    mu_hat,
                    svydesign)

      inf <- "not computed for bootstrap variance"

    }


    se <- sqrt(var)

    alpha <- control_inference$alpha
    z <- stats::qnorm(1-alpha/2)

    # confidence interval based on the normal approximation
    ci <- c(mu_hat - z*se, mu_hat + z*se)



    structure(
      list(mean = mu_hat,
           #VAR = var,
           SE = se,
           #variance_nonprob = ifelse(control_inference$var_method == "analytic", v_a, inf),
           #variance_prob = ifelse(control_inference$var_method == "analytic", v_b, inf),
           SE_nonprob = ifelse(control_inference$var_method == "analytic", se_nonprob, inf),
           SE_prob = ifelse(control_inference$var_method == "analytic", se_prob, inf),
           CI = ci
           #beta = model_nons_coefs
           ),
      class = "Mass imputation")

  }

  # inference based on mi method
  infer_nons <- nonprobMI_inference()

  infer_nons
}


#' nonprobMI_fit
#
#' nonprobMI_fit: Function for outcome variable estimation based on nonprobability sample and using model based approach
#'
#' @param outcome - `formula`, the outcome equation.
#' @param data - an optional `data.frame` with data from the nonprobability sample.
#' @param svydesign - an optional `svydesign` object (from the survey package) containing probability sample.
#' @param family_outcome - a `character` string describing the error distribution and link function to be used in the model. Default is "gaussian". Currently supports: gaussian with identity link, poisson and binomial.
#' @param control_outcome - a
#' @param start - a
#' @param weights - an optional `vector` of ‘prior weights’ to be used in the fitting process. Should be NULL or a numeric vector. It is assumed that this vector contains frequency or analytic weights
#' @param verbose - a
#' @param model - a
#' @param x - a
#' @param y - a
#'


nonprobMI_fit <- function(outcome,
                          data,
                          weights,
                          svydesign,
                          family_outcome,
                          control_outcome = controlOut(),
                          start = NULL,
                          verbose,
                          model,
                          x,
                          y) {


  family <- family_outcome

  if (is.character(family)) {
    family <- get(family, mode = "function", envir = parent.frame())
  }
  if (is.function(family)) {
    family <- family()
  }

  model_nons <- stats::glm.fit(x = x,
                               y = y,
                               weights = weights,
                               start = start,
                               control = list(control_outcome$epsilon,
                                              control_outcome$maxit,
                                              control_outcome$trace),
                               family = family)


  model_nons

}

#' nonprobMI_nn
#
#' nonprobMI_nn: Function for outcome variable estimation based on nonprobability sample and using predictive mean matching
#'
#' @param data - an optional `data.frame` with data from the nonprobability sample.
#' @param query - a
#' @param k - a
#' @param treetype - a
#' @param searchtype - a
#' @param radius - a
#' @param eps - a


nonprobMI_nn <- function(data,
                         query,
                         k,
                         treetype,
                         searchtype,
                         radius = 0,
                         eps = 0) {


  model_nn <- nn2(data = data,
                  query = query,
                  k = k,
                  treetype = treetype,
                  searchtype = searchtype,
                  radius = radius,
                  eps = eps)

  model_nn

}

#' mu_hatMI
#
#' mu_hatMI: Function for outcome variable estimation based on mass imputation
#' @param y - a
#' @param weights - a
#' @param N - a

mu_hatMI <- function(y, weights, N) {

  mu_hat <- 1/N * sum(weights * y)

  mu_hat


}



