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
#' @importFrom RANN nn2
#' @importFrom stats terms
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

  # model for outcome formula
  OutcomeModel <- model_frame(formula = outcome, data = data, svydesign = svydesign)
  X_nons <- OutcomeModel$X_nons
  X_rand <- OutcomeModel$X_rand
  nons_names <- OutcomeModel$nons_names
  y_nons <- OutcomeModel$y_nons

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

    # Estimation for outcome model
    model_out <- internal_outcome(X_nons,
                                  X_rand,
                                  y_nons,
                                  weights,
                                  family_outcome)

    y_rand_pred <- model_out$y_rand_pred
    y_nons_pred <- model_out$y_nons_pred
    model_nons_coefs <- model_out$model_nons_coefs

  } else if (control_outcome$method == "nn") {

    model_rand <- nonprobMI_nn(data = X_nons, query = X_rand,
                               k = control_outcome$k, treetype = "kd", searchtype = "standard")
    model_nons <- nonprobMI_nn(data = X_nons, query = X_nons,
                               k = control_outcome$k, treetype = "kd", searchtype = "standard")
    y_rand_pred <- vector(mode = "numeric", length = n_rand)
    y_nons_pred <- vector(mode = "numeric", length = n_nons)
    model_nons_coefs <- "Non-parametric method for outcome model"


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

  mu_hat <- mu_hatMI(y = y_rand_pred,
                     weights = weights_rand,
                     N = N_est_rand)

  # design based variance estimation based on approximations of the second-order inclusion probabilities

  if (control_inference$var_method == "analytic") {

    svydesign_mean <- survey::svymean(~y_hat_MI, svydesign)
    var_prob <- as.vector(attr(svydesign_mean, "var")) #probability component

    if (control_outcome$method == "nn") {

      method_selection <- "logit"
      optim_method = "NR"

      # Estimation for selection model
      model_sel <- internal_selection(X,
                                      X_nons,
                                      X_rand,
                                      weights,
                                      weights_rand,
                                      R,
                                      method_selection,
                                      optim_method)

      maxLik_nons_obj <- model_sel$maxLik_nons_obj
      maxLik_rand_obj <- model_sel$maxLik_rand_obj
      ps_nons <- maxLik_nons_obj$ps


      sigma_hat <- switch(family_outcome,
                          "gaussian" = mean((y_nons - y_nons_pred)^2),
                          "binomial" = y_nons_pred*(1 - y_nons_pred),
                          "poisson" = mean(y_nons_pred))
      N_est_nons <- sum(1/ps_nons)
      var_nonprob <- n_rand/N_est_nons^2 * sum((1 - ps_nons)/ps_nons * sigma_hat)

    } else if (control_outcome$method == "glm") {


      mx <- 1/N_est_rand * colSums(weights_rand * X_rand)
      c <- solve(1/n_nons * t(X_nons) %*% X_nons) %*% mx
      e <- y_nons - y_nons_pred

      # nonprobability component
      var_nonprob <- 1/n_nons^2 * t(as.matrix(e^2)) %*% (X_nons %*% c)^2
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
                  svydesign,
                  rep_type = control_inference$rep_type)
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
         CI = ci,
         beta = model_nons_coefs
         ),
    class = "Mass imputation")

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


  model_nn <- RANN::nn2(data = data,
                        query = query,
                        k = k,
                        treetype = treetype,
                        searchtype = searchtype,
                        radius = radius,
                        eps = eps)

  model_nn

}



