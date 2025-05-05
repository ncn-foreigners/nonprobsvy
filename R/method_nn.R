#' Mass imputation using nearest neighbours matching method
#'
#' @importFrom stats update
#' @importFrom survey svymean
#' @importFrom stats weighted.mean
#' @importFrom RANN nn2
#'
#'
#' @description Mass imputation using nearest neighbours approach as described in Yang et al. (2021).
#' The implementation is currently based on [RANN::nn2] function and thus it uses
#' Euclidean distance for matching units from \eqn{S_A} (non-probability) to \eqn{S_B} (probability).
#' Estimation of the mean is done using \eqn{S_B} sample.
#'
#'
#' @details Analytical variance
#'
#' The variance of the mean is estimated based on the following approach
#'
#' (a) non-probability part  (\eqn{S_A} with size \eqn{n_A}; denoted as `var_nonprob` in the result)
#'
#' This may be estimated using
#'
#' \deqn{
#' \hat{V}_1 = \frac{1}{N^2}\sum_{i=1}^{S_A}\frac{1-\hat{\pi}_B(\boldsymbol{x}_i)}{\hat{\pi}_B(\boldsymbol{x}_i)}\hat{\sigma}^2(\boldsymbol{x}_i),
#' }
#'
#' where \eqn{\hat{\pi}_B(\boldsymbol{x}_i)} is an estimator of propensity scores which
#' we currently estimate using \eqn{n_A/N} (constant) and \eqn{\hat{\sigma}^2(\boldsymbol{x}_i)} is
#' estimated using based on the average of \eqn{(y_i - y_i^*)^2}.
#'
#' Chlebicki et al. (2025, Algorithm 2) proposed non-parametric mini-bootstrap estimator
#' (without assuming that it is consistent) but with good finite population properties.
#' This bootstrap can be applied using `control_inference(nn_exact_se=TRUE)` and
#' can be summarized as follows:
#'
#' 1. Sample \eqn{n_A} units from \eqn{S_A} with replacement to create \eqn{S_A'} (if pseudo-weights are present inclusion probabilities should be proportional to their inverses).
#' 2. Match units from \eqn{S_B} to \eqn{S_A'} to obtain predictions \eqn{y^*}=\eqn{{k}^{-1}\sum_{k}y_k}.
#' 3. Estimate \eqn{\hat{\mu}=\frac{1}{N} \sum_{i \in S_B} d_i y_i^*}.
#' 4. Repeat steps 1-3 \eqn{M} times (we set \eqn{M=50} in our simulations; this is hard-coded).
#' 5. Estimate \eqn{\hat{V}_1=\text{var}({\hat{\boldsymbol{\mu}}})} obtained from simulations and save it as `var_nonprob`.
#'
#'
#' (b) probability part (\eqn{S_B} with size \eqn{n_B}; denoted as `var_prob` in the result)
#'
#'  This part uses functionalities of the `{survey}` package and the variance is estimated using the following
#'  equation:
#'
#' \deqn{
#' \hat{V}_2=\frac{1}{N^2} \sum_{i=1}^n \sum_{j=1}^n \frac{\pi_{i j}-\pi_i \pi_j}{\pi_{i j}}
#' \frac{y_i^*}{\pi_i} \frac{y_j^*}{\pi_j},
#' }
#'
#'where \eqn{y^*_i} and \eqn{y_j^*} are values imputed imputed as an average
#' of \eqn{k}-nearest neighbour, i.e. \eqn{{k}^{-1}\sum_{k}y_k}. Note that \eqn{\hat{V}_2} in principle can be estimated in various ways depending on the type of the design and whether population size is known or not.
#'
#'
#' @param y_nons target variable from non-probability sample
#' @param X_nons a `model.matrix` with auxiliary variables from non-probability sample
#' @param X_rand a `model.matrix` with auxiliary variables from non-probability sample
#' @param svydesign a svydesign object
#' @param weights case / frequency weights from non-probability sample
#' @param family_outcome a placeholder (not used in `method_nn`)
#' @param start_outcome a placeholder (not used in `method_nn`)
#' @param vars_selection whether variable selection should be conducted
#' @param pop_totals a placeholder (not used in `method_nn`)
#' @param pop_size population size from the `nonprob` function
#' @param control_outcome controls passed by the `control_out` function
#' @param control_inference controls passed by the `control_inf` function
#' @param verbose parameter passed from the main `nonprob` function
#' @param se whether standard errors should be calculated
#'
#' @returns an `nonprob_method` class which is a `list` with the following entries
#'
#' \describe{
#'   \item{model_fitted}{`RANN::nn2` object}
#'   \item{y_nons_pred}{predicted values for the non-probablity sample (query to itself)}
#'   \item{y_rand_pred}{predicted values for the probability sample}
#'   \item{coefficients}{coefficients for the model (if available)}
#'   \item{svydesign}{an updated `surveydesign2` object (new column `y_hat_MI` is added)}
#'   \item{y_mi_hat}{estimated population mean for the target variable}
#'   \item{vars_selection}{whether variable selection was performed (not implemented, for further development)}
#'   \item{var_prob}{variance for the probability sample component (if available)}
#'   \item{var_nonprob}{variance for the non-probability sample component}
#'   \item{var_tot}{total variance, if possible it should be `var_prob+var_nonprob` if not, just a scalar}
#'   \item{model}{model type (character `"nn"`)}
#'   \item{family}{placeholder for the `NN approach` information}
#' }
#'
#' @references
#'  Yang, S., Kim, J. K., & Hwang, Y. (2021). Integration of data from probability surveys and
#'  big found data for finite population inference using mass imputation.
#'  Survey Methodology, June 2021 29 Vol. 47, No. 1, pp. 29-58
#'
#'  Chlebicki, P., Chrostowski, Ł., & Beręsewicz, M. (2025). Data integration of non-probability
#' and probability samples with predictive mean matching. arXiv preprint arXiv:2403.13750.
#'
#' @examples
#'
#' data(admin)
#' data(jvs)
#' jvs_svy <- svydesign(ids = ~ 1,  weights = ~ weight, strata = ~ size + nace + region, data = jvs)
#'
#' res_nn <- method_nn(y_nons = admin$single_shift,
#'                     X_nons = model.matrix(~ region + private + nace + size, admin),
#'                     X_rand = model.matrix(~ region + private + nace + size, jvs),
#'                     svydesign = jvs_svy)
#'
#' res_nn
#'
#' @export
method_nn <- function(y_nons,
                      X_nons,
                      X_rand,
                      svydesign,
                      weights=NULL,
                      family_outcome=NULL,
                      start_outcome=NULL,
                      vars_selection=FALSE,
                      pop_totals=NULL,
                      pop_size=NULL,
                      control_outcome=control_out(),
                      control_inference=control_inf(),
                      verbose=FALSE,
                      se=TRUE) {

  if (missing(y_nons) | missing(X_nons)) {
    stop("`y_nons` and `X_nons`, `X_rand` are required.")
  }

  if (missing(svydesign) | is.null(svydesign)) {
    stop("The NN method is suited only for the unit-level data.")
  }
  if (is.null(weights)) weights <- rep(1, NROW(X_nons))
  if (is.null(pop_size)) pop_size <- sum(weights(svydesign))

  if (verbose) {
    message("Matching units between samples...")
  }

  model_fitted_nons <- RANN::nn2(
    data = X_nons,
    query = X_nons,
    #k = control_outcome$k,
    treetype = control_outcome$treetype,
    searchtype = control_outcome$searchtype
  )

  model_fitted <- RANN::nn2(
    data = X_nons,
    query = X_rand,
    #k = control_outcome$k,
    treetype = control_outcome$treetype,
    searchtype = control_outcome$searchtype
  )

  k_range <- 1:control_outcome$k

  y_rand_pred <- switch(control_outcome$pmm_weights, ## this should be changed to nn_weights
         "none" = apply(model_fitted$nn.idx[, k_range], 1, FUN = function(x) mean(y_nons[x])),
         "dist" = {
           # TODO:: these weights will need to be saved for variance estimation
           if (control_outcome$k == 1) {
             apply(model_fitted$nn.idx[, 1], 1, FUN = function(x) mean(y_nons[x]))
           } else {
             sapply(1:NROW(model_fitted$nn.idx),
                    FUN = function(x) {
                      w_scaled <- max(model_fitted$nn.dists[x, k_range]) - model_fitted$nn.dists[x, k_range]
                      w_scaled <- w_scaled/sum(w_scaled)
                      stats::weighted.mean(y_nons[model_fitted$nn.idx[x, k_range]],
                                    w = w_scaled)
                    })
           }

         }
  )

  y_nons_pred <- apply(model_fitted_nons$nn.idx[, k_range], 1, FUN = function(x) mean(y_nons[x]))

  svydesign_updated <- stats::update(svydesign, y_hat_MI = y_rand_pred)
  svydesign_mean <- survey::svymean( ~ y_hat_MI, svydesign_updated)
  y_mi_hat <- as.numeric(svydesign_mean)

  if (se) {

    var_prob <- as.vector(attr(svydesign_mean, "var"))

    sigma_hat <- mean((y_nons - y_nons_pred)^2)
    est_ps <- NROW(X_nons) / pop_size
    var_nonprob <- 1/ pop_size^2 * (1 - est_ps) / est_ps * sigma_hat
    var_total <- var_prob + var_nonprob

    if (control_inference$nn_exact_se) {

      message("The `nn_exact_se=TRUE` option used. Remember to set the seed for reproducibility.")

      loop_size <- 50

      if (verbose) {
        message("Estimating variance component using mini-bootstrap...")
        pb <- utils::txtProgressBar(min = 0, max = loop_size, style = 3)
      }

      dd <- numeric(loop_size)
      for (jj in 1:loop_size) {

        if (verbose) {
          utils::setTxtProgressBar(pb, jj)
        }

        boot_samp <- sample(1:NROW(X_nons), size = NROW(X_nons), replace = TRUE)
        y_nons_b <- y_nons[boot_samp]
        X_nons_b <- X_nons[boot_samp, , drop = FALSE]

        YY <- RANN::nn2(
          data = X_nons_b,
          query = X_rand,
          #k = control_outcome$k,
          treetype = control_outcome$treetype,
          searchtype = control_outcome$searchtype
        )

        k_range <- 1:control_outcome$k

        y_rand_pred_mini_boot <- switch(control_outcome$pmm_weights, ## this should be changed to nn_weights
                                        "none" = apply(model_fitted$nn.idx[, k_range], 1, FUN = function(x) mean(y_nons_b[x])),
                                        "dist" = {
                                          # TODO:: these weights will need to be saved for variance estimation
                                          if (control_outcome$k == 1) {
                                            apply(model_fitted$nn.idx[, 1], 1, FUN = function(x) mean(y_nons_b[x]))
                                          } else {
                                            sapply(1:NROW(model_fitted$nn.idx),
                                                   FUN = function(x) {
                                                     w_scaled <- max(model_fitted$nn.dists[x, k_range]) - model_fitted$nn.dists[x, k_range]
                                                     w_scaled <- w_scaled/sum(w_scaled)
                                                     stats::weighted.mean(y_nons_b[model_fitted$nn.idx[x, k_range]],
                                                                          w = w_scaled)
                                                   })}
                                          })

        dd[jj] <- stats::weighted.mean(y_rand_pred_mini_boot, weights(svydesign))
      }
      var_nonprob <- var(dd)
      var_total <- var_prob + var_nonprob

    }
  } else {
    var_prob <- var_nonprob <- var_total <- NA
  }

  return(
    structure(
      list(
        model_fitted = model_fitted,
        y_nons_pred = y_nons_pred,
        y_rand_pred = y_rand_pred,
        coefficients = NULL,
        svydesign = if (is.null(svydesign)) svydesign else svydesign_updated,
        y_mi_hat = y_mi_hat,
        vars_selection = vars_selection,
        var_prob = var_prob,
        var_nonprob = var_nonprob,
        var_total = var_total,
        model = "nn",
        family = "NN approach"
      ),
      class = "nonprob_method")
  )
}
