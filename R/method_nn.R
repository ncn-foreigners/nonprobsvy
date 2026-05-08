#' @title Mass Imputation Using Nearest Neighbours Matching Method
#'
#' @importFrom stats update
#' @importFrom survey svytotal
#' @importFrom stats weighted.mean
#' @importFrom RANN nn2
#'
#'
#' @description Mass imputation using nearest neighbours approach as described in Yang et al. (2021).
#' The implementation is currently based on [RANN::nn2] function and thus it uses
#' Euclidean distance for matching units from \eqn{S_A} (non-probability) to \eqn{S_B} (probability).
#' Matching ties are randomized before donor values are aggregated, so tied nearest neighbours
#' are not selected only by input row order. Estimation of the mean is done using \eqn{S_B} sample:
#' when `pop_size` is supplied this is the known-\eqn{N} Horvitz-Thompson mean,
#' otherwise it reduces to the usual ratio mean with \eqn{\hat{N} = \sum_{i\in S_B} d_i}.
#' The `pop_size` argument is not converted into a finite population correction;
#' if an fpc is needed, it should be supplied in `svydesign`, where it is handled
#' by the `{survey}` variance routines.
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
#' estimated using based on the average of \eqn{(y_i - y_i^*)^2}. The \eqn{y_i^*}
#' values used in this proxy are obtained by leave-one-out matching in \eqn{S_A},
#' so a unit is not used as its own donor.
#'
#' Chlebicki et al. (2025, Algorithm 2) proposed non-parametric mini-bootstrap estimator
#' (without assuming that it is consistent) but with good finite population properties.
#' This bootstrap can be applied using `control_inference(nn_exact_se=TRUE)` and
#' can be summarized as follows:
#'
#' 1. Sample \eqn{n_A} units from \eqn{S_A} with replacement to create \eqn{S_A'}.
#'    If non-constant pseudo-weights are supplied through `weights`, sampling probabilities
#'    are proportional to their inverses; equal weights use uniform resampling.
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
#' @param weights case / frequency weights from non-probability sample. If `nn_exact_se=TRUE`,
#'   non-constant weights also define mini-bootstrap sampling probabilities proportional
#'   to their inverses.
#' @param family_outcome a placeholder (not used in `method_nn`)
#' @param start_outcome a placeholder (not used in `method_nn`)
#' @param vars_selection whether variable selection should be conducted
#' @param pop_totals a placeholder (not used in `method_nn`)
#' @param pop_size population size from the `nonprob` function. If `NULL`, the
#'   method uses `sum(weights(svydesign))`. If supplied, it is used as the
#'   known-\eqn{N} denominator for the mean and variance scaling, but it does not
#'   modify the finite population correction of `svydesign`.
#' @param control_outcome controls passed by the `control_out` function
#' @param control_inference controls passed by the `control_inf` function
#' @param verbose parameter passed from the main `nonprob` function
#' @param se whether standard errors should be calculated
#' @param nn_matches optional precomputed nearest-neighbour search results for
#'   internal reuse across outcomes. If supplied, it should be a list with
#'   `rand` and `nons` entries from [RANN::nn2()].
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
                      se=TRUE,
                      nn_matches=NULL) {

  if (missing(y_nons) | missing(X_nons)) {
    stop("`y_nons` and `X_nons`, `X_rand` are required.")
  }

  if (missing(svydesign) | is.null(svydesign)) {
    stop("The NN method is suited only for the unit-level data.")
  }
  X_nons <- as.matrix(X_nons)
  X_rand <- as.matrix(X_rand)
  if (is.null(weights)) weights <- rep(1, NROW(X_nons))
  if (is.null(pop_size)) pop_size <- sum(weights(svydesign))
  boot_prob <- if (length(unique(weights)) == 1) NULL else 1 / weights

  randomize_returned_ties <- function(idx_mat, dist_mat) {
    for (i in seq_len(NROW(idx_mat))) {
      idx <- idx_mat[i, ]
      dists <- dist_mat[i, ]
      keep <- !is.na(idx)
      if (!any(keep)) next
      if (!any(duplicated(dists[keep]))) next

      ord <- order(dists[keep], stats::runif(sum(keep)))
      idx_mat[i, keep] <- idx[keep][ord]
      dist_mat[i, keep] <- dists[keep][ord]
    }

    list(idx = idx_mat, dists = dist_mat)
  }

  exact_matches <- function(data, query, k, exclude_self = FALSE, response = NULL) {
    data <- as.matrix(data)
    query <- as.matrix(query)
    n_data <- NROW(data)
    k <- min(k, n_data - as.integer(exclude_self))
    if (k <= 0) {
      return(list(
        idx = matrix(NA_integer_, nrow = NROW(query), ncol = 1),
        dists = matrix(NA_real_, nrow = NROW(query), ncol = 1)
      ))
    }

    idx_mat <- matrix(NA_integer_, nrow = NROW(query), ncol = k)
    dist_mat <- matrix(NA_real_, nrow = NROW(query), ncol = k)
    tol <- sqrt(.Machine$double.eps)

    for (i in seq_len(NROW(query))) {
      dists_all <- sqrt(rowSums(sweep(data, 2, query[i, ], FUN = "-")^2))
      if (exclude_self && i <= length(dists_all)) dists_all[i] <- Inf
      remaining <- which(is.finite(dists_all))
      selected <- integer()
      selected_dists <- numeric()

      while (length(selected) < k && length(remaining)) {
        d_min <- min(dists_all[remaining])
        candidates <- remaining[abs(dists_all[remaining] - d_min) <= tol * max(1, d_min)]
        need <- k - length(selected)
        if (length(candidates) > need &&
            (is.null(response) || length(unique(response[candidates])) > 1)) {
          take <- sample(candidates, need)
        } else {
          take <- utils::head(candidates, need)
        }
        selected <- c(selected, take)
        selected_dists <- c(selected_dists, dists_all[take])

        if (length(selected) >= k) break
        remaining <- setdiff(remaining, candidates)
      }

      if (length(selected)) {
        idx_mat[i, seq_along(selected)] <- selected
        dist_mat[i, seq_along(selected)] <- selected_dists
      }
    }

    list(idx = idx_mat, dists = dist_mat)
  }

  use_exact_matches <- function(data, query) {
    n_query <- NROW(query)
    as.numeric(NROW(data)) * as.numeric(n_query) <= 2e6
  }

  predict_from_matches <- function(nn_fit,
                                   response,
                                   weighting = "none",
                                   exclude_self = FALSE,
                                   data = NULL,
                                   query = NULL,
                                   k = NULL) {
    if (!is.null(data) && !is.null(query) && !is.null(k) &&
        use_exact_matches(data, query)) {
      matches <- exact_matches(data = data, query = query, k = k,
                               exclude_self = exclude_self, response = response)
      idx_mat <- matches$idx
      dist_mat <- matches$dists
      exclude_self <- FALSE
    } else {
      idx_mat <- nn_fit$nn.idx[, seq_len(NCOL(nn_fit$nn.idx)), drop = FALSE]
      dist_mat <- nn_fit$nn.dists[, seq_len(NCOL(nn_fit$nn.dists)), drop = FALSE]
      randomized <- randomize_returned_ties(idx_mat, dist_mat)
      idx_mat <- randomized$idx
      dist_mat <- randomized$dists
    }

    vapply(seq_len(NROW(idx_mat)), FUN.VALUE = numeric(1), FUN = function(i) {
      idx <- idx_mat[i, ]
      dists <- dist_mat[i, ]
      keep <- !is.na(idx)
      idx <- idx[keep]
      dists <- dists[keep]

      if (exclude_self) {
        keep <- idx != i
        idx <- idx[keep]
        dists <- dists[keep]
      }

      if (!length(idx)) {
        return(response[i])
      }

      if (weighting == "dist" && length(idx) > 1) {
        w_scaled <- max(dists) - dists
        if (sum(w_scaled) == 0) {
          return(mean(response[idx]))
        }

        return(stats::weighted.mean(response[idx], w = w_scaled / sum(w_scaled)))
      }

      mean(response[idx])
    })
  }

  if (verbose) {
    message("Matching units between samples...")
  }

  model_fitted_nons <- if (!is.null(nn_matches$nons)) {
    nn_matches$nons
  } else {
    RANN::nn2(
      data = X_nons,
      query = X_nons,
      k = min(control_outcome$k + 1L, NROW(X_nons)),
      treetype = control_outcome$treetype,
      searchtype = control_outcome$searchtype
    )
  }

  model_fitted <- if (!is.null(nn_matches$rand)) {
    nn_matches$rand
  } else {
    RANN::nn2(
      data = X_nons,
      query = X_rand,
      k = control_outcome$k,
      treetype = control_outcome$treetype,
      searchtype = control_outcome$searchtype
    )
  }

  y_rand_pred <- switch(control_outcome$pmm_weights, ## this should be changed to nn_weights
                        "none" = predict_from_matches(model_fitted, y_nons,
                                                       data = X_nons, query = X_rand,
                                                       k = control_outcome$k),
                        "dist" = predict_from_matches(model_fitted, y_nons, weighting = "dist",
                                                       data = X_nons, query = X_rand,
                                                       k = control_outcome$k)
  )

  svydesign_updated <- stats::update(svydesign, y_hat_MI = y_rand_pred)
  # svytotal() honors any fpc already stored in svydesign; pop_size is only the mean denominator.
  svydesign_total <- survey::svytotal( ~ y_hat_MI, svydesign_updated)
  y_mi_hat <- as.numeric(svydesign_total) / pop_size
  y_nons_pred <- NULL

  if (se) {

    var_prob <- as.vector(attr(svydesign_total, "var")) / pop_size^2
    var_nonprob <- 0
    var_total <- var_prob

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

        boot_samp <- sample(1:NROW(X_nons), size = NROW(X_nons), replace = TRUE, prob = boot_prob)
        y_nons_b <- y_nons[boot_samp]
        X_nons_b <- X_nons[boot_samp, , drop = FALSE]

        boot_matches <- RANN::nn2(
          data = X_nons_b,
          query = X_rand,
          k = control_outcome$k,
          treetype = control_outcome$treetype,
          searchtype = control_outcome$searchtype
        )

        y_rand_pred_mini_boot <- switch(control_outcome$pmm_weights, ## this should be changed to nn_weights
                                        "none" = predict_from_matches(boot_matches, y_nons_b,
                                                                      data = X_nons_b, query = X_rand,
                                                                      k = control_outcome$k),
                                        "dist" = predict_from_matches(boot_matches, y_nons_b, weighting = "dist",
                                                                      data = X_nons_b, query = X_rand,
                                                                      k = control_outcome$k))

        dd[jj] <- sum(weights(svydesign) * y_rand_pred_mini_boot) / pop_size
      }
      var_nonprob <- var(dd)
      var_total <- var_prob + var_nonprob

    } else {
      # Use leave-one-out matching for the non-probability variance proxy.
      y_nons_pred <- predict_from_matches(model_fitted_nons, y_nons, exclude_self = TRUE,
                                          data = X_nons, query = X_nons,
                                          k = control_outcome$k)

      sigma_hat <- mean((y_nons - y_nons_pred)^2)
      est_ps <- NROW(X_nons) / pop_size
      var_nonprob <- 1/ pop_size^2 * (1 - est_ps) / est_ps * sigma_hat
      var_total <- var_prob + var_nonprob
    }
  } else {
    var_prob <- var_nonprob <- var_total <- NA
  }

  if (is.null(y_nons_pred)) {
    y_nons_pred <- predict_from_matches(model_fitted_nons, y_nons, exclude_self = TRUE,
                                        data = X_nons, query = X_nons,
                                        k = control_outcome$k)
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
