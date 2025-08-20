#' @title Mass Imputation Using Predictive Mean Matching Method
#'
#' @description
#' Model for the outcome for the mass imputation estimator. The implementation is currently based on [RANN::nn2] function and thus it uses Euclidean distance for matching units from \eqn{S_A} (non-probability) to \eqn{S_B} (probability) based on predicted values from model \eqn{\boldsymbol{x}_i} based
#' either on `method_glm` or `method_npar`. Estimation of the mean is done using \eqn{S_B} sample.
#'
#' This implementation extends Yang et al. (2021) approach as described in Chlebicki et al. (2025), namely:
#'
#' \describe{
#'  \item{pmm_weights}{if k>1 weighted aggregation of the mean for a given unit is used. We use distance
#'  matrix returned by [RANN::nn2] function (`pmm_weights` from the [control_out()] function)}
#'  \item{nn_exact_se}{if the non-probability sample is small we recommend using a mini-bootstrap
#'  approach to estimate variance from the non-probability sample  (`nn_exact_se` from the [control_inf()] function)}
#'  \item{pmm_k_choice}{the main `nonprob` function allows for dynamic selection of `k` neighbours based on the
#'  variance minimization procedure (`pmm_k_choice` from the [control_out()] function)}
#' }
#'
#' @details
#'
#' Matching
#'
#' In the package we support two types of matching:
#'
#' 1. \eqn{\hat{y} - \hat{y}} matching (default; `control_out(pmm_match_type = 1)`).
#' 2. \eqn{\hat{y} - y} matching (`control_out(pmm_match_type = 2)`).
#'
#' Analytical variance
#'
#' The variance of the mean is estimated based on the following approach

#' (a) non-probability part  (\eqn{S_A} with size \eqn{n_A}; denoted as `var_nonprob` in the result) is currently estimated using the non-parametric mini-bootstrap estimator proposed by
#' Chlebicki et al. (2025, Algorithm 2). It is not proved to be consistent but with good finite population properties.
#' This bootstrap can be applied using `control_inference(nn_exact_se=TRUE)` and
#' can be summarized as follows:
#'
#' 1. Sample \eqn{n_A} units from \eqn{S_A} with replacement to create \eqn{S_A'} (if pseudo-weights are present inclusion probabilities should be proportional to their inverses).
#' 2. Estimate regression model \eqn{\mathbb{E}[Y|\boldsymbol{X}]=m(\boldsymbol{X}, \cdot)} based on \eqn{S_{A}'} from step 1.
#' 3. Compute \eqn{\hat{\nu}'(i,t)} for \eqn{t=1,\dots,k, i\in S_{B}} using estimated \eqn{m(\boldsymbol{x}', \cdot)} and \eqn{\left\lbrace(y_{j},\boldsymbol{x}_{j})| j\in S_{A}'\right\rbrace}.
#' 4. Compute \eqn{\displaystyle\frac{1}{k}\sum_{t=1}^{k}y_{\hat{\nu}'(i)}} using \eqn{Y} values from \eqn{S_{A}'}.
#' 5. Repeat steps 1-4 \eqn{M} times (we set (hard-coded) \eqn{M=50} in our code).
#' 6. Estimate \eqn{\hat{V}_1=\text{var}({\hat{\boldsymbol{\mu}}})} obtained from simulations and save it as `var_nonprob`.
#'
#' (b) probability part (\eqn{S_B} with size \eqn{n_B}; denoted as `var_prob` in the result)
#'
#'  This part uses functionalities of the `{survey}` package and the variance is estimated using the following
#'  equation:
#'
#' \deqn{
#' \hat{V}_2=\frac{1}{N^2} \sum_{i=1}^{n_B} \sum_{j=1}^{n_B} \frac{\pi_{i j}-\pi_i \pi_j}{\pi_{i j}}
#' \frac{m(\boldsymbol{x}_i; \hat{\boldsymbol{\beta}})}{\pi_i} \frac{m(\boldsymbol{x}_i; \hat{\boldsymbol{\beta}})}{\pi_j}.
#' }
#'
#' Note that \eqn{\hat{V}_2} in principle can be estimated in various ways depending on the type of the design and whether population size is known or not.
#'
#' @param y_nons target variable from non-probability sample
#' @param X_nons a `model.matrix` with auxiliary variables from non-probability sample
#' @param X_rand a `model.matrix` with auxiliary variables from non-probability sample
#' @param svydesign a svydesign object
#' @param weights case / frequency weights from non-probability sample
#' @param family_outcome family for the glm model
#' @param start_outcome start parameters
#' @param vars_selection whether variable selection should be conducted
#' @param pop_totals a place holder (not used in `method_pmm`)
#' @param pop_size population size from the `nonprob` function
#' @param control_outcome controls passed by the `control_out` function
#' @param control_inference controls passed by the `control_inf` function
#' @param verbose parameter passed from the main `nonprob` function
#' @param se whether standard errors should be calculated
#'
#' @returns an `nonprob_method` class which is a `list` with the following entries
#'
#' \describe{
#'   \item{model_fitted}{fitted model either an `glm.fit` or `cv.ncvreg` object}
#'   \item{y_nons_pred}{predicted values for the non-probablity sample}
#'   \item{y_rand_pred}{predicted values for the probability sample or population totals}
#'   \item{coefficients}{coefficients for the model (if available)}
#'   \item{svydesign}{an updated `surveydesign2` object (new column `y_hat_MI` is added)}
#'   \item{y_mi_hat}{estimated population mean for the target variable}
#'   \item{vars_selection}{whether variable selection was performed}
#'   \item{var_prob}{variance for the probability sample component (if available)}
#'   \item{var_nonprob}{variance for the non-probability sampl component}
#'   \item{model}{model type (character `"pmm"`)}
#'   \item{family}{depends on the method selected for estimating E(Y|X)}
#' }
#'
#' @examples
#'
#' data(admin)
#' data(jvs)
#' jvs_svy <- svydesign(ids = ~ 1,  weights = ~ weight, strata = ~ size + nace + region, data = jvs)
#'
#' res_pmm <- method_pmm(y_nons = admin$single_shift,
#'                       X_nons = model.matrix(~ region + private + nace + size, admin),
#'                       X_rand = model.matrix(~ region + private + nace + size, jvs),
#'                       svydesign = jvs_svy)
#'
#' res_pmm
#'
#' @export
method_pmm <- function(y_nons,
                       X_nons,
                       X_rand,
                       svydesign,
                       weights=NULL,
                       family_outcome="gaussian",
                       start_outcome=NULL,
                       vars_selection=FALSE,
                       pop_totals=NULL,
                       pop_size=NULL,
                       control_outcome=control_out(),
                       control_inference=control_inf(),
                       verbose=FALSE,
                       se=TRUE) {

  ## passing arguments to the specified method of estimation E(Y|X)
  method_results <- switch(control_outcome$pmm_reg_engine,
                           "glm" = method_glm(y_nons=y_nons,
                                              X_nons=X_nons,
                                              X_rand=X_rand,
                                              svydesign=svydesign,
                                              weights=weights,
                                              family_outcome=family_outcome,
                                              start_outcome=start_outcome,
                                              vars_selection=vars_selection,
                                              pop_size=pop_size,
                                              control_outcome=control_outcome,
                                              control_inference=control_inference,
                                              verbose=verbose,
                                              se=FALSE),
                           "loess" = method_npar(y_nons=y_nons,
                                                 X_nons=X_nons,
                                                 X_rand=X_rand,
                                                 svydesign=svydesign,
                                                 weights=weights,
                                                 family_outcome=family_outcome,
                                                 vars_selection=vars_selection,
                                                 pop_size=pop_size,
                                                 control_outcome=control_outcome,
                                                 control_inference=control_inference,
                                                 verbose=verbose,
                                                 se=FALSE))

  ## Disable `nn_exact_se` during the estimation for PMM
  control_inference_ <- control_inference
  if (isTRUE(control_inference_$nn_exact_se))  control_inference_$nn_exact_se <- FALSE

  ## passing results to method_nn depending on how the matching should be done
  ## 1 - yhat - yhat matching (y_nons_pred, y_rand_pred)
  ## 2 - yhat - y matching (y_nons, y_rand_pred)
  pmm_results <- switch(control_outcome$pmm_match_type,
                        "1" = method_nn(y_nons=y_nons,
                                        X_nons=method_results$y_nons_pred,
                                        X_rand=method_results$y_rand_pred,
                                        svydesign=svydesign,
                                        weights=weights,
                                        vars_selection=vars_selection,
                                        pop_size=pop_size,
                                        control_outcome=control_outcome,
                                        control_inference=control_inference_,
                                        verbose=verbose,
                                        se=se),
                        "2" =  method_nn(y_nons=y_nons,
                                         X_nons=y_nons,
                                         X_rand=method_results$y_rand_pred,
                                         svydesign=svydesign,
                                         weights=weights,
                                         vars_selection=vars_selection,
                                         pop_size=pop_size,
                                         control_outcome=control_outcome,
                                         control_inference=control_inference_,
                                         verbose=verbose,
                                         se=se))


  if (se) {
    svydesign_mean <- svymean(~y_hat_MI, pmm_results$svydesign)
    var_prob <- as.vector(attr(svydesign_mean, "var"))
    var_nonprob <- 0

    if (control_inference$nn_exact_se) {

      message("The `nn_exact_se=TRUE` option used. Remember to set the seed for reproducibility.")

      loop_size <- 50

      if (verbose) {
        message("Estimating nonprob variance component using mini-bootstrap...")
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

        method_results_boot <- switch(control_outcome$pmm_reg_engine,
                                 "glm" = method_glm(y_nons=y_nons_b,
                                                    X_nons=X_nons_b,
                                                    X_rand=X_rand,
                                                    svydesign=svydesign,
                                                    weights=weights,
                                                    family_outcome=family_outcome,
                                                    start_outcome=start_outcome,
                                                    vars_selection=vars_selection,
                                                    pop_size=pop_size,
                                                    control_outcome=control_outcome,
                                                    control_inference=control_inference,
                                                    verbose=FALSE,
                                                    se=FALSE),
                                 "loess" = method_npar(y_nons=y_nons_b,
                                                       X_nons=X_nons_b,
                                                       X_rand=X_rand,
                                                       svydesign=svydesign,
                                                       weights=weights,
                                                       family_outcome=family_outcome,
                                                       vars_selection=vars_selection,
                                                       pop_size=pop_size,
                                                       control_outcome=control_outcome,
                                                       control_inference=control_inference,
                                                       verbose=FALSE,
                                                       se=FALSE))

        pmm_results_boot <- switch(control_outcome$pmm_match_type,
                              "1" = method_nn(y_nons=y_nons_b,
                                              X_nons=method_results_boot$y_nons_pred,
                                              X_rand=method_results_boot$y_rand_pred,
                                              svydesign=svydesign,
                                              weights=weights,
                                              vars_selection=vars_selection,
                                              pop_size=pop_size,
                                              control_outcome=control_outcome,
                                              control_inference=control_inference_,
                                              verbose=FALSE,
                                              se=FALSE),
                              "2" =  method_nn(y_nons=y_nons_b,
                                               X_nons=y_nons_b,
                                               X_rand=method_results_boot$y_rand_pred,
                                               svydesign=svydesign,
                                               weights=weights,
                                               vars_selection=vars_selection,
                                               pop_size=pop_size,
                                               control_outcome=control_outcome,
                                               control_inference=control_inference_,
                                               verbose=FALSE,
                                               se=FALSE))

        dd[jj] <- pmm_results_boot$y_mi_hat

      }
      var_nonprob <- var(dd)

    }
    var_total <- var_prob + var_nonprob

    } else {
    var_prob <- var_nonprob <- var_total <- NA
  }

  return(
    structure(
      list(
        model_fitted = method_results$model_fitted,
        y_nons_pred = method_results$y_nons_pred,
        y_rand_pred = method_results$y_rand_pred,
        coefficients = NULL,
        svydesign = pmm_results$svydesign,
        y_mi_hat = pmm_results$y_mi_hat,
        vars_selection = vars_selection,
        var_prob = var_prob,
        var_nonprob = var_nonprob,
        var_total = var_total,
        model = "pmm",
        family = pmm_results$family
      ), class = "nonprob_method")
  )
}
