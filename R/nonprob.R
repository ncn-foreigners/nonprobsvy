#' @rdname nonprob
#' @export
nonprob <- function(data,
                    selection = NULL,
                    outcome = NULL,
                    target = NULL,
                    svydesign = NULL,
                    pop_totals = NULL,
                    pop_means = NULL,
                    pop_size = NULL,
                    overlap = FALSE,
                    method_selection = c("logit", "cloglog", "probit"),
                    method_outcome = c("glm", "nn"),
                    family_selection = "binomial",
                    family_outcome = c("gaussian", "binomial", "poisson"),
                    subset = NULL,
                    strata = NULL,
                    weights = NULL,
                    na_action = NULL,
                    control_selection = controlSel(),
                    control_outcome = controlOut(),
                    control_inference = controlInf(),
                    start = NULL,
                    verbose = FALSE,
                    x = TRUE,
                    y = TRUE,
                     ...){

  call <- match.call()

  var_selection <- control_inference$vars_selection

  if (!is.data.frame(data)) {
    data <- data.frame(data)
  }

  if (is.null(weights)) weights <- rep(1, nrow(data))

  if(missing(method_selection)) method_selection <- "logit"
  if(missing(family_outcome)) family_outcome <- "gaussian"
  if(missing(method_outcome)) method_outcome <- "glm"
  if(!(method_outcome %in% c("glm", "nn"))) stop("Invalid method for outcome variable.")

  if(!(method_selection %in% c("logit", "cloglog", "probit"))) stop("Invalid method for selection formula.")
  if(!(family_outcome %in% c("gaussian", "binomial", "poisson"))) stop("Invalid family for outcome formula.")
  if (!is.null(control_selection$key)) {
    if (!(control_selection$key %in% colnames(data)) || !(control_selection$key %in% colnames(svydesign$variables))) {
      stop("key variable for overlapping units must be defined with this same name in prob and nonprob sample.")
    }
  }

  ## basic checkers
  if (is.null(selection) & is.null(outcome)) {
    stop("Please provide selection or outcome formula.")
  }
  if (inherits(selection, "formula") && inherits(target, "formula") && (is.null(outcome) || inherits(outcome, "formula") == FALSE)) {
    ifelse(var_selection == FALSE, model_used <- "P", model_used <- "Psel")
  }

  if (inherits(outcome, "formula") && (is.null(selection) || inherits(selection, "formula") == FALSE)) {
    ifelse(var_selection == FALSE, model_used <- "M", model_used <- "Msel")
  }

  if (inherits(selection, "formula") && inherits(outcome, "formula")) {

    ifelse(var_selection == FALSE, model_used <- "DR", model_used <- "DRsel")
  }

  ## validate data

  ## model estimates
  model_estimates <- switch(model_used,
    P = nonprobIPW(selection,
                   target,
                   data,
                   svydesign,
                   pop_totals,
                   pop_means,
                   pop_size,
                   overlap,
                   method_selection,
                   family_selection,
                   subset,
                   strata,
                   weights,
                   na_action,
                   control_selection,
                   control_inference,
                   start,
                   verbose,
                   x,
                   y,
                   ...),
    M = nonprobMI(outcome,
                  data,
                  svydesign,
                  pop_totals,
                  pop_means,
                  pop_size,
                  method_outcome,
                  family_outcome,
                  subset,
                  strata,
                  weights,
                  na_action,
                  control_outcome,
                  control_inference,
                  start,
                  verbose,
                  x,
                  y,
                  ...),
    DR = nonprobDR(selection,
                   outcome,
                   data,
                   svydesign,
                   pop_totals,
                   pop_means,
                   pop_size,
                   overlap,
                   method_selection,
                   method_outcome,
                   family_selection,
                   family_outcome,
                   subset,
                   strata,
                   weights,
                   na_action,
                   control_selection,
                   control_outcome,
                   control_inference,
                   start,
                   verbose,
                   x,
                   y,
                   ...),
    Psel = nonprobSelP(selection,
                      target,
                      data,
                      svydesign,
                      pop_totals,
                      pop_means,
                      pop_size,
                      method_selection,
                      method_outcome,
                      family_selection,
                      family_outcome,
                      subset,
                      strata,
                      weights,
                      na_action,
                      control_selection,
                      control_outcome,
                      control_inference,
                      start,
                      verbose,
                      x,
                      y,
                     ...),
    Msel = nonprobSelM(outcome,
                      data,
                      svydesign,
                      pop_totals,
                      pop_means,
                      pop_size,
                      method_outcome,
                      family_outcome,
                      subset,
                      strata,
                      weights,
                      na_action,
                      control_outcome,
                      control_inference,
                      start,
                      verbose,
                      x,
                      y,
                      ...),
    DRsel = nonprobSel(selection,
                       outcome,
                       data,
                       svydesign,
                       pop_totals,
                       pop_means,
                       pop_size,
                       method_selection,
                       method_outcome,
                       family_selection,
                       family_outcome,
                       subset,
                       strata,
                       weights,
                       na_action,
                       control_selection,
                       control_outcome,
                       control_inference,
                       start,
                       verbose,
                       x,
                       y,
                       ...)
  )

  names <- names(model_estimates)
  res <- append(model_estimates, call, after = 0)
  names(res) <- c("call", names)
  structure(res, class = class(model_estimates))
}
