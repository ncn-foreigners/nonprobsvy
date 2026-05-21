#' @rdname nonprob
#' @importFrom stats na.omit
#' @export
nonprob <- function(data,
                    selection = NULL,
                    outcome = NULL,
                    target = NULL,
                    svydesign = NULL,
                    pop_totals = NULL,
                    pop_means = NULL,
                    pop_size = NULL,
                    method_selection = c("logit", "cloglog", "probit"),
                    method_outcome = c("glm", "nn", "pmm", "npar"),
                    family_outcome = c("gaussian", "binomial", "poisson"),
                    subset = NULL,
                    strata = NULL,
                    case_weights = NULL,
                    na_action = na.omit,
                    control_selection = control_sel(),
                    control_outcome = control_out(),
                    control_inference = control_inf(),
                    start_selection = NULL,
                    start_outcome = NULL,
                    verbose = FALSE,
                    se = TRUE,
                    ...) {

  call <- match.call()
  method_selection <- match.arg(method_selection)
  family_outcome <- match.arg(family_outcome)
  method_outcome <- match.arg(method_outcome)
  #stopifnot("We currently support `family_outcome` with a single entry." = NROW(family_outcome) == 1)

  if (!identical(na_action, na.omit)) {
    stop("We currently support only `na_action=na.omit`.")
  }
  data <- if (!is.data.frame(data)) data.frame(data) else data
  data <- subset(data, subset = if (is.null(subset)) TRUE else subset)
  data <- na_action(data)

  case_weights <- if (is.null(case_weights)) rep(1, nrow(data)) else case_weights

  # Check data
  if (is.data.frame(data) && ncol(data) == 0) {
    stop("The `data` argument cannot be an empty data frame.")
  }

  if (method_outcome == "npar") {
    num_of_vars <- NROW(attr(terms(outcome), "term.labels"))
    if (num_of_vars >= 4) stop("Non-parametric methods are suited for smaller number of variables (less than 4)")
  }

  # Validation checks for formulas
  if (!is.null(selection) && !inherits(selection, "formula")) {
    stop("The `selection` argument must be a formula.")
  }
  if (!is.null(outcome) && !inherits(outcome, "formula")) {
    stop("The `outcome` argument must be a formula.")
  }
  if (!is.null(target) && !inherits(target, "formula")) {
    stop("The `target` argument must be a formula.")
  }

  if (!is.null(outcome) & !is.null(selection) & method_outcome != "glm") {
    stop("Currently we only support `method_glm` for doubly robust estimators.")
  }

  # Validation checks for totals and means
  if (!is.null(pop_totals) && !is.vector(pop_totals)) {
    stop("The `pop_totals` argument must be a vector.")
  }
  if (!is.null(pop_means) && !is.vector(pop_means)) {
    stop("The `pop_means` argument must be a vector.")
  }

  if (!is.null(pop_size) &&
      (!is.numeric(pop_size) || length(pop_size) != 1 ||
       is.na(pop_size) || !is.finite(pop_size) || pop_size <= 0)) {
    stop("The `pop_size` argument must be a positive numeric scalar.")
  }

  if (!is.null(pop_totals) && !is.null(pop_means)) {
    stop("Specify one of the `pop_totals` or `pop_means` arguments, not both.")
  }
  if (!is.null(pop_means) && is.null(pop_size)) {
    stop("The `pop_size` argument must be supplied when `pop_means` is supplied.")
  }
  if (!is.null(pop_size) && pop_size < nrow(data)) {
    stop("The `pop_size` argument cannot be smaller than sample size.")
  }

  ## for case_weights
  if (!is.null(case_weights) && !is.numeric(case_weights)) {
    stop("The `case_weights` argument must be a numeric vector.")
  }
  if (!is.null(case_weights) && length(case_weights) != nrow(data)) {
    stop("Length of the `case_weights` argument must match the number of rows in data.")
  }
  if (!is.null(case_weights) && anyNA(case_weights)) {
    stop("The `case_weights` argument cannot contain missing values.")
  }
  if (!is.null(case_weights) && !all(is.finite(case_weights))) {
    stop("The `case_weights` argument must contain only finite values.")
  }
  if (!is.null(case_weights) && any(case_weights <= 0)) {
    stop("The `case_weights` argument must contain only positive values.")
  }

  ## selection and outcome should be specified
  if (is.null(selection) && is.null(outcome)) {
    stop("Please provide the `selection` or `outcome` argument.")
  }

  if (!is.null(svydesign)) {
    if ("svyrep.design" %in% class(svydesign)) {
      stop("We do not currently support the `svyrep.design` class. Provide the survey data in the `survey.design2` class.")
    }
    if ("pps" %in% class(svydesign)) {
      stop("The `as.svrepdesign` function does not allow `pps` designs. For more details, see the `survey` package.")
    }
  }

  ## overlap handling between the two samples is not yet implemented
  if (isTRUE(control_selection$dependence) || !is.null(control_selection$key)) {
    stop("Handling of overlapping units between the samples (`control_sel(dependence = TRUE)` / `key = ...`) is not yet implemented.")
  }


  ## specification of the methods and validation
  if (inherits(selection, "formula")) {
    if (is.null(outcome) || !inherits(outcome, "formula")) {
      # Case: IPW method
      if (!inherits(target, "formula")) {
        stop("Please provide the `target` argument as a formula.")
      }

      # Extract variables from target and selection formulas
      target_vars <- all.vars(target)
      selection_vars <- all.vars(selection)

      # Check for overlapping variables
      common_vars <- intersect(target_vars, selection_vars)
      if (length(common_vars) > 0) {
        stop(sprintf(
          "The following variables in target should not be present in selection: %s",
          paste(common_vars, collapse = ", ")
        ))
      }
      estimator <- "ipw"
    } else {
      # Case: DR method
      estimator <- "dr"

      if (isTRUE(control_inference$bias_correction) & isFALSE(control_inference$vars_combine) & isTRUE(control_inference$vars_selection)) {
        stop("Bias correction (joint estimation) with variable selection can only be used when variables
         are combined after they are selected (change `vars_combine` to `TRUE`).")
      }

      if (isTRUE(control_inference$bias_correction) & is.null(svydesign)) {
        stop("Bias correction (joint estimation) requires individual-level probability sample data;
         supply a `svydesign` argument instead of `pop_totals` / `pop_means`.")
      }

    }
  } else if (inherits(outcome, "formula")) {
    # Case: MI method
    estimator <- "mi"
  }

  ## outcome/target variable(s) must be numeric or logical: factor/character
  ## outcomes silently misbehave (e.g. the estimate is returned as NA). Checked
  ## after the estimator-specific formula validation so that more specific
  ## messages (e.g. target/selection overlap) take precedence.
  response_vars <- unique(c(
    if (!is.null(target)) all.vars(target),
    if (!is.null(outcome) && length(outcome) == 3L) all.vars(outcome[[2L]])
  ))
  for (.response in response_vars) {
    if (.response %in% names(data) &&
        (is.factor(data[[.response]]) || is.character(data[[.response]]))) {
      stop(sprintf(
        "The outcome/target variable `%s` is a %s; only numeric or logical outcomes are supported. Convert it before calling `nonprob()` (e.g. to 0/1 for `family_outcome = \"binomial\"`).",
        .response, if (is.factor(data[[.response]])) "factor" else "character"))
    }
  }

  if (!is.null(pop_totals) && (is.null(names(pop_totals)) || names(pop_totals)[1] != "(Intercept)")) {
    stop("`pop_totals` must be a named numeric vector whose first element is named `(Intercept)` (the population size N), followed by the covariate population totals.")
  }

  if (isTRUE(control_inference$vars_selection) && !is.null(pop_totals) && is.null(svydesign) && estimator %in% c("ipw", "dr")) {
    warning("Variable selection for the selection model with population totals only (no `svydesign`) is experimental. ",
            "The cross-validation used to pick the penalty (`lambda`) is unreliable in this setting and may select an ",
            "over-sparse model (dropping informative covariates), which biases the point estimate. ",
            "Prefer supplying a unit-level probability sample via `svydesign`, or fix the penalty with `control_sel(lambda = ...)`.")
  }

  pop_size_fixed <- !is.null(pop_size) | (!is.null(pop_totals) && names(pop_totals)[1] == "(Intercept)")  ## for variance estimation

  ## processing data for all methods
  if (is.null(pop_totals) && !is.null(pop_means) && !is.null(pop_size)) {
    pop_totals <- c(pop_size, pop_size * pop_means)
    names(pop_totals) <- c("(Intercept)", names(pop_means))
  } else if (!is.null(pop_totals)) {
    pop_size <- unname(pop_totals["(Intercept)"])
  } else if (is.null(pop_size)) {
    ## estimated population size (we may consider for future to use fpc of the svydesing2 object)
    pop_size <- sum(weights(svydesign))
  }

  ## merge formulas prior dr if user wants joint

  if (estimator == "dr" & control_inference$vars_combine & !control_inference$vars_selection) {
    dr_formulas <- merge_formulas(outcome, selection)
    outcome <- dr_formulas[[1]]
    selection <- dr_formulas[[2]]
  }

  ## model estimates
  model_estimates <- switch(estimator,
    ipw = nonprob_ipw(
      selection = selection,
      target = target,
      data = data,
      svydesign = svydesign,
      pop_totals = pop_totals,
      pop_means = pop_means,
      pop_size = pop_size,
      method_selection = method_selection,
      strata = strata,
      case_weights = case_weights,
      na_action = na_action,
      control_selection = control_selection,
      control_inference = control_inference,
      start_selection = start_selection,
      verbose = verbose,
      se = se,
      pop_size_fixed=pop_size_fixed,
      ...
    ),
    mi = nonprob_mi(
      outcome = outcome,
      data = data,
      svydesign = svydesign,
      pop_totals = pop_totals,
      pop_means = pop_means,
      pop_size = pop_size,
      method_outcome = method_outcome,
      family_outcome = family_outcome,
      strata = strata,
      case_weights = case_weights,
      na_action = na_action,
      control_outcome = control_outcome,
      control_inference = control_inference,
      start_outcome = start_outcome,
      verbose = verbose,
      se = se,
      pop_size_fixed=pop_size_fixed,
      ...
    ),
    dr = nonprob_dr(
      selection = selection,
      outcome = outcome,
      data = data,
      svydesign = svydesign,
      pop_totals = pop_totals,
      pop_means = pop_means,
      pop_size = unname(pop_size),
      method_selection = method_selection,
      method_outcome = method_outcome,
      family_outcome = family_outcome,
      strata = strata,
      case_weights = case_weights,
      na_action = na_action,
      control_selection = control_selection,
      control_outcome = control_outcome,
      control_inference = control_inference,
      start_selection = start_selection,
      start_outcome = start_outcome,
      verbose = verbose,
      se = se,
      pop_size_fixed=pop_size_fixed,
      ...
    )
  )
  names <- names(model_estimates)
  res <- append(model_estimates, call, after = 0)
  names(res) <- c("call", names)

  res$estimator <- estimator
  res$outcome_formula <- outcome
  res$selection_formula <- selection

  ## add information about the type of the model
  res$estimator_method <- if (!is.null(outcome)) {
    ifelse(method_outcome == "nn", method_outcome, paste0(method_outcome, " (", family_outcome,")"))
    } else paste0(method_selection, " (", control_selection$est_method,")")



  return(structure(res,
                   class = class(model_estimates)))
}
