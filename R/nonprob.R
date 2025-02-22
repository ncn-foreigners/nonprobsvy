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
                    method_selection = c("logit", "cloglog", "probit"),
                    method_outcome = c("glm", "nn", "pmm", "npar"),
                    family_outcome = c("gaussian", "binomial", "poisson"),
                    subset = NULL,
                    strata = NULL,
                    weights = NULL,
                    na_action = NULL,
                    control_selection = control_sel(),
                    control_outcome = control_out(),
                    control_inference = control_inf(),
                    start_selection = NULL,
                    start_outcome = NULL,
                    verbose = FALSE,
                    x = TRUE,
                    y = TRUE,
                    se = TRUE,
                    ...) {

  call <- match.call()
  data <- if (!is.data.frame(data)) data.frame(data) else data
  weights <- if (is.null(weights)) rep(1, nrow(data)) else weights
  num_of_vars <- NROW(attr(terms(outcome), "term.labels"))



  # Method defaults
  method_selection <- if (missing(method_selection)) "logit" else method_selection
  family_outcome <- if (missing(family_outcome)) "gaussian" else family_outcome
  method_outcome <- if (missing(method_outcome)) "glm" else method_outcome

  # Check data
  if (is.data.frame(data) && ncol(data) == 0) {
    stop("The `data` argument cannot be an empty data frame.")
  }

  # Validation checks for methods
  if (!(method_outcome %in% c("glm", "nn", "pmm", "npar"))) {
    stop("Invalid method for the `outcome` variable. Choose from 'glm', 'nn', 'pmm', 'npar'")
  }
  if (method_outcome == "npar" & num_of_vars >= 4) {
    stop("Non-parametric methods are suited for smaller number of variables (less than 4)")
  }

  if (!(method_selection %in% c("logit", "cloglog", "probit"))) {
    stop("Invalid method for the `selection` formula. Choose from 'logit', 'cloglog', 'probit'.")
  }

  if (!(family_outcome %in% c("gaussian", "binomial", "poisson"))) {
    stop("Invalid family for the `outcome` formula. Choose from 'gaussian', 'binomial', 'poisson'.")
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

  # Validation checks for totals and means
  if (!is.null(pop_totals) && !is.vector(pop_totals)) {
    stop("The `pop_totals` argument must be a vector.")
  }
  if (!is.null(pop_means) && !is.vector(pop_means)) {
    stop("The `pop_means` argument must be a vector.")
  }

  if (!is.null(pop_size) && (!is.numeric(pop_size) || pop_size <= 0)) {
    stop("The `pop_size` argument must be a positive numeric scalar.")
  }

  if (!is.null(pop_totals) && !is.null(pop_means)) {
    stop("Specify one of the `pop_totals` or `pop_means` arguments, not both.")
  }
  if (!is.null(pop_size) && pop_size < nrow(data)) {
    stop("The `pop_size` argument cannot be smaller than sample size.")
  }

  ## for weights
  if (!is.null(weights) && !is.numeric(weights)) {
    stop("The `weights` argument must be a numeric vector.")
  }
  if (!is.null(weights) && length(weights) != nrow(data)) {
    stop("Length of the `weights` argument must match the number of rows in data.")
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

  ## this check is for future development
  if (!is.null(control_selection$key)) {
    if (!(control_selection$key %in% colnames(data)) || !(control_selection$key %in% colnames(svydesign$variables))) {
      stop("Key variable for overlapping units must be defined with the same name in prob and nonprob sample.")
    }
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
    }
  } else if (inherits(outcome, "formula")) {
    # Case: MI method
    estimator <- "mi"
  }

  ## processing data for all methods
  if (is.null(pop_totals) && !is.null(pop_means) && !is.null(pop_size)) {
    pop_totals <- c(pop_size, pop_size * pop_means)
    names(pop_totals) <- c("(Intercept)", names(pop_means))
  } else if (!is.null(pop_totals)) {
    pop_size <- unname(pop_totals["(Intercept)"])
  } else {
    pop_size <- sum(weights(svydesign))
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
      subset = subset,
      strata = strata,
      weights = weights,
      na_action = na_action,
      control_selection = control_selection,
      control_inference = control_inference,
      start_selection = start_selection,
      verbose = verbose,
      x = x,
      y = y,
      se = se,
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
      subset = subset,
      strata = strata,
      weights = weights,
      na_action = na_action,
      control_outcome = control_outcome,
      control_inference = control_inference,
      start_outcome = start_outcome,
      verbose = verbose,
      x = x,
      y = y,
      se = se,
      ...
    ),
    dr = nonprob_dr(
      selection = selection,
      outcome = outcome,
      data = data,
      svydesign = svydesign,
      pop_totals = pop_totals,
      pop_means = pop_means,
      pop_size = pop_size,
      method_selection = method_selection,
      method_outcome = method_outcome,
      family_outcome = family_outcome,
      subset = subset,
      strata = strata,
      weights = weights,
      na_action = na_action,
      control_selection = control_selection,
      control_outcome = control_outcome,
      control_inference = control_inference,
      start_selection = start_selection,
      start_outcome = start_outcome,
      verbose = verbose,
      x = x,
      y = y,
      se = se,
      ...
    )
  )
  names <- names(model_estimates)
  res <- append(model_estimates, call, after = 0)
  names(res) <- c("call", names)

  res$estimator <- estimator

  return(structure(res,
                   class = class(model_estimates)))
}
