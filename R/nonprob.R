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
                    method_outcome = c("glm", "nn", "pmm"),
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

  if (!is.data.frame(data)) {
    data <- data.frame(data)
  }

  if (is.null(weights)) weights <- rep(1, nrow(data))

  if (missing(method_selection)) method_selection <- "logit"
  if (missing(family_outcome)) family_outcome <- "gaussian"
  if (missing(method_outcome)) method_outcome <- "glm"
  if (!(method_outcome %in% c("glm", "nn", "pmm"))) stop("Invalid method for the `outcome` variable.")
  if (!is.null(svydesign)) {
    if ("svyrep.design" %in% class(svydesign)) stop("We do not currently support the `svyrep.design` class. Provide the survey data in the `survey.design2` class.")
    if ("pps" %in% class(svydesign)) stop("The `as.svrepdesign` function does not allow `pps` designs. For more details, see the `survey` package.")
  }
  if (!is.null(pop_totals)) {
    if (!is.vector(pop_totals)) stop("The `pop_totals` argument must be a vector.")
  }
  if (!is.null(pop_means)) {
    if (!is.vector(pop_means)) stop("The `pop_means` argument must be a vector.")
  }
  if (!(method_selection %in% c("logit", "cloglog", "probit"))) stop("Invalid method for the `selection` formula.")
  if (!(family_outcome %in% c("gaussian", "binomial", "poisson"))) stop("Invalid family for the `outcome` formula.")
  if (!is.null(control_selection$key)) {
    if (!(control_selection$key %in% colnames(data)) || !(control_selection$key %in% colnames(svydesign$variables))) {
      stop("Key variable for overlapping units must be defined with the same name in prob and nonprob sample.")
    }
  }

  ## basic checkers
  if (is.null(selection) & is.null(outcome)) {
    stop("Please provide thw `selection` or `outcome` argument.")
  }

  # Check formula inputs
  if (!is.null(selection) && !inherits(selection, "formula")) {
    stop("The `selection` argument must be a formula.")
  }
  if (!is.null(outcome) && !inherits(outcome, "formula")) {
    stop("The `outcome` argument must be a formula.")
  }
  if (!is.null(target) && !inherits(target, "formula")) {
    stop("The `target` argument must be a formula.")
  }

  if (inherits(selection, "formula") && (is.null(outcome) || inherits(outcome, "formula") == FALSE)) {
    if (inherits(target, "formula") == FALSE) stop("Please provide the `target` argument.")
    model_used <- "P"
  }

  if (inherits(outcome, "formula") && (is.null(selection) || inherits(selection, "formula") == FALSE)) {
    model_used <- "M"
  }

  if (inherits(selection, "formula") && inherits(outcome, "formula")) {
    model_used <- "DR"
  }

  # Check numeric inputs
  if (!is.null(pop_size) && !is.numeric(pop_size)) {
    stop("The `pop_size` argument must be numeric scalar.")
  }
  if (!is.null(weights) && !is.numeric(weights)) {
    stop("The `weights` argument must be a numeric vector.")
  }

  # Check weights length
  if (!is.null(weights) && length(weights) != nrow(data)) {
    stop("Length of the `weights` argument must match the number of rows in data.")
  }

  if (!is.null(pop_totals) && !is.null(pop_means)) {
    stop("Specify one of the `pop_totals' or `pop_means' arguments, not both.")
  }

  if (!is.null(pop_size)) {
    if (pop_size <= 0) {
      stop("The `pop_size` argument must be positive.")
    }
    if (pop_size < nrow(data)) {
      stop("The `pop_size` argument cannot be smaller than sample size.")
    }
  }

  ## model estimates
  model_estimates <- switch(model_used,
    P = nonprobIPW(
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
    M = nonprobMI(
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
    DR = nonprobDR(
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
  structure(res, class = class(model_estimates))
}
