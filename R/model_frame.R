# An internal function to create an object with model frames and matrices to preprocess
#
# @param formula a formula for the model
# @param data data frame
# @param weights weights for the data
# @param svydesign survey design object
# @param pop_totals population totals
#
# @return list with model frames and matrices to preprocess
# @noRd
model_frame <- function(formula,
                        data,
                        weights = NULL,
                        svydesign = NULL,
                        pop_totals = NULL) {

  update_formula <- function(formula, nons_names, outcome_name) {
    dot_check <- sapply(formula, function(x) x == ".")
    if (length(formula) == 2)
      nons_names <- nons_names[-1]
    if (any(dot_check)) {
      formula <- as.formula(paste(outcome_name, "~", paste(nons_names, collapse = "+")))
    }
    return(formula)
  }

  model_frame_svydesign <- function(formula, data, svydesign) {
    model_frame_nons <- model.frame(formula, data)
    y_nons <- model.response(model_frame_nons)
    outcome_name <- names(model_frame_nons)[1]

    mt <- terms(formula)
    nons_names <- all.vars(as.formula(paste("~", paste(attr(mt, "term.labels"), collapse = " + "))))

    design_to_frame <- svydesign$variables

    if (outcome_name %in% colnames(design_to_frame)) {
      design_to_frame[, outcome_name][is.na(design_to_frame[, outcome_name])] <- 0
      names_rand <- all.vars(formula)
      model_frame_nons_rand <- design_to_frame[, names_rand, drop = FALSE]
      model_frame_rand <- intersect(names_rand, colnames(design_to_frame))
    } else {
      terms_object <- terms(formula)
      names_rand <- all.vars(as.formula(paste("~", paste(
        attr(terms_object, "term.labels"), collapse = " + "
      ))))
      model_frame_nons_rand <- design_to_frame[, names_rand, drop = FALSE]
      model_frame_rand <- intersect(names_rand, colnames(design_to_frame))
    }

    if (all(nons_names %in% model_frame_rand)) {
        formula <- update_formula(formula, nons_names, outcome_name)
        X_rand <- model.matrix(delete.response(terms(formula)), svydesign$variables[, nons_names, drop = FALSE])
        frame_nons <- model.frame(formula, data)
        X_nons <- model.matrix(frame_nons, data)
      } else {
        stop("The names of the variables in the data and in svydesign do not match.")
    }


    list(
      X_nons = X_nons,
      X_rand = X_rand,
      y_nons = y_nons,
      pop_totals = NULL,
      nons_names = nons_names,
      total_names = NULL,
      outcome_name = outcome_name,
      model_frame_rand = model_frame_nons_rand
    )
  }

  model_frame_pop <- function(formula, data, pop_totals) {
    model_frame_nons <- model.frame(formula, data)
    X_nons <- model.matrix(model_frame_nons, data)
    total_names <- colnames(X_nons)

    if (all(total_names %in% names(pop_totals))) {
        pop_totals <- pop_totals[total_names]
      } else {
        stop("Selection and population totals have different names.")
      }

    y_nons <- model.response(model_frame_nons)
    outcome_name <- names(model_frame_nons)[1]

    list(
      X_nons = X_nons,
      X_rand = NULL,
      y_nons = y_nons,
      pop_totals = pop_totals,
      nons_names = NULL,
      total_names = total_names,
      outcome_name = outcome_name,
      model_frame_rand = NULL
    )
  }

  if (!is.null(svydesign)) {
    return(model_frame_svydesign(formula, data, svydesign, flag))
  } else if (!is.null(pop_totals)) {
    return(model_frame_pop(formula, data, pop_totals, flag))
  } else {
    stop("Either svydesign or pop_totals must be provided.")
  }
}
