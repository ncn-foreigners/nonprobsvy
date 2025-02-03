#'
#' Internal function to create an object with model frames and matrices to preprocess
#'
#' @param formula a formula for the model
#' @param data data frame
#' @param weights weights for the data
#' @param svydesign survey design object
#' @param pop_totals population totals
#' @param pop_size population size
#' @param flag logical value indicating whether to check if the names of the variables in the data and in svydesign match
#'
#' @return list with model frames and matrices to preprocess
#' @noRd
model_frame <- function(formula,
                        data,
                        weights = NULL,
                        svydesign = NULL,
                        pop_totals = NULL,
                        pop_size = NULL,
                        flag = TRUE) {
  if (!is.null(svydesign)) {

    ##### Model frame for non-probability sample #####

    model_frame_nons <- model.frame(formula, data)
    y_nons <- model.response(model_frame_nons)
    outcome_name <- names(model_frame_nons)[1]

    mt <- terms(formula) # attr(model_frame_nons, "terms")
    nons_names <- all.vars(as.formula(paste("~", paste(attr(mt, "term.labels"), collapse = " + "))))
    ##### Model frame for probability sample #####
    if (outcome_name %in% colnames(svydesign$variables)) {
      # TODO think about implementation below - probablby can be done in a better way
      design_to_frame <- svydesign$variables
      design_to_frame[, outcome_name][is.na(design_to_frame[, outcome_name])] <- 0 # replace NA in dependent outcome with 0
      names_rand <- all.vars(formula)
      model_frame_nons_rand <- design_to_frame[, names_rand, drop = FALSE]

      model_frame_rand <- intersect(names_rand, colnames(design_to_frame))
    } else {
      design_to_frame <- svydesign$variables
      ##
      terms_object <- terms(formula)
      names_rand <- all.vars(as.formula(paste("~", paste(attr(terms_object, "term.labels"), collapse = " + "))))

      model_frame_nons_rand <- design_to_frame[, names_rand, drop = FALSE]
      model_frame_rand <- intersect(names_rand, colnames(design_to_frame))
    }
    # TODO think out this condition
    if (all(nons_names %in% model_frame_rand)) { # colnames(svydesign$variables)
      dot_check <- sapply(formula, FUN = function(x) {
        x == "."
      })
      if (length(formula) == 2) nons_names <- nons_names[-1]
      if (any(dot_check)) {
        xx <- paste("~", paste(nons_names, collapse = "+"))
        formula <- as.formula(paste(outcome_name, xx))
        X_rand <- model.matrix(delete.response(terms(formula)), svydesign$variables[, nons_names])
      } else {
        X_rand <- model.matrix(delete.response(terms(formula)), svydesign$variables) # matrix of probability sample with intercept
      }
      frame_nons <- model.frame(formula, data)
      X_nons <- model.matrix(frame_nons, data) # matrix for nonprobability sample with intercept
      # old version
      # if (outcome) {
      #  xx <- paste("~", paste(nons_names[2:length(nons_names)], collapse = "+"))
      #  formula <- as.formula(paste(formula[2], xx))
      #  X_rand <- model.matrix(delete.response(terms(formula)), svydesign$variables[, nons_names])
      # } else {
      #  xx <- paste("~", paste(nons_names, collapse = "+"))
      #  formula <- as.formula(xx)
      #  X_rand <- model.matrix(formula, svydesign$variables[, nons_names])# matrix of probability sample with intercept
      #  }
    } else {
      stop("The names of the variables in the data and in svydesign do not match.")
    }

    list(
      X_nons = X_nons,
      X_rand = X_rand,
      nons_names = nons_names,
      y_nons = y_nons,
      outcome_name = outcome_name,
      model_frame_rand = model_frame_nons_rand
    )
  } else if (!is.null(pop_totals)) {
    model_frame_nons <- model.frame(formula, data)
    X_nons <- model.matrix(model_frame_nons, data)
    mt <- attr(model_frame_nons, "terms")
    total_names <- colnames(X_nons)
    if (flag) {
      if (all(total_names %in% names(pop_totals))) { # TODO verify whether this warning works well.. pop_totals, pop_means defined such as in `calibrate` function
        pop_totals <- pop_totals[total_names]
      } else {
        warning("Selection and population totals have different names.")
      }
    }
    y_nons <- model.response(model_frame_nons)
    outcome_name <- names(model_frame_nons)[1]

    list(
      X_nons = X_nons,
      pop_totals = pop_totals,
      total_names = total_names,
      y_nons = y_nons,
      outcome_name = outcome_name,
      X_rand = NULL
    )
  }
}
