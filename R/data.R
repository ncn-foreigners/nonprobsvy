#' Job Vacancy Survey
#'
#' A subset of Job Vacancy Survey from Poland containing 6,523 units and 6 columns
#'
#' \itemize{
#'   \item \code{id} identifier of an entity (company: legal or local).
#'   \item \code{private} whether the company is a private (1) or public (0) entity.
#'   \item \code{size} the size of the entity: S -- small (to 9 employees), M -- medium (10-49) or L -- large (over 49).
#'   \item \code{nace} the main NACE code for a given entity: C, D.E, F, G, H, I, J, K.L, M, N, O, P, Q or R.S (14 levels, 3 combined: D and E, K and L, and R and S).
#'   \item \code{region} the region of Poland (16 levels: 02, 04, ..., 32).
#'   \item \code{weight} the final (calibrated) weight (w-weight). We do not have access to design weights (d-weights).
#' }
#'
#' @docType data
#' @keywords datasets
#' @name jvs
#' @rdname jvs
#' @format A single data.frame with 6,523 rows and 6 columns
"jvs"

#' Admin data (non-probability survey)
#'
#' A subset of Central Job Vacancy Database: a voluntary admin data on
#'
#' \itemize{
#'   \item \code{id} identifier of an entity (company: legal or local).
#'   \item \code{private} whether the company is a private (1) or public (0) entity.
#'   \item \code{size} the size of the entity: S -- small (to 9 employees), M -- medium (10-49) or L -- large (over 49).
#'   \item \code{nace} the main NACE code for a given entity: C, D.E, F, G, H, I, J, K.L, M, N, O, P, Q or R.S (14 levels, 3 combined: D and E, K and L, and R and S).
#'   \item \code{region} the region of Poland (16 levels: 02, 04, ..., 32).
#'   \item \code{single_shift} whether an entity seeks employees on a single shift.
#' }
#'
#' @docType data
#' @keywords datasets
#' @name admin
#' @rdname admin
#' @format A single data.frame with 9,344 rows and 6 columns
"admin"
