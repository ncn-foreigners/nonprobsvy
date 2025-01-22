#' Job Vacancy Survey
#'
#' @description
#' This is a subset of the Job Vacancy Survey from Poland (for one quarter).
#' The data has been subject to slight manipulation, but the relationships in the data have been preserved.
#' For further details on the JVS, please refer to the following link:
#' \url{https://stat.gov.pl/obszary-tematyczne/rynek-pracy/popyt-na-prace/zeszyt-metodologiczny-popyt-na-prace,3,1.html}.
#'
#'
#' @format A single data.frame with 6,523 rows and 6 columns
#'
#' \describe{
#'   \item{\code{id}}{Identifier of an entity (company: legal or local).}
#'   \item{\code{private}}{Whether the company is a private (1) or public (0) entity.}
#'   \item{\code{size}}{The size of the entity: S -- small (to 9 employees), M -- medium (10-49) or L -- large (over 49).}
#'   \item{\code{nace}}{The main NACE code for a given entity: C, D.E, F, G, H, I, J, K.L, M, N, O, P, Q or R.S (14 levels, 3 combined: D and E, K and L, and R and S).}
#'   \item{\code{region}}{The region of Poland (16 levels: 02, 04, ..., 32).}
#'   \item{\code{weight}}{The final (calibrated) weight (w-weight). We do not have access to design weights (d-weights).}
#' }
#'
#' @docType data
#' @keywords datasets
#' @name jvs
#' @rdname jvs
#' @examples
#'
#' data("jvs")
#' head(jvs)
#'
"jvs"

#' Admin data (non-probability survey)
#'
#' @description
#' This is a subset of the Central Job Offers Database, a voluntary administrative data set (non-probability sample).
#' The data was slightly manipulated to ensure the relationships were preserved, and then aligned.
#' For more information about the CBOP, please refer to: \url{https://oferty.praca.gov.pl/}.
#'
#' @format A single data.frame with 9,344 rows and 6 columns
#'
#' \describe{
#'   \item{\code{id}}{Identifier of an entity (company: legal or local).}
#'   \item{\code{private}}{Whether the company is a private (1) or public (0) entity.}
#'   \item{\code{size}}{The size of the entity: S -- small (to 9 employees), M -- medium (10-49) or L -- large (over 49).}
#'   \item{\code{nace}}{The main NACE code for a given entity: C, D.E, F, G, H, I, J, K.L, M, N, O, P, Q or R.S (14 levels, 3 combined: D and E, K and L, and R and S).}
#'   \item{\code{region}}{The region of Poland (16 levels: 02, 04, ..., 32).}
#'   \item{\code{single_shift}}{Whether an entity seeks employees on a single shift.}
#' }
#'
#' @docType data
#' @keywords datasets
#' @name admin
#' @rdname admin
#' @examples
#'
#' data("admin")
#' head(admin)
"admin"
