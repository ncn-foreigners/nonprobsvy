#' @title Plots the estimated mean(s) and their confidence interval(s)
#' @description
#' Simple plotting method that compares the estimated mean(s) and CI(s) with the naive (uncorrected) estimates.
#'
#' @param x the \code{nonprob} class object
#' @param ... other arguments passed to the plot method (currently not supported)
#'
#' @method plot nonprob
#'
#' @examples
#'
#' data(admin)
#' data(jvs)
#'
#' jvs_svy <- svydesign(ids = ~ 1,  weights = ~ weight,
#' strata = ~ size + nace + region, data = jvs)
#'
#' ipw_est1 <- nonprob(selection = ~ region + private + nace + size,
#' target = ~ single_shift,
#' svydesign = jvs_svy,
#' data = admin, method_selection = "logit")
#'
#' plot(ipw_est1)
#'
#' @importFrom graphics arrows
#' @importFrom graphics axis
#' @importFrom graphics grid
#' @importFrom graphics legend
#' @importFrom graphics points
#'
#' @exportS3Method
plot.nonprob <- function(x, ...) {

  est <- confint(x)
  est$mean <- extract(x)
  est$naive <- unname(sapply(x$y, mean))

  x_positions <- 1:nrow(est)
  y_range <- c(min(c(est$lower_bound), est$naive) * 0.95,
               max(c(est$upper_bound, est$naive)) * 1.05)

  plot(x_positions, est$naive, ylim = y_range,
       xlim = c(0.5, length(x_positions) + 0.5),
       pch = 17, col = "#0072B2", cex = 1.5,
       xlab = "", ylab = "Estimates", xaxt = "n",
       main = "Comparison of estimates")

  axis(1, at = x_positions, labels = est$target)
  points(x_positions, est$mean, pch = 19, col = "#D55E00", cex = 1.5)
  arrows(x_positions, est$lower_bound, x_positions, est$upper_bound,
         code = 3, angle = 90, length = 0.1, lwd = 2, col = "#D55E00")

  legend("topright",
         legend = c("Mean with 95% CI", "Naive"),
         pch = c(19, 17),
         col = c("#D55E00", "#0072B2"),
         pt.cex = 1.5,
         bty = "n")

  grid(nx = NA, ny = NULL, lty = 2, col = "gray")

}


