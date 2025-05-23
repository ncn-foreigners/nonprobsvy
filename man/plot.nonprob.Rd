% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plot.R
\name{plot.nonprob}
\alias{plot.nonprob}
\title{Plots the estimated mean(s) and their confidence interval(s)}
\usage{
\method{plot}{nonprob}(x, ...)
}
\arguments{
\item{x}{the \code{nonprob} class object}

\item{...}{other arguments passed to the plot method (currently not supported)}
}
\description{
Simple plotting method that compares the estimated mean(s) and CI(s) with the naive (uncorrected) estimates.
}
\examples{

data(admin)
data(jvs)

jvs_svy <- svydesign(ids = ~ 1,  weights = ~ weight,
strata = ~ size + nace + region, data = jvs)

ipw_est1 <- nonprob(selection = ~ region + private + nace + size,
target = ~ single_shift,
svydesign = jvs_svy,
data = admin, method_selection = "logit")

plot(ipw_est1)

}
