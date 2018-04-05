\name{coxrw}
\alias{coxrw}
\alias{print.coxr}
\alias{predict.coxr}

\title{Fit Robustly Proportional Hazards Regression Model, with Sampling Weights}

\description{
Fits efficiently and robustly Cox proportional hazards regression model in its
basic form, where explanatory variables are time independent with one event
per subject.  Method is based on a smooth modification of the partial
likelihood. Option to include sampling weights was added to original coxr function.
}

\usage{

coxrw(formula, data, subset, na.action, trunc = 0.95,
     f.weight = c("linear", "quadratic", "exponential"),
     singular.ok = TRUE, model = FALSE, weights)

}

\arguments{

\item{formula}{
a formula object, with the response on the left of a \code{~}
operator, and the terms on the right.  The response must be a
survival object as returned by the \code{\link{Surv}} function.}

\item{data}{a data frame in which to interpret the variables
named in the \code{formula}, or in the \code{subset}.}

\item{subset}{expression saying that only a subset of the rows of the data
should be used in the fit.}

\item{na.action}{a missing-data filter function, applied to the model.frame,
after any subset argument has been used.}

\item{trunc}{roughtly, quantile of the sample \eqn{T_i exp(\beta'Z_i)},
it determines the trimming level for the robust estimator}

\item{f.weight}{type of weighting function, default is \code{"quadratic"}}

\item{singular.ok}{logical value indicating how to handle collinearity in the
model matrix. If \code{TRUE}, the program will automatically skip
over columns of the X matrix that are linear combinations of
earlier columns.  In this case the coefficients for such
columns will be \code{NA}, and the variance matrix will contain
zeros.  For ancillary calculations, such as the linear
predictor, the missing coefficients are treated as zeros.}

\item{model}{a logical value indicating whether model frame should be
included as a component of the returned value.}

\item{weights}{a vector of sampling weights to be used in fitting the model, with length equal to the number of observations in the data.}
}

\details{
The method consists in maximization of an objective function which is
a smooth modification of the partial likelihood. Observations with excessive
values of \eqn{\Lambda(T) exp(\beta'Z)}, where \eqn{\Lambda} is the cumulated
hazard, \eqn{\beta} vector of parameters, \eqn{Z} explanatory variables and
\eqn{T} possibly censored survival time, are down-weighted.  Both \eqn{\Lambda}
and \eqn{\beta} are iteratively robustly estimated. Sampling weights are incorporated after iterative estimation is complete.

Numerical results are supported by a graphical tool \code{plot}, which in
a series of 5 graphs let us compare how well data are explained by the estimated
proportional hazards model with non-robust (black color) and robust method
(green color).  The first graph shows standardized difference of two estimated
survival functions; one via the Cox model and the other via Kaplan Meier
estimator.  The following four graphs show the same differences for four
strata, defined by the quartiles of the estimated linear predictor.
Comparison of estimation results along with analysis of the graphs leads
frequently to a very detailed information about the model fit (see examples).
}

\value{
An object of class \code{coxr}. See \code{\link{coxr.object}} for details.
}

\references{

Bednarski, T. (1993). Robust estimation in Cox's regression model.
Scandinavian Journal of Statistics. Vol. 20, 213--225.

Bednarski, T. (1989). On sensitivity of Cox's estimator.
Statistics and Decisions. 7, 215--228.

Grzegorek, K.(1993). On robust estimation of baseline hazard under the Cox model and via
Frechet differentiability. Preprint of the Institute of Mathematics of the Polish Academy
of Sciences.518.

Minder, C.E. & Bednarski, T. (1996). A robust method for proportional hazards
regression. Statistics in Medicine Vol. 15, 1033--1047.

Sitlani, C.M., Lumley, T, ... & Delaney, J.A.C.. (2018). Incorporating sampling weights into robust estimation of Cox proportional hazards regression model. Submitted to Statistics in Medicine.

}

\seealso{
\code{\link{coxph}}, \code{\link{coxr.object}}, \code{\link{gen_data}},
\code{\link{plot.coxr}}, \code{\link{Surv}}
}

\examples{
# Create a simple test data set using the attached function gen_data
a <- gen_data(200, c(1, 0.1, 2), cont = 0.05, p.censor = 0.30)
w1 <- rep(1,dim(a)[1])
w2 <- c(rep(0.5,200/4),rep(1,200/2),rep(2,200/4))

result1 <- coxrw(Surv(time, status) ~ X1 + X2 + X3, data = a , trunc = 0.9, weights=w1)
result1
plot(result1)

result2 <- coxrw(Surv(time, status) ~ X1 + X2 + X3, data = a , trunc = 0.9, weights=w2)
result2
plot(result2)


\keyword{robust}
\keyword{survival}
