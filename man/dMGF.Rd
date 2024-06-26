% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/dMGF.R
\name{dMGF}
\alias{dMGF}
\title{Moment generating functions (MGF) of standard normal distribution.}
\usage{
dMGF(tab, t, coef = TRUE)
}
\arguments{
\item{tab}{a dataframe with the first column contain indices of components
of a multivariate random vector \eqn{\bold{X}}, and the second column is the
order derivatives with respect to that components.}

\item{t}{vector in \eqn{\mathbb{R}^p}.}

\item{coef}{take \code{TRUE} or \code{FALSE} value to
obtain only polynomial or whole expression by multiplying the
polynomial term with the exponent term \eqn{\exp(.5 t't)}.}
}
\value{
Value of derivatives.
}
\description{
Get the polynomial term in the expression of derivatives of  moment
generating function of \eqn{N_p(0, I_p)}, with
respect to a given component and its exponent. Up to eighth order.
}
\details{
For a standard multivariate normal random variables \eqn{Y \sim N_p(0, I_p)}
\deqn{
\mathbb{E}\left(Y_1^{k_1} ... Y_p^{k_p} \exp(t'X)\right) =
\dfrac{\partial^{k_1}\dots
\partial^{k_p}}{t_1^{k_1} \dots t_p^{k_p}} \exp(t't/2) =
\mu^{(k_1)} (t_1) ... \mu^{(k_p)}(t_p) \exp(t't/2)
}
For example,
\eqn{
\mathbb{E}Y_2^4 \exp(t'Y) = \dfrac{\partial^4}{\partial t_2^4} \exp(t't/2)
= \mu^{(4)}(t_2) \exp(t't/2).
}
}
\examples{
#Calculation of above example
t <- rep(.2, 7)
tab <- data.frame(j = 2, exponent = 4)
dMGF(tab, t = t)
dMGF(tab, t = t, coef = FALSE)

}
