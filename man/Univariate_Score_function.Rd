% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/cox.R, R/score_plot.R
\name{Univariate_Score_function}
\alias{Univariate_Score_function}
\alias{cox}
\alias{score_plot1D}
\title{Graphical plots to assess the univarite noramality assumption of data.}
\usage{
cox(x, P = NULL, lambda = 0.5, x.dist = NULL)

score_plot1D(x, P = NULL, lambda = 0.5, x.dist = NULL, ori.index = NULL)
}
\arguments{
\item{x}{univariate data.}

\item{P}{vector of weight.}

\item{lambda}{smoothing parameter, default is \eqn{0.5}.}

\item{x.dist}{the minimum distance between two data points in vector x.}

\item{ori.index}{original index of vector x, default is \code{NULL}
when index is just the order.}
}
\value{
\code{cox} returns the estimate of score function.
\itemize{
\item{\code{x}: The updated univariate data if merging happens.}
\item{\code{a}: Score value estimated at \code{x}.}
\item{\code{P}: Updated weight (if merging happens).}
\item{\code{slt}: Index of merged data point
(is \code{NULL} if \code{x.dist = NULL}).}
}

\code{score_plot1D} returns score functions together with
2-sigma bands for outlier detection.
\itemize{
\item{\code{plot}: plot of estimate score function and its band.}
\item{\code{outlier}: index of outliers.}
}
}
\description{
Score function of a univariate normal distribution is
a straight line. A non-linear graph of score function estimator shows
evidence of non-normality.

Outliers are detected using the 2-sigma bands method.
}
\details{
To avoid the singularity of coefficient matrices in spline method, points
with distance less than \code{x.dist} are merged and weight of the
representative points is updated by the summation of weight of
discarded points.

Under null hypothesis, a unbiased estimator score function of a
given data point \eqn{x_k} is
\deqn{
\hat{\psi}(x_k) = \dfrac{n - 4}{n - 2} \dfrac{x_k - \bar{X}_{-k}}{S_{-k}^2}
}
and if \eqn{a_{k}} is the estimate score from function \code{cox} at
the point \eqn{x_k}, then
\deqn{a_k\in \hat{\psi}(x_k) \pm 2 \sqrt{\hat{\text{Var}}(\hat{\psi}(x_k))}.}
Hence points outside the 2-sigma bands are outliers.
}
\examples{
set.seed(1)
x <- rnorm(100, 2, 4)
re <- cox(sort(x))
plot(re$x, re$a, xlab = "x", ylab = "Estimated Score",
 main = "Estimator of score function")
abline(0, 1)

set.seed(1)
x <- rnorm(100, 2, 4)
score_plot1D(sort(x))

}
\references{
\insertRef{ref:sp_estimate}{PlotNormTest}
}
