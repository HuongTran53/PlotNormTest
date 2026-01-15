# -- Filename: score_plot.R --#
#' @rdname Univariate_Score_function
#' @description Outliers are detected using the 2-sigma bands method.
#' @param ori.index original index of vector x, default is \code{NULL}
#' when index is just the order.
#' @return \code{score_plot1D} returns score functions together with
#' 2-sigma bands for outlier detection.
#' \itemize{
#' \item{\code{plot}: plot of estimate score function and its band.}
#' \item{\code{outlier}: index of outliers.}
#' }
#' @details Under null hypothesis, a unbiased estimator score function of a
#'  given data point \eqn{x_k} is
#' \deqn{
#' \hat{\psi}(x_k) = \dfrac{n - 4}{n - 2} \dfrac{x_k - \bar{X}_{-k}}{S_{-k}^2}
#' }
#' and if \eqn{a_{k}} is the estimate score from function \code{cox} at
#' the point \eqn{x_k}, then
#' \deqn{a_k\in \hat{\psi}(x_k) \pm 2 \sqrt{\hat{\text{Var}}(\hat{\psi}(x_k))}.}
#' Hence points outside the 2-sigma bands are outliers.
#' @import ggplot2
#' @importFrom stats weighted.mean
#' @importFrom rlang .data
#' @export
#' @examples
#' set.seed(1)
#' x <- rnorm(100, 2, 4)
#' score_plot1D(sort(x))
#'
score_plot1D <- function(x, P = NULL,
                         lambda = NULL, x.dist = NULL,ori.index = NULL){
  nx <- length(x)
  if (is.unsorted(x) == T) {
    stop("x should be sort")
  }
  if (is.null(ori.index)) {
    ori.index <- 1:nx
  } else {
    if (length(ori.index) != nx) {
      stop("Uncomformable length of index and vector")
    }
  }
  if (is.null(P)) {
    P <- rep(1/nx, nx)
  }
  else {
    if (length(P) < nx) {
      stop("Incorrect dimensions")
    }
  }
  index <- ori.index
  recox <- cox(x, P = P, lambda = lambda, x.dist = x.dist)
  ######################3
  ## Adjust to make a having balance tails
  a_raw <- as.numeric(recox$a)
  w <- as.numeric(recox$P)
  w <- w / sum(w)

  a_bar <- sum(w * a_raw)
  var_a <- sum(w * (a_raw - a_bar)^2) / (1 - sum(w^2))  # unbiased weighted var
  a <- (a_raw - a_bar) / sqrt(var_a)
  #####################3##
  a <- as.numeric(unlist(recox$a))
  x <- as.numeric(unlist(recox$x))
  weight <- as.numeric(unlist(recox$P))
  slt <- as.numeric(unlist(recox$slt))
  if (length(slt) != 0) {
    index <- index[-slt]
  }
  nx <- length(x)
  f <- c()
  upper <- c()
  lower <- c()
  for (i in 1:nx) {
    x.temp <- x[-i]
    weight.temp <- weight[-i]
    # xbar <- stats::weighted.mean(x.temp, weight.temp)
    xbar <- sum(weight.temp*x.temp)
    ## Added for weighted variance
    varx <- sum(weight.temp * (x.temp - xbar)^2)/(1 - sum(weight.temp^2))
    # varx1 <- sum(weight.temp * (x.temp - xbar)^2)
    fi <- ((nx - 4)/(nx - 2)) * ((x[i] - xbar)/varx)
    bi <- 2 * (nx - 4)/((nx - 2)^2) * fi^2 + (nx - 4)/((nx -
                                                          1) * (nx - 2) * varx)
    upperi <- fi + 2 * sqrt(bi)
    loweri <- fi - 2 * sqrt(bi)
    f <- c(f, fi)
    upper <- c(upper, upperi)
    lower <- c(lower, loweri)
  }
  df <- data.frame(x = x, a = a, f = f, upper = upper, lower = lower)
  ind <- which(df$a <= df$lower | df$a >= df$upper)
  # 20250113 - HT - update outliers based on symmetric distance
  # df$dev <- df$a - df$f
  # df$crit <- 2*sqrt(bi)
  # df$Type <- abs(df$dev) > df$crit
  # ind <- which(abs(df$dev) > df$crit)
  out <- index[ind]
  #####
  df$Type <- (df$a <= df$lower | df$a >= df$upper)
  # df$Type <- as.factor(ifelse(df$a < lower | df$a > upper,
  #                             "Outlier", "Non outlier"))
  df$Type <- factor(df$Type,
                    levels = c(T, F),
                    labels = c("Outlier", "Non-outlier"))
  g <- ggplot2::ggplot(data = df) +
    geom_line(mapping = aes(x, f),
              colour = "darkorange3", linewidth = 0.7) +
    ggplot2::geom_line(mapping = aes(x, upper), linetype = "twodash", linewidth = 0.7, colour = "darkorange3") +
    ggplot2::geom_line(mapping = aes(x, lower), linetype = "twodash",
                       linewidth = 0.7, colour = "darkorange3") +
    ggplot2::geom_point(aes(x, a, color = .data$Type, shape = .data$Type, size = .data$Type)) +
    ggplot2::scale_shape_manual(values = c("Non-outlier" = 1,"Outlier" =  17)) +
    ggplot2::scale_color_manual(values = c("Non-outlier" = "steelblue3",
                                           "Outlier" = "seagreen")) +
    ggplot2::scale_size_manual(values = c("Non-outlier" =  1.5, "Outlier" =2)) +
    ggplot2::labs(y = expression(hat(psi))) +
    ggplot2::coord_fixed() +
    ggplot2::theme_bw() + ggplot2::theme(aspect.ratio = 1/1,
                                         panel.grid = element_blank(), axis.line = element_line(colour = "black"),
                                         axis.text = element_text(size = 12), axis.title = element_text(size = 14,
                                                                                                        face = "bold"), legend.background = element_rect(size = 0.5,
                                                                                                                                                         linetype = "solid"), legend.text = element_text(size = 12)) +
    ylab("Score function") + ggplot2::theme(legend.position = "top")
  return(list(plot = g, outlier = unique(out)))
}

