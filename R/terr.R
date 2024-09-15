#' The terr function is to generate data based on different distributions.
#'
#' @param n is the total number of observations
#' @param p is the number of variables
#' @param K is the number of subsets
#' @param nr is the number of outliers
#' @param alpha is the significance level
#' @param sigma1 is a parameter related to the distribution
#' @param sigma2 is a parameter related to the distribution
#' @param beta is a coefficient vector
#' @param dist is the type of distribution
#'
#' @return A list containing the design matrix X and the response vector y
#' @export
#'
#' @importFrom stats qt rnorm rt
#' @examples
#' # Example usage of the terr function
#' n <- 100
#' p <- 5
#' K <- 3
#' nr <- 10
#' alpha <- 0.05
#' sigma1 <- 2
#' sigma2 <- 5
#' beta <- c(1, 2, 3, 4, 5)
#' dist <- "Student t Distribution"
#' result <- terr(n, p, K, nr, alpha, sigma1, sigma2, beta, dist)
#' X <- result$X
#' y <- result$y

terr <- function(n, p, K, nr, alpha, sigma1, sigma2, beta, dist) {
  if (dist == "Student t Distribution") {
    data1 <- matrix(rnorm((n - nr) * p, 0, sigma1), nrow = n - nr, ncol = p)
    data2 <- matrix(rt(nr * p, df = 10), nrow = nr, ncol = p)
  } else if (dist == "Student t Distribution: Univariate with loc and scal") {
    data1 <- matrix(rnorm((n - nr) * p, 0, sigma1), nrow = n - nr, ncol = p)
    data2 <- matrix(rt(nr * p, df = 10), nrow = nr, ncol = p)
  } else if (dist == "Student t Distribution: Precision Parameterization") {
    data1 <- matrix(rnorm((n - nr) * p, 0, sigma1), nrow = n - nr, ncol = p)
    data2 <- matrix(rt(nr * p, df = 10), nrow = nr, ncol = p)
  } else if (dist == "Skew-t Distribution") {
    data1 <- matrix(rnorm((n - nr) * p, 0, sigma1), nrow = n - nr, ncol = p)
    data2 <- matrix(rt(nr * p, df = 10), nrow = nr, ncol = p)
  } else if (dist == "Generalized Hyperbolic Student-t distribution") {
    data1 <- matrix(rnorm((n - nr) * p, 0, sigma1), nrow = n - nr, ncol = p)
    data2 <- matrix(rt(nr * p, df = 10), nrow = nr, ncol = p)
  } else if (dist == "Standardized generalized hyperbolic Student-t Distribution") {
    data1 <- matrix(rnorm((n - nr) * p, 0, sigma1), nrow = n - nr, ncol = p)
    data2 <- matrix(rt(nr * p, df = 10), nrow = nr, ncol = p)
  }

  X <- rbind(data1, data2)

  e1 <- rnorm(n - nr, 0, sigma1)
  y1 <- X[1:(n - nr), ] %*% beta + e1
  e2 <- rnorm(nr, 0, sigma2)
  y2 <- X[(n - nr + 1):n, ] %*% beta + e2

  y <- c(y1, y2)

  return(list(X = X, y = y))
}
