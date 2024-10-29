#' terr function is used to generate a dataset where the error term follows a T-distribution
#'
#' This terr function generates a dataset with a specified number of observations and predictors, along with a response vector that has an error term following a T-distribution.
#' @param n is the number of observations
#' @param nr is the number of observations with a different error T distribution
#' @param p is the dimension of the observation
#' @param dist_type is the type where the error term obeys a T-distribution
#' @param ... is additional arguments for the T-distribution function
#'
#' @return X,Y,e
#' @export
#'

#' @examples
#' set.seed(12)
#' data <- terr(n = 1200, nr = 200, p = 5, dist_type = "student_t")
#' str(data)

terr <- function(n, nr, p, dist_type, ...) {
  beta <- sort(runif(p, 1, 5))
  X <- matrix(runif(n * p, 0, 1), ncol = p)

  if (dist_type == "student_t") {
    x2 <- runif(n, 0, 1)
    df <- exp(exp(0.5 - x2))
    e1 <- rt(n - nr, df)
    e2 <- rt(nr, df)
  } else if (dist_type == "student_t_loc_scale") {
    e1 <- rst(n - nr, 0, 1, 10, ...)
    e2 <- rst(nr, 0, 5, 7, ...)
  } else if (dist_type == "student_t_precision") {
    e1 <- rstp(n - nr, 0, 1, 10, ...)
    e2 <- rstp(nr, 0, 5, 7, ...)
  } else if (dist_type == "skew_t") {
    xi <- 5
    omega <- 1
    alpha <- -5
    beta <- 8
    e1 <- rst(n - nr, xi, omega, alpha)
    alpha <- -2
    beta <- 4
    e2 <- rst(nr, xi, omega, alpha)
  } else if (dist_type == "gen_hyperbolic_stud_t") {
    beta <- 0.1
    delta <- 1
    mu <- 0
    nu <- 10
    e1 <- rght(n - nr, beta = beta, delta = delta, mu = mu, nu = nu, ...)
    e2 <- rght(nr, beta = beta, delta = delta, mu = mu, nu = nu, ...)
  } else if (dist_type == "std_gen_hyperbolic_stud_t") {
    beta <- 0.1
    delta <- 1
    mu <- 0
    nu <- 10
    e1 <- rsght(n - nr, beta = beta, delta = delta, mu = mu, nu = nu, ...)
    e2 <- rsght(nr, beta = beta, delta = delta, mu = mu, nu = nu, ...)
  } else {
    stop("Unknown distribution type: ", dist_type)
  }

  e <- c(e1, e2)
  Y <- X %*% beta + e

  return(list(X = X, Y = Y, e = e))
}
