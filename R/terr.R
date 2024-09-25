#' Generate Data with T-distributed Errors
#'
#' @param n Number of observations.
#' @param nr Number of observations with different error distribution.
#' @param p Number of predictors.
#' @param dist_type Type of distribution for the error terms.
#' @param ... Additional parameters for specific distributions.
#'
#' @return A list containing the design matrix X, the response vector Y, and the error vector e.
#' @export
#'
#' @importFrom stats runif
#' @importFrom stats rt
#' @importFrom LaplacesDemon rst rstp
#' @importFrom fBasics rght rsght
#'
#' @examples
#' set.seed(12)
#' n <- 1200
#' nr <- 200
#' p <- 5
#' data <- terr(n, nr, p, dist_type = "student_t")
#' print(data$X)
#' print(data$Y)
#' print(data$e)
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
