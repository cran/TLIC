#' The TLIC function is to perform a specific analysis based on the input data.
#'
#' @param X is a design matrix
#' @param y is a random response vector of observed values
#' @param alpha is the significance level
#' @param K is the number of subsets
#'
#' @return A list containing various results of the analysis
#' @export
#'
#' @importFrom stats qt rnorm rt
#' @examples
#' # Example usage of the TLIC function
#' set.seed(12)
#' X <- matrix(data = sample(1:3, 1200 * 5, replace = TRUE), nrow = 1200, ncol = 5)
#' b <- sample(1:3, 5, replace = TRUE)
#' e <- rnorm(1200, 0, 1)
#' Y <- X %*% b + e
#' alpha <- 0.05
#' K <- 10
#' result <- TLIC(X, Y, alpha, K)
#' MUopt <- result$MUopt
#' Bopt <- result$Bopt
#' MAEMUopt <- result$MAEMUopt
#' MSEMUopt <- result$MSEMUopt
#' opt <- result$opt
#' Yopt <- result$Yopt

TLIC <- function(X, y, alpha, K) {
  n <- nrow(X)
  p <- ncol(X)
  N <- L1 <- c(1:K)
  Rm <- matrix(rep(0, n * K), ncol = K)
  mr <- matrix(rep(0, K * n), ncol = n)
  for (i in 1:K) {
    mr[i, ] <- sample(1:n, n, replace = T)
    r <- matrix(c(1:n, mr[i, ]), ncol = n, byrow = T)
    Rm[, i] <- r[2,]
    R <- matrix(rep(0, n * n), ncol = n)
    R[t(r)] <- 1
    X1 <- R %*% X
    Y1 <- R %*% y
    Hr <- X1 %*% solve(crossprod(X1)) %*% t(X1)
    I1 <- diag(rep(1, n))
    SX <- (t(Y1) %*% (I1 - Hr) %*% Y1) / (n - p)
    SY <- sqrt(t(Y1) %*% (I1 - Hr) %*% Y1) / (n - p)
    C1 <- sum(diag(X1 %*% solve(crossprod(X1)) %*% t(X1))) / n
    L1[i] <- 2 * SY * C1 * qt(1 - alpha / 2, n - p)
    N[i] <- det(t(X1) %*% X1)
  }
  opt1 <- Rm[, which.min(L1)]
  opt2 <- Rm[, which.max(N)]
  opt <- intersect(opt1, opt2)
  Yopt <- y[opt]
  Xopt <- X[opt,]
  Bopt <- solve(crossprod(Xopt)) %*% t(Xopt) %*% Yopt
  MUopt <- Xopt %*% Bopt
  Nopt <- length(Yopt)
  E5 <- (t(Yopt - MUopt) %*% (Yopt - MUopt)) / Nopt
  A5 <- sum(abs(Yopt - MUopt)) / Nopt
  return(list(MUopt = MUopt, Bopt = Bopt, MAEMUopt = A5, MSEMUopt = E5, opt = opt, Yopt = Yopt))
}

