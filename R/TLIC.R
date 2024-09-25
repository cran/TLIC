#' TLIC function based on LIC with T-distributed errors
#'
#' @param X is a design matrix
#' @param Y is a random response vector of observed values
#' @param alpha is the significance level
#' @param K is the number of subsets
#' @param nk is the sample size of subsets
#' @param dist_type Type of distribution for the error terms.
#'
#' @return MUopt, Bopt, MAEMUopt, MSEMUopt, opt, Yopt
#' @export
#'
#' @examples
#' set.seed(12)
#' n <- 1200
#' nr <- 200
#' p <- 5
#' data <- terr(n, nr, p, dist_type = "student_t")
#' TLIC(data$X, data$Y, alpha = 0.05, K = 10, nk = n / 10, dist_type = "student_t")
#'
#' @importFrom stats qt
#' @importFrom stats runif
#' @importFrom stats rt
#' @importFrom LaplacesDemon rst
#' @importFrom LaplacesDemon rstp
#' @importFrom fBasics rght
#' @importFrom fBasics rsght
TLIC <- function(X, Y, alpha = 0.05, K = 10, nk = NULL, dist_type = "student_t") {
  if (is.null(nk)) {
    nk <- nrow(X) / K
  }

  n <- nrow(X)
  p <- ncol(X)

  N <- L1 <- c(1:K)
  Rm <- matrix(rep(0, nk * K), ncol = K)
  mr <- matrix(rep(0, K * nk), ncol = nk)

  for (i in 1:K) {
    mr[i, ] <- sample(1:n, nk, replace = FALSE)
    r <- matrix(c(1:nk, mr[i, ]), ncol = nk, byrow = TRUE)
    Rm[, i] <- r[2, ]
    R <- matrix(rep(0, nk * n), ncol = n)
    R[t(r)] <- 1

    X1 <- R %*% X
    Y1 <- R %*% Y

    Hr <- X1 %*% solve(crossprod(X1)) %*% t(X1)
    I1 <- diag(rep(1, nk))

    SX <- (t(Y1) %*% (I1 - Hr) %*% Y1) / (nk - p)
    SY <- sqrt(t(Y1) %*% (I1 - Hr) %*% Y1) / (nk - p)
    C1 <- sum(diag(X1 %*% solve(crossprod(X1)) %*% t(X1))) / nk
    L1[i] <- 2 * SY * C1 * qt(1 - alpha / 2, nk - p)
    N[i] <- det(t(X1) %*% X1)
  }

  opt1 <- Rm[, which.min(L1)]
  opt2 <- Rm[, which.max(N)]
  opt <- intersect(opt1, opt2)
  Yopt <- Y[opt]
  Xopt <- X[opt, ]
  Bopt <- solve(crossprod(Xopt)) %*% t(Xopt) %*% Yopt
  MUopt <- Xopt %*% Bopt
  Nopt <- length(Yopt)
  E5 <- (t(Yopt - MUopt) %*% (Yopt - MUopt)) / Nopt
  A5 <- sum(abs(Yopt - MUopt)) / Nopt

  return(list(MUopt = MUopt, Bopt = Bopt, MAEMUopt = A5, MSEMUopt = E5, opt = opt, Yopt = Yopt))
}
