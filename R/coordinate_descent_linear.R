#' @export
coord.desc.lin <- function(X, y, Omega.inv, sig.sq, eps = 10^(-12), max.iter = 1000, print.iter = TRUE,
                           start.beta = rep(0, ncol(X))) {

  p <- ncol(X)
  beta <- start.beta
  Xbeta <- crossprod(t(X), beta)
  objs <- rep(NA, max.iter)
  for (i in 1:max.iter) {
    for (j in 1:p) {
      if (print.iter) {
        cat("  j = ", j, "\n")
      }
      if (sum(range(X[, j]) == c(0, 0)) == 2) {
        beta[j] <- 0
      } else {
        Xbeta.noj <- Xbeta - X[, j]*beta[j]
        beta[j] <- (Omega.inv[j, j] + sum(X[, j]^2)/sig.sq)^(-1)*((sum(X[, j]*y) - sum(X[, j]*Xbeta.noj))/sig.sq - sum(Omega.inv[j, -j]*beta[-j]))
        Xbeta <- Xbeta.noj + X[, j]*beta[j]
      }
    }
    objs[i] <- (sum((y - Xbeta)^2)/sig.sq + crossprod(beta, crossprod(Omega.inv, beta)))/length(y)
    if (print.iter) {
      cat("  i = ", i, "\n")
      cat("obj = ", objs[i], "\n")
    }
    if (i > 1) {
      if (abs(objs[i] - objs[i - 1]) < eps) {
        break
      }
    }
  }
  return(list("beta" = beta, "objs" = objs[1:i]))
}

## Test code
# n <- 10
# p <- 5
# X <- matrix(rnorm(p*n), nrow = n, ncol = p)
# beta <- rnorm(p)
# y <- X%*%beta + rnorm(n)
# Omega.inv <- (1 - 0.25)*diag(p) + 0.25*matrix(1, nrow = p, ncol = p)
# sig.sq <- 2
#
# lm(y~X-1)$coef
# coord.desc.lin(X = X, y = y, sig.sq = sig.sq, Omega.inv = 0*Omega.inv, print.iter = FALSE)$beta
#
# solve(crossprod(X)/sig.sq + Omega.inv)%*%crossprod(X, y)/sig.sq
# coord.desc.lin(X = X, y = y, sig.sq = sig.sq, Omega.inv = Omega.inv, print.iter = FALSE)$beta
