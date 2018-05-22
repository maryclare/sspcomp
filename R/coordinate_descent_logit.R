# # Pars for testing
# n <- 100
# p <- 10
# X <- matrix(rnorm(n*p), nrow = n, ncol = p)
# beta <- rnorm(p)
# Xbeta <- crossprod(t(X), beta)
# y <- rbinom(n, 1, exp(Xbeta)/(1 + exp(Xbeta)))
# Omega.inv <- diag(p)

grad <- function(X, y, Omega.inv, beta, Xbeta = crossprod(t(X), beta), j) {
  # beta[j]*Omega.inv[j, j] + sum(Omega.inv[j, -j]*beta[-j]) + sum(X[, j]*((1 - y) - 1/(1 + exp(Xbeta))))
  beta[j]*Omega.inv[j, j] + sum(Omega.inv[j, -j]*beta[-j]) + sum(X[, j]*((-y) + 1/(1 + exp(-Xbeta))))
}

hess <- function(X, Omega.inv, beta, Xbeta = crossprod(t(X), beta), j) {
  Omega.inv[j, j] + sum(X[, j]^2*exp(Xbeta)/(1 + exp(Xbeta))^2)
}

obj <- function(X, y, Omega.inv, beta, Xbeta = crossprod(t(X), beta)) {
  sum((1 - y)*Xbeta + log(1 + exp(-Xbeta))) + crossprod(beta, crossprod(Omega.inv, beta))/2
}

#' @export
coord.desc <- function(X, y, Omega.inv, eps = 10^(-12), max.iter = 1000, print.iter = TRUE,
                       start.beta = rep(0, ncol(X))) {

  p <- ncol(X)
  beta <- start.beta
  Xbeta <- crossprod(t(X), beta)
  objs <- rep(NA, max.iter)
  for (i in 1:max.iter) {
    for (j in 1:p) {
      if (print.iter) {
        # cat("  j = ", j, "\n")
      }
      if (sum(range(X[, j]) == c(0, 0)) == 2) {
        beta[j] <- 0
      } else {

        beta.old <- beta
        diff <- Inf
        while(abs(diff) > eps) {
          # print("Grad, Hess")
          # print(grad(X = X, y = y, Omega.inv = Omega.inv, beta = beta, Xbeta = Xbeta, j = j))
          # print(hess(X = X, Omega.inv = Omega.inv, beta = beta, Xbeta = Xbeta, j = j))
          # print(beta.old[j])
          beta[j] <- beta[j] - grad(X = X, y = y, Omega.inv = Omega.inv, beta = beta, Xbeta = Xbeta, j = j)/hess(X = X, Omega.inv = Omega.inv, beta = beta, Xbeta = Xbeta, j = j)
          # print(beta[j])
          # print(diff)
          diff <- abs(beta[j] - beta.old[j])
          Xbeta <- Xbeta + X[, j]*(beta[j] - beta.old[j])
          beta.old <- beta
        }
      }
    }
    objs[i] <- obj(X = X, y = y, Omega.inv = Omega.inv, beta = beta, Xbeta = Xbeta)/length(y)
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

# cd <- coord.desc(X = X, y = y, Omega.inv = Omega.inv)
# plot(cd$objs)
# cbind(cd$beta, glm(y~X-1, family = binomial(link = "logit"))$coef)
