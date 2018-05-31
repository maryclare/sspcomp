# # # Pars for testing
# n <- 20
# p <- 5
# X <- matrix(rnorm(n*p), nrow = n, ncol = p)
# beta <- rnorm(p)
# Xbeta <- crossprod(t(X), beta)
# y <- rbinom(n, 1, exp(Xbeta)/(1 + exp(Xbeta)))
# Omega.inv <- diag(p)

#' @export
coord.desc.logit <- function(X, y, Omega.inv, eps = 10^(-12), max.iter = 1000, print.iter = TRUE,
                       start.beta = rep(0, ncol(X)), eps.inner = 10^(-12)) {

  if (is.vector(Omega.inv)) {
    beta <- coordDescLogitVec(y = y, X = X, OmegaInv = Omega.inv,
                              eps = eps, maxit = max.iter, betaStart = start.beta, epsInner = eps.inner)
  } else if (is.matrix(Omega.inv)) {
    beta <- coordDescLogitMat(y = y, X = X, OmegaInv = Omega.inv,
                              eps = eps, maxit = max.iter, betaStart = start.beta, epsInner = eps.inner)
  }


  return(list("beta" = beta))
}

# Omega.inv <- rWishart(1, p + 2, diag(p))[, , 1] # diag(abs(rnorm(ncol(X))))
# # Omega.inv <- diag(Omega.inv)
# is.vector(Omega.inv)
# system.time(cd <- coord.desc.logit(X = X, y = y, Omega.inv = Omega.inv, print.iter = FALSE))
# if (is.vector(Omega.inv)) {
#   system.time(cd.cp <- coordDescLogitVec(y = y, X = X, OmegaInv = Omega.inv, eps = 10^(-14), maxit = 1000, betaStart = rep(0, p)))
# } else if (is.matrix(Omega.inv)) {
#   system.time(cd.cp <- coordDescLogitMat(y = y, X = X, OmegaInv = Omega.inv, eps = 10^(-14), maxit = 1000, betaStart = rep(0, p)))
# }
# print(cbind(cd.cp,
#             cd$beta, glm(y~X-1, family = binomial(link = "logit"))$coef))

