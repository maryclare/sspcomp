# # Pars for testing
# # rm(list = ls())
# n <- 100
# p <- 10
# X <- matrix(rnorm(n*p), nrow = n, ncol = p)
# beta <- rnorm(p)
# Xbeta <- crossprod(t(X), beta)
# y <- rbinom(n, 1, exp(Xbeta)/(1 + exp(Xbeta)))

grad <- function(X, y, Omega.inv, beta, Xbeta = crossprod(t(X), beta), j) {
  not.j <- (1:ncol(X))[!1:ncol(X) %in% j]
  if (is.matrix(Omega.inv)) {
    as.vector(crossprod(Omega.inv[j, j], beta[j]) + crossprod(Omega.inv[not.j, j], beta[not.j]) + crossprod(X[, j], ((-y) + 1/(1 + exp(-Xbeta)))))
  } else if (is.vector(Omega.inv)) {
    as.vector(beta[j]*Omega.inv[j] + crossprod(X[, j], ((-y) + 1/(1 + exp(-Xbeta)))))
  }
}

hess <- function(X, Omega.inv, beta, Xbeta = crossprod(t(X), beta), j) {
  if (is.matrix(Omega.inv)) {
    Omega.inv[j, j] + crossprod(X[, j], crossprod(diag(as.vector(exp(Xbeta)/(1 + exp(Xbeta))^2)), X[, j]))
  } else if (is.vector(Omega.inv)) {
    diag(Omega.inv[j], nrow = length(j), ncol = length(j)) + crossprod(X[, j], crossprod(diag(as.vector(exp(Xbeta)/(1 + exp(Xbeta))^2)), X[, j]))
  }
}

# obj <- function(X, y, Omega.inv, beta, Xbeta = crossprod(t(X), beta)) {
#   if (is.matrix(Omega.inv)) {
#     sum((1 - y)*Xbeta + log(1 + exp(-Xbeta))) + crossprod(beta, crossprod(Omega.inv, beta))/2
#   } else if (is.vector(Omega.inv)) {
#     sum((1 - y)*Xbeta + log(1 + exp(-Xbeta))) + sum(beta^2*Omega.inv)/2
#   }
# }

#' @export
coord.desc.logit <- function(X, y, Omega.inv, eps = 10^(-12), max.iter = 1000, print.iter = TRUE,
                             start.beta = rep(0, ncol(X)), nr = TRUE,
                             joint.beta = as.list(1:ncol(X))) {

  maxs <- colSums(abs(X))
  zero.ind <- which(maxs == 0)
  nz.ind <- which(maxs != 0)
  p <- ncol(X)
  beta <- start.beta
  Xbeta <- crossprod(t(X), beta)
  objs <- rep(NA, max.iter)
  for (i in 1:max.iter) {
    for (j in joint.beta) {
      # if (print.iter) {
      #   cat("  j = ", j, "\n")
      # }

      beta[intersect(zero.ind, j)] <- 0

      j <- intersect(nz.ind, j)
      beta.old <- beta


      # Bisection
      if (!nr & length(j) == 1 & is.vector(Omega.inv)) {
        if (Omega.inv[j] != 0) {

          aa <- -maxs[j]/Omega.inv[j]
          bb <- maxs[j]/Omega.inv[j]

          # Check endpoints
          Xbeta <- Xbeta + X[, j]*(aa - beta.old[j])
          beta[j] <- aa
          a.val <- grad(X = X, y = y, Omega.inv = Omega.inv, beta = beta, Xbeta = Xbeta, j = j)
          if (abs(a.val) < eps) {
            break;
          }
          Xbeta <- Xbeta + X[, j]*(bb - aa)
          beta[j] <- bb
          b.val <- grad(X = X, y = y, Omega.inv = Omega.inv, beta = beta, Xbeta = Xbeta, j = j)
          if (abs(b.val) < eps) {
            break;
          }
          inc <- ifelse(b.val > 0 & a.val < 0, 1, ifelse(a.val > 0 & b.val < 0, 0, NA))
          beta[j] <- (aa + bb)/2
          Xbeta <- Xbeta + X[, j]*(beta[j] - bb)
          grad.val <- grad(X = X, y = y, Omega.inv = Omega.inv, beta = beta, Xbeta = Xbeta, j = j)
          while (abs(grad.val) > 10^(-14)) {
            if (inc) {
              if (grad.val < 0) {
                aa <- beta[j]
              } else {
                bb <- beta[j]
              }
            } else if (!inc) {
              if (grad.val > 0) {
                aa <- beta[j]
              } else {
                bb <- beta[j]
              }
            }
            beta.old[j] <- beta[j]
            beta[j] <- (aa + bb)/2
            Xbeta <- Xbeta + X[, j]*(beta[j] - beta.old[j])
            grad.val <- grad(X = X, y = y, Omega.inv = Omega.inv, beta = beta, Xbeta = Xbeta, j = j)
          }
        }
        } else {

          # Newton Raphson
          diff <- Inf
          while(abs(diff) > eps) {
            # print("Grad, Hess")
            # print(grad(X = X, y = y, Omega.inv = Omega.inv, beta = beta, Xbeta = Xbeta, j = j))
            # print(hess(X = X, Omega.inv = Omega.inv, beta = beta, Xbeta = Xbeta, j = j))
            # print(beta.old[j])
            hess <- hess(X = X, Omega.inv = Omega.inv, beta = beta, Xbeta = Xbeta, j = j)
            hess.inv <- try(solve(hess), silent = TRUE)
            if (grepl("Error", hess.inv[1])) {
              hess.inv <- diag(nrow(hess))
            }
            beta[j] <- beta[j] - crossprod(hess.inv,
                                           grad(X = X, y = y, Omega.inv = Omega.inv, beta = beta, Xbeta = Xbeta, j = j))
            # print(beta[j])
            # print(diff)
            diff <- mean(abs(beta[j] - beta.old[j]))
            Xbeta <- Xbeta + crossprod(t(X[, j]), (beta[j] - beta.old[j]))
            beta.old <- beta
          }
        }
    }
    if (is.matrix(Omega.inv)) {
      objs[i] <- objMat(y = y, OmegaInv = Omega.inv, beta = beta, Xbeta = Xbeta)/length(y)
    } else {
      objs[i] <- objVec(y = y, OmegaInv = Omega.inv, beta = beta, Xbeta = Xbeta)/length(y)
    }
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

# Omega.inv <- rWishart(1, p + 2, diag(p))[, , 1]
# # Omega.inv <- Omega.inv*0
# Omega.inv <- diag(Omega.inv)
# # Omega.inv <- rep(0, p)
# joint.beta <- vector("list", length = 2)
# joint.beta[[1]] <- 1:5; joint.beta[[2]] <- 6:10
# cd.nr <- coord.desc.logit(X = X, y = y, Omega.inv = Omega.inv, joint.beta = joint.beta)
# cd.bi <- coord.desc.logit(X = X, y = y, Omega.inv = Omega.inv, nr = FALSE)
# plot(cd.nr$objs)
# plot(cd.bi$objs)
# print(cbind(cd.nr$beta,
#       cd.bi$beta, glm(y~X-1, family = binomial(link = "logit"))$coef))

