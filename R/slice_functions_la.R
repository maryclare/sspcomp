#' @export
samp.Omega.inv <- function(Beta, pr.V.inv = diag(ncol(Beta)),
                           pr.df = ncol(Beta) + 2, str = "uns") {
  p <- ncol(Beta)
  if (str == "uns") {
    V.inv <- crossprod(Beta) + pr.V.inv
    df <- nrow(Beta) + pr.df
    return(rWishart(1, df, solve(V.inv))[, , 1])
  } else if (str == "het") {
    b <- apply(Beta, 2, function(x) {sum(x^2)})/(2) + diag(pr.V.inv)/2
    a <- rep(nrow(Beta), ncol(Beta))/2 + pr.df/2
    return(diag(rgamma(p, shape = a, rate = b)))
  } else if (str == "con") {
    b <- sum(apply(Beta, 2, function(x) {sum(x^2)}))/(2) + sum(diag(pr.V.inv))/2
    # I'm a little worried about the code below if 'Beta' is a matrix wtih more than 1 column,
    # Should check. I think it works!
    a <- sum(rep(nrow(Beta), ncol(Beta)))/2 + p*pr.df/2
    return(rgamma(1, shape = a, rate = b)*diag(p))
  }
}

#' @export
samp.beta <- function(XtX, Xty, Omega.inv, sig.sq) {

  V.inv <- XtX/sig.sq + Omega.inv
  V.inv.eig <- eigen(V.inv/2 + Matrix::t(V.inv)/2)
  V.rt <- Matrix::tcrossprod(Matrix::tcrossprod(V.inv.eig$vectors,
                                                diag(sqrt(ifelse(1/V.inv.eig$values > 0, 1/V.inv.eig$values, 0)), nrow = length(V.inv.eig$values), ncol = length(V.inv.eig$values))),
                             V.inv.eig$vectors)
  V <- Matrix::tcrossprod(Matrix::tcrossprod(V.inv.eig$vectors,
                                             diag(ifelse(1/V.inv.eig$values > 0, 1/V.inv.eig$values, 0), nrow = length(V.inv.eig$values), ncol = length(V.inv.eig$values))),
                          V.inv.eig$vectors)
  m <- Matrix::crossprod(V, Xty/sig.sq)
  return(m + Matrix::crossprod(V.rt, rnorm(length(Xty))))
}


#' @export
make.ar.mat <- function(p, rho, inv) {
  if (inv) {
    ARMat <- diag(p)
    for (l in 1:p) {
      if (l == 1 | l == p) {
        ARMat[l, l] <- (1 - rho^2)^(-1)
      } else {
        ARMat[l, l] <- (1 + rho^2)*(1 - rho^2)^(-1)
      }
      if (l < p) {
        ARMat[l, l + 1] <- -rho*(1 - rho^2)^(-1)
        ARMat[l + 1, l] <- -rho*(1 - rho^2)^(-1)
      }
    }
  } else {
    ARMat <- diag(p)
    for (i in 1:p) {
      for (j in 1:p) {
        ARMat[i, j] <- rho^(abs(i - j))
      }
    }
  }

  return(ARMat)
}

# Snagged from tensr
mat <- function (A, k) {
  Ak <- t(apply(A, k, "c"))
  if (nrow(Ak) != dim(A)[k]) {
    Ak <- t(Ak)
  }
  return(Ak)
}

f.deltas <- function(deltas, c) {
  alpha <- c/2
  sin(alpha*deltas)^(alpha/(1 - alpha))*sin((1 - alpha)*deltas)/sin(deltas)^(1/(1 - alpha))
}

g.delta <- function(deltas, xi, c) {

  f.val <- f.deltas(deltas = deltas, c = c)
  log(f.val) - xi*f.val
}

# Performs univariate slice sampling for an arbitrary RV, requires
slice <- function(x.tilde, ll.fun, var.lim,
                  ll.args) {
  z <- do.call(ll.fun, c(x.tilde, ll.args)) - rexp(1)
  a <- var.lim[1]
  b <- var.lim[2]
  d <- runif(1, a, b)
  ll.d <- do.call(ll.fun, c(d, ll.args))
  while (z > ll.d) {

    if (d < x.tilde) {
      a <- d
    } else {
      b <- d
    }
    d <- runif(1, a, b)
    ll.d <- do.call(ll.fun, c(d, ll.args))
    # If we end up in a bad spot, reinitalize interval choice
    if (is.nan(ll.d)) {
      a <- var.lim[1]
      b <- var.lim[2]
      d <- runif(1, a, b)
      ll.d <- do.call(ll.fun, c(d, ll.args))
    }

  }
  x.tilde.p <- d
  return(x.tilde.p)
}

#' @export
cond.rho.log <- function(rho, B, pr.a, pr.b, j) {
  p <- dim(B)
  B.mat <- mat(B, j)

  c1 <- -(prod(p[-j])*(p[j] - 1)/2)*log(1 - rho^2)

  O.i <- make.ar.mat(p = p[j], rho = rho, inv = TRUE)

  c2 <- -sum(diag(crossprod(B.mat, crossprod(O.i, B.mat))))/2

  c3 <- dbeta((rho + 1)/2, pr.a, pr.b, log = TRUE)

  # cat("rho=", rho, "\n")
  # cat("c1=", c1, "\n")
  # cat("c2=", c2, "\n")
  # cat("c3=", c3, "\n")
  return(c1 + c2 + c3)

}

cond.tilde.c.log <- function(S, tilde.c, pr.shape, pr.rate) {

  c1 <- (length(c(S))*exp(tilde.c) + 1 + (pr.shape - 1))*tilde.c
  if (exp(tilde.c) < 10^(-307)) {
    c2 <- -length(c(S))*log(Inf)
  } else if (exp(tilde.c) > 171) {
    c2 <- -length(c(S))*log(Inf)
  } else {
    c2 <- -length(c(S))*log(gamma(exp(tilde.c)))
  }
  c3 <- (exp(tilde.c) - 1)*sum(log(c(S^2)))
  c4 <- -(pr.rate + sum(c(S^2)))*exp(tilde.c)

  return(c1 + c2 + c3 + c4)

}

sample.c <- function(S, c.old, tune) {

  tilde.c.old <- log(c.old)
  tilde.c.new <- tilde.c.old + tune*rnorm(1)

  llik.old <- cond.tilde.c.log(S = S, tilde.c = tilde.c.old, pr.shape = 0, pr.rate = 0)
  llik.new <- cond.tilde.c.log(S = S, tilde.c = tilde.c.new, pr.shape = 0, pr.rate = 0)
  # if (exp(tilde.c.new) > 100 | exp(tilde.c.new) < 10^(-14)) {
  #   llik.new <- -Inf
  # }

  acc <- 0
  u <- runif(1)
  if (u < exp(llik.new - llik.old)) {
    c.old <- exp(tilde.c.new)
    acc <- 1
  }
  return(list("c" = c.old, "acc" = acc))

}


V.factor <- function(X, U, max.it = 1) { # This DOES NOT WORK WELL when we use more iterations
  n <- dim(X)[1]
  p <- dim(X)[-1]
  q <- ncol(U)
  V <- V.inv <- V.half <- lapply(c(q, p, 1), diag)
  V.inv[[1]] <- crossprod(U)
  V.half[[1]] <- sym.sq.root.inv(V.inv[[1]])
  V[[1]] <- crossprod(V.half[[1]])

  if (n < prod(p)) {
    svd <- svd(mat(X, 1))
    R <- svd$v
    d <- c(svd$d)
  } else {
    svd <- svd(crossprod(mat(X, 1)))
    R <- svd$v
    d <- sqrt(c(svd$d))
  }

  X.tilde <- tcrossprod(diag(sqrt(1 + d^2)), R)
  V.half[[length(V.half)]] <- sqrt(rep(1, prod(p)) - apply(R, 1, function(x) {sum(x^2*(d^2/(1 + d^2)))}))
  V.inv[[length(V.inv)]] <- sqrt(diag(crossprod(X.tilde)))
  X.tilde.arr <- array(c(tcrossprod(X.tilde, diag(V.half[[length(V.half)]]))), c(nrow(X.tilde), p))

  for (m in 1:max.it) {
    for (k in 1:length(p)) {

      arr.1 <- X.tilde.arr
      arr.2 <- X.tilde.arr
      for (l in (1:length(p))[-k]) {
        arr.1 <- amprod.mc(arr.1, rep(1, p[l])%*%V[[l + 1]], l + 1)
        arr.2 <- amprod.mc(arr.2, t(rep(1, p[l])), l + 1)
      }
      V.mat <- crossprod(mat(arr.1, 1), mat(arr.2, 1))
      if (sum(V.mat != t(V.mat)) > 0) {
        V.mat <- (V.mat + t(V.mat))/2
      }
      V.mat.val <- eigen(V.mat)$values
      if (min(V.mat.val) < 0) {
        V.mat <- V.mat + (abs(min(V.mat.val)) + 10^(-8))*diag(p[k])
      }
      V[[k + 1]] <- cov2cor(solve(V.mat))
      V.inv[[k + 1]] <- solve(V[[k + 1]])
      V.half[[k + 1]] <- sym.sq.root(V[[k + 1]])
    }
  }
  return(list("V.half" = V.half, "V.inv" = V.inv))
}

# Get Kronecker indices
get.kron.row <- function(i, Omega) {
  p <- unlist(lapply(Omega, function(x) {nrow(x)}))
  r <- numeric(length(p))
  # Lowest index moves fastest!
  r[1] <- ifelse(i%%p[1] == 0, p[1], i%%p[1])
  for (k in 2:length(p)) {
    r[k] <- floor((i - 1)/prod(p[1:(k - 1)])) + 1
    r[k] <- ifelse(r[k]%%p[k] == 0, p[k], r[k]%%p[k])
  }
  # Can improve this
  row <- rep(Omega[[1]][r[1], ], times = prod(p[-1]))
  for (k in 2:length(p)) {
    row <- rep(Omega[[k]][r[k], ], each = prod(p[1:(k - 1)]))*row
  }
  return(row)
}

#' @export
sym.sq.root <- function(A) {
  A.eig <- eigen((A + t(A))/2)
  crossprod(t(A.eig$vectors), tcrossprod(diag(sqrt(ifelse(A.eig$values > 0, A.eig$values, 0)),
                                              nrow = nrow(A), ncol = ncol(A)), A.eig$vectors))
}

#' @export
sym.sq.root.inv <- function(A) {
  A.eig <- eigen(A)
  crossprod(t(A.eig$vectors), tcrossprod(diag(sqrt(ifelse(A.eig$values > 0, 1/A.eig$values, 0)),
                                              nrow = nrow(A), ncol = ncol(A)), A.eig$vectors))
}

#' @export
ei.inv <- function(A) {
  A.eig <- eigen(A)
  crossprod(t(A.eig$vectors), tcrossprod(diag(ifelse(A.eig$values > 0, 1/A.eig$values, 0),
                                              nrow = nrow(A), ncol = ncol(A)), A.eig$vectors))
}

#' @export
amprod.mc <- function (A, M, k) {
  K <- length(dim(A))
  AM <- crossprod(t(M), mat(A, k))
  AMA <- array(AM, dim = c(dim(M)[1], dim(A)[-k]))
  return(aperm(AMA, match(1:K, c(k, (1:K)[-k]))))
}

#' @export
atrans.mc <- function(A, B) {
  X <- A
  for (k in 1:length(B)) {
    X <- amprod.mc(X, B[[k]], k)
  }
  return(X)
}

h.log <- function(theta, d0, d1, X, U, y, beta, Omega.half, beta.tilde, V.inv, Omega.inv, sig.sq, reg) {
  beta = beta.tilde + d0*sin(theta) + d1*cos(theta)
  n <- length(y)
  p <- unlist(lapply(Omega.inv, function(x) {nrow(x)}))
  q <- ncol(U)
  gamma <- beta[1:q]
  B <- array(beta[(q + 1):length(beta)], dim = p)
  Xbeta <- crossprod(t(X), c(B))
  Ugamma <- crossprod(t(U), gamma)

  if (reg == "logit") {
    enXbeUga <- exp(-Xbeta-Ugamma)
    enXbeUga[is.infinite(enXbeUga)] <- 10^(308)
    c1 <- sum(-(1 - y)*(Xbeta+Ugamma) - log(1 + enXbeUga))
  } else if (reg == "linear") {
    c1 <- -sum((y - Ugamma - Xbeta)^2/sig.sq)/2
  }

  c2 <- -sum(c(atrans.mc(B, Omega.inv))*c(B))/2

  if (is.list(V.inv)) {

    q <- nrow(V.inv[[1]])
    p <- unlist(lapply(V.inv[2:(length(V.inv) - 1)], nrow))
    gamma.diff <- (beta[1:q] - beta.tilde[1:q])
    beta.diff <- (beta[(q + 1):length(beta)] - beta.tilde[(q + 1):length(beta.tilde)])
    c3 <- sum(crossprod(gamma.diff, crossprod(V.inv[[1]], gamma.diff)))/2 +
      sum(c(atrans.mc(array(beta.diff*V.inv[[length(V.inv)]], p), V.inv[2:(length(V.inv) - 1)]))*beta.diff/2)

  } else if (is.matrix(V.inv)) {

    c3 <- sum(crossprod(beta - beta.tilde, crossprod(V.inv, beta - beta.tilde)))/2

  } else if (class(V.inv) == "dsyMatrix" | class(V.inv) == "dpoMatrix") {

    c3 <- sum(Matrix::crossprod(beta - beta.tilde, Matrix::crossprod(V.inv, beta - beta.tilde)))/2

  }
  else if (is.vector(V.inv)) {
    c3 <- sum((beta - beta.tilde)^2*V.inv)/2
  }
  comp.sum <- c(c1, c2, c3)

  return(sum(comp.sum))
}

h.log.r.sng <- function(eta, d0, d1, Omega.inv, beta, c, r.tilde, V.r.inv) {
  r = d0*sin(eta) + d1*cos(eta)
  s.sq <- 1/r^2
  p <- unlist(lapply(Omega.inv, function(x) {nrow(x)}))
  c1 <- -log(abs(r)) + (2*c - 1)*log(sqrt(s.sq))
  c2 <- -c*s.sq
  # Need to improve how this is coded

  c3 <- 0
  for (i in 1:length(r)) {
    c3 <- c3 - sum((get.kron.row(i, Omega.inv)*abs(r)*beta)[-i]*abs(r[i])*beta[i])/2
  }
  comp.sum <- c(c1, c2, c3)
  return(sum(comp.sum))

}

h.log.r.spb <- function(eta, d0, d1, Omega.inv, beta, c, r.tilde, V.r.inv, deltas) {
  r = d0*sin(eta) + d1*cos(eta)
  p <- unlist(lapply(Omega.inv, function(x) {nrow(x)}))
  alpha <- c/2
  s.sq <- 1/r^2
  c1 <- -log(abs(r)) + ((1 + alpha)/(1 - alpha) - 1)*log(sqrt(s.sq))
  c2 <- -(((2*gamma(3/c))/gamma(1/c))^(alpha/(1 - alpha))*f.deltas(deltas = deltas, c = c))*(s.sq)^(alpha/(1 - alpha))

  c3 <- 0
  for (i in 1:length(r)) {
    c3 <- c3 - sum((get.kron.row(i, Omega.inv)*abs(r)*beta)[-i]*abs(r[i])*beta[i])/2
  }
  comp.sum <- c(c1, c2, c3)
  return(sum(comp.sum))

}

h.log.r.spn <- function(eta, d0, d1, Omega.inv, beta, Psi.inv, r.tilde, V.r.inv) {
  r = d0*sin(eta) + d1*cos(eta)
  p <- unlist(lapply(Omega.inv, function(x) {nrow(x)}))

  c1 <- -log(r^2)/2
  c2 <- -sum(c(atrans.mc(array(1/r, dim = p), Psi.inv))*(1/r))/2
  c3 <- 0
  for (i in 1:length(r)) {
    c3 <- c3 - sum((get.kron.row(i, Omega.inv)*r*beta)[-i]*r[i]*beta[i])/2
  }
  comp.sum <- c(c1, c2, c3)

  return(sum(comp.sum))

}

sample.d <- function(theta, delta, V.half) {

  if (is.list(V.half)) {
    q <- nrow(V.half[[1]])
    p <- unlist(lapply(V.half[-1], nrow))
    gamma <- crossprod(V.half[[1]], rnorm(q))
    beta <- c(atrans.mc(array(rnorm(prod(p)), p), V.half[2:(length(V.half) - 1)]))*V.half[[length(V.half)]]
    d <- c(gamma, beta)
  } else if (is.matrix(V.half)) {
    d <- crossprod(V.half, rnorm(nrow(V.half)))
  } else if (class(V.half) == "Cholesky") {
    d <- Matrix::crossprod(V.half, rnorm(nrow(V.half)))
  } else if (is.vector(V.half)) {
    d <- V.half*rnorm(length(V.half))
  }
  d0 <- delta*sin(theta) + d*cos(theta)
  d1 <- delta*cos(theta) - d*sin(theta)
  return(list("d0" = d0, "d1" = d1))
}

sample.r.eta <- function(r, Omega.inv, beta, c = NULL, eta, r.tilde, V.r.inv, V.r.half, prior, Psi.inv = NULL,
                         q = NULL, deltas = NULL) {
  delta <- r - r.tilde
  d <- sample.d(V.half = V.r.half, theta = eta, delta = delta)
  if (prior == "sng") {
    ll.args <- list("d0" = d$d0,
                    "d1" = d$d1,
                    "Omega.inv" = Omega.inv,
                    "beta" = beta,
                    "c" = c,
                    "r.tilde" = r.tilde,
                    "V.r.inv" = V.r.inv)
    ll.fun <- "h.log.r.sng"
  } else if (prior == "spn") {
    ll.args <- list("d0" = d$d0,
                    "d1" = d$d1,
                    "Omega.inv" = Omega.inv,
                    "beta" = beta,
                    "Psi.inv" = Psi.inv,
                    "r.tilde" = r.tilde,
                    "V.r.inv" = V.r.inv)
    ll.fun <- "h.log.r.spn"
  } else if (prior == "spb") {
    ll.args <- list("d0" = d$d0,
                    "d1" = d$d1,
                    "Omega.inv" = Omega.inv,
                    "beta" = beta,
                    "c" = c,
                    "deltas" = deltas,
                    "r.tilde" = r.tilde,
                    "V.r.inv" = V.r.inv)
    ll.fun <- "h.log.r.spb"
  }
  eta <- slice(x.tilde = eta, ll.fun = ll.fun,
               var.lim = c(0, 2*pi),
               ll.args = ll.args)
  delta <- d$d0*sin(eta) + d$d1*cos(eta)
  return(list("eta" = eta,
              "r" = r.tilde + delta))
}

sample.beta.theta <- function(X, U, y, V.half, beta, theta, beta.tilde, Omega.inv, V.inv,
                              sig.sq, reg) {
  delta <- beta - beta.tilde
  d <- sample.d(V.half = V.half, theta = theta, delta = delta)
  theta <- slice(x.tilde = theta, ll.fun = "h.log", var.lim = c(0, 2*pi),
                 ll.args = list("X" = X,
                                "U" = U,
                                "y" = y,
                                "d0" = d$d0,
                                "d1" = d$d1,
                                "beta.tilde" = beta.tilde,
                                "Omega.inv" = Omega.inv,
                                "V.inv" = V.inv,
                                "sig.sq" = sig.sq,
                                "reg" = reg))
  delta <- d$d0*sin(theta) + d$d1*cos(theta)
  return(list("theta" = theta,
              "beta" = beta.tilde + delta))
}

#' @export
sampler <- function(X, y, Omega.half = NULL,
                    U = NULL, # Matrix of Unpenalized covariates
                    num.samp = 100,
                    print.iter = TRUE,
                    max.iter = 1000,
                    eps = 10^(-12),
                    diag.app = FALSE,
                    burn.in = 0, prior = "sno", c = 1, Psi.half = NULL, sig.sq = ifelse(reg == "logit", 1, NULL), reg = "linear",
                    fix.beta = FALSE, beta.fix = rep(0, prod(dim(X)[-1]) + ifelse(is.null(U), 1, ncol(U))),
                    rho = 0, pr.rho.a = 1, pr.rho.b = 1, tune = 0.5,
                    from.prior = FALSE,
                    do.svd = TRUE, slice = TRUE, joint.beta = list(prod(dim(X)[-1]) + ifelse(is.null(U), 1, ncol(U))),
                    str = "uns", thin = 1,
                    pr.Omega.V.inv = lapply(dim(X)[-1], diag),
                    pr.Psi.V.inv = lapply(dim(X)[-1], diag),
                    pr.Omega.df = lapply(dim(X)[-1], function(x) {x + 2}),
                    pr.Psi.df = lapply(dim(X)[-1], function(x) {x + 2}),
                    pr.sig.sq.shape = 3/2,
                    pr.sig.sq.rate = 1/2) {

  # Record some quantities and set up objects to save results in
  pr.Omega.V.inv <- pr.Omega.V.inv
  pr.Psi.V.inv <- pr.Psi.V.inv

  pr.Omega.df <- pr.Omega.df
  pr.Psi.df <- pr.Psi.df

  n <- length(y)
  p <- dim(X)[-1]
  if (!fix.beta) {
    beta.fix <- rep(0, prod(dim(X)[-1]) + ifelse(is.null(U), 1, ncol(U)))
  }
  X.arr <- X
  X <- t(apply(X, 1, "c"))

  # Set up indicators for null arguments, will be used to decide whether or not to resample
  null.Omega.half <- unlist(lapply(Omega.half, function(x) {is.null(x)}))


  # Set starting values
  null.rho <- is.null(rho)
  null.sig.sq <- is.null(sig.sq)
  if (null.rho) {
    rho <- 0
    rho.psi <- 0
  }
  if (is.null(sig.sq)) {
    sig.sq <- 1
  }
  # Specification of rho for first matrix overrides specification of Omega.half
  null.Omega.half[[1]] <- null.rho
  null.U <- is.null(U)
  if (max(null.Omega.half) == 1) {
    res.Omega <- vector("list", length(p))
    res.Sigma <- vector("list", length(p))
    if (prior == "spn") {
      res.Psi <- vector("list", length(p))
    }
    for (i in which(null.Omega.half)) {
      if (i != 1) {
        Omega.half[[i]] <- diag(p[i])
      } else {
        Omega.half[[i]] <- sym.sq.root(make.ar.mat(p = p[i], rho = rho, inv = FALSE))
      }
      res.Omega[[i]] <- array(dim = c(num.samp, p[i], p[i]))
      res.Sigma[[i]] <- array(dim = c(num.samp, p[i], p[i]))
      if (prior == "spn") {
        if (i != 1) {
          Psi.half[[i]] <- diag(p[i])
        } else {
          Psi.half[[i]] <- sym.sq.root(make.ar.mat(p = p[i], rho = rho.psi, inv = FALSE))
        }
        res.Psi[[i]] <- array(dim = c(num.samp, p[i], p[i]))
      }
    }
  }
  Omega <- lapply(Omega.half, function(x) {crossprod(x)})
  Omega.inv <- lapply(Omega, function(x) {ei.inv(x)})
  Omega.half.inv <- lapply(Omega.half, function(x) {ei.inv(x)})

  if (prior == "spn") {
    Psi <- lapply(Psi.half, function(x) {crossprod(x)})
    Psi.inv <- lapply(Psi, function(x) {ei.inv(x)})
  }


  # No intercept if U is null
  if (null.U) {
    U <- matrix(0, nrow = n, ncol = 1)
  }
  q <- ncol(U)

  # Results objects
  # acc.c <- array(dim = c(num.samp, 1))
  # res.c <- array(dim = c(num.samp, 1))
  res.rho <- array(dim = c(num.samp, 1))
  res.rho.psi <- array(dim = c(num.samp, 1))
  res.sig.sq <- array(dim = c(num.samp, 1))
  res.B <- array(dim = c(num.samp, prod(p)))
  res.S <- array(dim = c(num.samp, prod(p)))
  res.theta <- numeric(num.samp)
  res.ome <- array(dim = c(num.samp, n))
  res.eta <- numeric(num.samp)
  res.gamma <- array(dim = c(num.samp, q))
  res.D <- array(dim = c(num.samp, prod(p)))

  null.c <- is.null(c)
  if (null.c) {
    c <- 1
  }
  S <- array(1, dim = p)
  gamma <- beta.fix[1:q]
  beta <- beta.fix[(q + 1):length(beta.fix)]
  B <- array(beta, dim = p)
  Z <- array(0, dim = p)
  R <- array(1, dim = p)
  if (slice) {
    theta <- pi
  } else {
    if (reg == "logit" | reg == "nb") {
      ome <- rep(1, n)
      omeD <- Matrix(0, nrow = n, ncol = n, sparse = TRUE)
      diag(omeD) <- ome
      if (reg == "logit") {
        offset <- rep(1/2, n)
      }
    } else {
      offset <- rep(0, n)
    }
  }
  eta <- pi

  X.arr.s <- X.arr
  X.arr.s <- array(c(crossprod(apply(X.arr, 1, "c"), diag(c(S)))), dim = c(n, p))
  for (l in 1:length(p)) {
    X.arr.s <- amprod.mc(X.arr.s, Omega.half[[l]], l + 1)
  }
  W <- t(apply(X.arr.s, 1, "c"))
  UW <- cbind(U, W)

  penC <- c(rep(0, ncol(U)), rep(1, ncol(W)))

  for (i in 1:(burn.in + thin*num.samp)) {
    if (print.iter) {cat("i=", i, "\n")}

    if (!fix.beta) {
      sample.beta <- c(gamma, c(Z))
      if (print.iter) {cat("Sample Beta\n")}
      if (slice) {

        # Set mean for proposal distribution
        if (from.prior) {
          z.tilde <- rep(0, ncol(UW))
        } else if ((i == 1 & prior == "sno" & max(null.Omega.half[-1]) == 0 & (!null.rho & !null.Omega.half[1]) & reg == "linear" & !null.sig.sq) |
                   (i >= 1 & (prior != "sno" | max(null.Omega.half[-1]) != 0 | (null.rho | null.Omega.half[1]) | reg != "linear" | null.sig.sq))) {
          if (reg == "logit") {
            if (print.iter) {cat("Get Mode\n")}

            z.tilde <- coord.desc.logit(y = y, X = UW, Omega.inv = penC,
                                        print.iter = FALSE, max.iter = max.iter, eps = eps,
                                        start.beta = rep(0, ncol(UW)),
                                        joint.beta = joint.beta)$beta
          } else if (reg == "linear") {
            if (print.iter) {cat("Get Mode\n")}
            z.tilde <- coord.desc.lin(y = y, X = UW, sig.sq = sig.sq, Omega.inv = penC,
                                      print.iter = FALSE, max.iter = max.iter, eps = eps,
                                      start.beta = rep(0, ncol(UW)))$beta
          }
        }
        if (from.prior) {
          V.half <- c(rep(10^(12), q), rep(1, length(z.tilde) - q))
          V.inv <- 1/V.half^2
        } else {
          if  (print.iter) {cat("Get Pieces for Covariance Matrix\n")}
          UWz.tilde <- crossprod(t(UW), z.tilde)[, 1]
          if (reg == "logit") {


            AA <- diag(exp(UWz.tilde)/(1 + exp(UWz.tilde))^2)
            AAU <- crossprod(sqrt(AA), U)
            if (!diag.app) {
              BB <- crossprod(AA, UW)
              if (do.svd) {
                AAX <- mat(amprod.mc(X.arr.s, sqrt(AA), 1), 1)
              }
            }

          } else {
            # AA <- diag(length(UWz.tilde))
            AAU <- U
            if (!diag.app) {
              BB <- UW
              if (do.svd) {
                AAX <- mat(X.arr.s, 1)
              }
            }
          }
          if (diag.app) {
            if (do.svd) {
              if (n < prod(p)) {

                svd.A <- svd(AAX)
                R.A <- svd.A$v
                d.A <- svd.A$d
                V.half <- sqrt(c(diag(solve(crossprod(AAU)/sig.sq)), rep(1, prod(p)) - apply(R.A, 1, function(x) {sum(x^2*(d.A^2/(1 + d.A^2)))})/sig.sq))

              } else {
                V.half <- sqrt(c(diag(solve(crossprod(AAU)/sig.sq)), diag(solve(diag(prod(p)) + crossprod(AAX)/sig.sq))))
              }
            } else {
              V.half <- sqrt(c(diag(solve(crossprod(AAU)/sig.sq)), rep(1, prod(p))))
            }

            V.inv <- 1/V.half^2

          } else {
            if (print.iter) {cat("Get Covariance Matrix\n")}
            UWtBB <- crossprod(UW, BB)/sig.sq
            V.inv <- UWtBB + diag(penC)
            V.half <- sym.sq.root.inv(V.inv)
          }

        }

        sample <- sample.beta.theta(X = W, U = U, y = y, V.half = V.half, beta = c(gamma, c(Z)),
                                    theta = theta, beta.tilde = z.tilde,
                                    Omega.inv = lapply(p, function(x) {diag(x)}), V.inv = V.inv,
                                    sig.sq = sig.sq, reg = reg)
        sample.beta <- sample$beta
        theta <- sample$theta
      } else {
        if (reg == "logit") {
          ome <- BayesLogit::rpg(n, offset*2, crossprod(t(UW), c(gamma, c(Z))))
          diag(omeD) <- ome
        }
        if (length(joint.beta) == 1) {
          if (reg == "logit") {

            UWtUW <- Matrix::crossprod(UW, Matrix::crossprod(omeD, UW))

          } else if ((i == 1 & reg == "linear" & max(null.Omega.half) == 0 & prior == "sno") |
                     !(reg == "linear") | !(max(null.Omega.half) == 0) | !(prior == "sno")) {
            UWtUW <- crossprod(UW)

          }
          UWty <- crossprod(UW, y - offset)

          if (reg == "linear") {


            if (min(U) == 0 & max(U) == 0) {
              sample.beta[(q + 1):nrow(UWty)] <- samp.beta(XtX = UWtUW[(q + 1):nrow(UWtUW), (q + 1):ncol(UWtUW)],
                                                           Xty = UWty[(q + 1):nrow(UWty)],
                                                           Omega.inv = diag(penC[(q + 1):nrow(UWty)]), sig.sq = sig.sq)
            } else {
              sample.beta <- samp.beta(XtX = UWtUW, Xty = UWty,
                                       Omega.inv = diag(penC), sig.sq = sig.sq)
            }


          } else {
            if (min(U) == 0 & max(U) == 0) {

              sample.beta[(q + 1):nrow(UWty)] <- samp.beta(XtX = UWtUW[(q + 1):nrow(UWtUW), (q + 1):ncol(UWtUW)],
                                                           Xty = UWty[(q + 1):nrow(UWty)],
                                                           Omega.inv = diag(penC[(q + 1):nrow(UWty)]), sig.sq = 1)
            } else {
              sample.beta <- samp.beta(XtX = UWtUW, Xty = UWty,
                                       Omega.inv = diag(penC), sig.sq = 1)
            }


          }

        } else {

          for (block in joint.beta) {

            not.block <- (1:ncol(UW))[!1:ncol(UW) %in% block]
            if (reg == "logit") {

              UWtUW <- Matrix::crossprod(UW[, block], Matrix::crossprod(omeD, UW[, block]))

              UWtNUWZ <- Matrix::crossprod(Matrix::t(Matrix::crossprod(UW[, block],
                                                                       Matrix::crossprod(omeD, UW[, not.block]))), sample.beta[ not.block])

            } else if ((i == 1 & reg == "linear" & max(null.Omega.half) == 0 & prior == "sno") |
                       !(reg == "linear") | !(max(null.Omega.half) == 0) | !(prior == "sno")) {
              UWtUW <- crossprod(UW[, block])
              UWtNUWZ <- crossprod(UW[, block], crossprod(t(UW[, not.block]), sample.beta[not.block]))

            }


            UWty <- crossprod(UW[, block], y - offset)

            if (reg == "linear") {
              sample.beta[block] <- samp.beta(XtX = UWtUW, Xty = UWty - UWtNUWZ,
                                              Omega.inv = diag(penC[block]), sig.sq = sig.sq)
            } else {
              sample.beta[block] <- samp.beta(XtX = UWtUW, Xty = UWty - UWtNUWZ,
                                              Omega.inv = diag(penC[block]), sig.sq = 1)
            }
          }

        }

      }

      gamma <- sample.beta[1:q]
      Z <- array(sample.beta[(q + 1):length(sample.beta)], dim = p)
      S <- array(S, p)
      Z <- array(Z, p)
      B <- S*atrans.mc(Z, Omega.half)


    }

    if (prior == "sng" | prior == "spn" | prior == "spb") {

      r.tilde <- rep(0, prod(p))
      Omega.inv.diag <- rep(diag(Omega.inv[[1]]), times = prod(p[-1]))
      for (k in 2:length(p)) {
        Omega.inv.diag <- rep(diag(Omega.inv[[k]]), each = prod(p[1:(k - 1)]))*Omega.inv.diag
      }
      V.r.inv <- as.vector(Omega.inv.diag*c(B^2))
      V.r.half <- sqrt(1/V.r.inv)

      if (prior == "sng") {
        if (print.iter) {cat("Set Sampling Values for R\n")}
        if (print.iter) {cat("Sample R\n")}
        sample <- sample.r.eta(r = c(R), Omega.inv = Omega.inv, beta = c(B), c = c, eta = eta,
                               r.tilde = r.tilde, V.r.inv = V.r.inv, V.r.half = V.r.half, prior = prior)

      } else if (prior == "spn") {
        if (print.iter) {cat("Set Sampling Values for R\n")}
        if (print.iter) {cat("Sample R\n")}
        sample <- sample.r.eta(r = c(R), Omega.inv = Omega.inv, beta = c(B), Psi.inv = Psi.inv, eta = eta,
                               r.tilde = r.tilde, V.r.inv = V.r.inv, V.r.half = V.r.half, c = c, prior = prior)
      } else if (prior == "spb") {
        if (print.iter) {cat("Sample Delta\n")}
        if (i == 1) {
          deltas <- runif(prod(p), 0, pi)
        } else {
          for (ii in 1:length(deltas)) {
            alpha <- c/2
            xi <- (2*gamma(3/(2*alpha))*c(S)[ii]^2/gamma(1/(2*alpha)))^(alpha/(1 - alpha))
            deltas[ii] <- slice(x.tilde = deltas[ii], ll.fun = "g.delta",
                                var.lim = c(0, pi), ll.args = list("c" = c, "xi" = xi))
          }
        }
        if (print.iter) {cat("Set Sampling Values for R\n")}
        if (print.iter) {cat("Sample R\n")}

        sample <- sample.r.eta(r = c(R), Omega.inv = Omega.inv, beta = c(B), deltas = deltas, eta = eta,
                               r.tilde = r.tilde, V.r.inv = V.r.inv, V.r.half =
                                 V.r.half, c = c, prior = prior)
      }
      r <- sample$r
      R <- array(r, p)
      eta <- sample$eta
      S <- array(1/r, dim = p)
      if (!fix.beta) {
        Z <- atrans.mc(B/S, Omega.half.inv)
        B <- S*atrans.mc(Z, Omega.half)
      }

    }

    for (k in which(null.Omega.half)) {
      Covar <- B/S
      for (l in 1:length(p)) {
        if (l != k) {
          Covar <- amprod.mc(Covar, Omega.inv[[l]], l)
        }
      }
      if (k == 1 & null.rho) {
        if (print.iter) {cat("Sample rho\n")}
        rho <- slice(x.tilde = rho, ll.fun = "cond.rho.log", var.lim = c(-1, 1),
                     ll.args = list("B" = Covar,
                                    "pr.a" = pr.rho.a,
                                    "pr.b" = pr.rho.b,
                                    "j" = k))
        Omega.inv[[k]] <- make.ar.mat(p = p[k], rho = rho, inv = TRUE)

      } else {

        # Covar <- t(apply(Covar, k, "c"))
        # Omega.inv[[k]] <- rWishart(1, prod(p[-k]) + p[k] + 2, solve(tcrossprod(Covar) + diag(p[k])))[, , 1]
        Covar <- apply(Covar, k, "c")
        Omega.inv[[k]] <- samp.Omega.inv(Beta = Covar, str = str,
                                         pr.V.inv = pr.Omega.V.inv[[k]],
                                         pr.df = pr.Omega.df[[k]])
      }
      Omega[[k]] <- ei.inv(Omega.inv[[k]])
      if (i > burn.in & (i - burn.in)%%thin == 0) {
        res.Omega[[k]][(i - burn.in)/thin, , ] <- Omega[[k]]
      }

      Omega.half[[k]] <- sym.sq.root(Omega[[k]])
      Omega.half.inv[[k]] <- sym.sq.root(Omega.inv[[k]])

      if (prior == "spn") {

        Covar <- S
        for (l in 1:length(p)) {
          if (l != k) {
            Covar <- amprod.mc(Covar, Psi.inv[[l]], l)
          }
        }

        if (k == 1 & null.rho) {
          if (print.iter) {cat("Sample rho psi\n")}
          rho.psi <- slice(x.tilde = rho.psi, ll.fun = "cond.rho.log", var.lim = c(-1, 1),
                           ll.args = list("B" = Covar,
                                          "pr.a" = pr.rho.a,
                                          "pr.b" = pr.rho.b,
                                          "j" = k))
          # print(rho.psi)
          Psi.inv[[k]] <- make.ar.mat(p = p[k], rho = rho.psi, inv = TRUE)
        } else {
          # Covar <- t(apply(Covar, k, "c"))
          # Psi.inv[[k]] <- rWishart(1, prod(p[-k]) + p[k] + 2, solve(tcrossprod(Covar) + diag(p[k])))[, , 1]
          Covar <- apply(Covar, k, "c")
          Psi.inv[[k]] <- samp.Omega.inv(Beta = Covar, str = str,
                                         pr.V.inv = pr.Psi.V.inv[[k]],
                                         pr.df = pr.Psi.df[[k]])
        }
        Psi[[k]] <- ei.inv(Psi.inv[[k]])
        if (i > burn.in & (i - burn.in)%%thin == 0) {
          res.Psi[[k]][(i - burn.in)/thin, , ] <- Psi[[k]]
        }
        Psi.half[[k]] <- sym.sq.root(Psi[[k]])
      }

      if (reg == "linear" & null.sig.sq) {
        resid <- y - crossprod(t(U), gamma) - crossprod(t(X), c(B))
        sig.sq <- 1/rgamma(1, shape = pr.sig.sq.shape + length(resid)/2, rate = pr.sig.sq.rate + sum(resid^2)/2)
      }

      if (i > burn.in & (i - burn.in)%%thin == 0) {
        if (prior == "sno") {
          res.Sigma[[k]][(i - burn.in)/thin, , ] <- Omega[[k]]
        } else if (prior == "sng") {
          es <- c^(-1/2)*gamma(c + 1/2)/gamma(c)
          els <- (1 - es^2)*diag(p[k]) + es^2*matrix(1, nrow = p[k], ncol = p[k])
          res.Sigma[[k]][(i - burn.in)/thin, , ] <- Omega[[k]]*els
        } else if (prior == "spb") {
          es <- sqrt(pi/2)*gamma(2/c)/(sqrt(gamma(1/c)*gamma(3/c)))
          els <- (1 - es^2)*diag(p[k]) + es^2*matrix(1, nrow = p[k], ncol = p[k])
          res.Sigma[[k]][(i - burn.in)/thin, , ] <- Omega[[k]]*els
        } else if (prior == "spn") {
          res.Sigma[[k]][(i - burn.in)/thin, , ] <- Omega[[k]]*Psi[[k]]
          # print(Omega[[k]]*Psi[[k]])
        }
      }
    }


    if (!fix.beta & max(null.Omega.half) == 1) {
      Z <- atrans.mc(B/S, Omega.half.inv)
      B <- S*atrans.mc(Z, Omega.half)
    }


    if (prior != "sno" | max(null.Omega.half) == 1) {
      X.arr.s <- array(c(crossprod(apply(X.arr, 1, "c"), diag(c(S)))), dim = c(n, p))
      for (l in 1:length(p)) {
        X.arr.s <- amprod.mc(X.arr.s, Omega.half[[l]], l + 1)
      }
      W <- t(apply(X.arr.s, 1, "c"))
      UW <- cbind(U, W)

    }

    if (prior == "sng") {
      if (null.c) {
        sample.c <- sample.c(S = S, c.old = c, tune = tune)
        c <- sample.c$c
        a.c <- sample.c$acc
      }
    }

    if (i > burn.in & (i - burn.in)%%thin == 0) {

      res.B[(i - burn.in)/thin, ] <- c(B)
      res.gamma[(i - burn.in)/thin, ] <- gamma
      if (slice) {
        res.theta[(i - burn.in)/thin] <- theta
      }
      res.S[(i - burn.in)/thin, ] <- c(S)
      res.eta[(i - burn.in)/thin] <- eta
      # if (prior == "sng" & null.c) {
      #   res.c[(i - burn.in)/thin] <- c
      #   acc.c[(i - burn.in)/thin] <- a.c
      # }
      if (prior == "spb") {
        res.D[(i - burn.in)/thin, ] <- c(deltas)
      }
      if (null.rho & null.Omega.half[1]) {
        res.rho[(i - burn.in)/thin] <- rho
        if (prior == "spn") {
          res.rho.psi[(i - burn.in)/thin] <- rho.psi
        }
      }
      if (reg == "linear" & null.sig.sq) {
        res.sig.sq[(i - burn.in)/thin] <- sig.sq
      }
    }
  }

  res.list <- list("Bs" = res.B, "gammas" = res.gamma, "etas" = res.eta,
                   "Ss" = res.S)
  if (slice) {
    res.list[["thetas"]] <- res.theta
  } else {
    res.list[["omes"]] <- res.ome
  }
  if (prior == "spb") {
    res.list[["Ds"]] <- res.D
  }
  if (max(null.Omega.half) == 1) {
    res.list[["Omegas"]] <- res.Omega
    res.list[["Sigmas"]] <- res.Sigma
    res.list[["rhos"]] <- res.rho
  }
  if (max(null.Omega.half) == 1 & prior == "spn") {
    res.list[["Psis"]] <- res.Psi
    res.list[["rho.psis"]] <- res.rho.psi
  }
  if (reg == "linear" & null.sig.sq) {
    res.list[["sig.sqs"]] <- res.sig.sq
  }

  return(res.list)
}



#' @export
em.est <- function(X, y, Omega.half,
                   U = NULL, # Matrix of Unpenalized covariates
                   num.samp = 100, #
                   print.iter = TRUE,
                   max.iter.slice = 1000,
                   eps.slice = 10^(-12),
                   max.iter.em = NULL,
                   eps.em = 10^(-3),
                   diag.app = FALSE,
                   burn.in = 0, prior = "sno", c = 1, Psi.half = NULL, sig.sq = NULL, reg = "linear",
                   rho = 0) {

  W <- t(apply(X, 1, "c"))

  if (is.null(U)) {
    UW <- cbind(rep(0, nrow(W)), W)
  } else {
    UW <- cbind(U, W)
  }

  penC <- matrix(0, nrow = ncol(UW), ncol = ncol(UW))

  Omega.inv <- lapply(Omega.half, function(x) {ei.inv(crossprod(x))})
  O.i <- matrix(1, nrow = 1, ncol = 1)
  for (i in 1:length(Omega.inv)) {
    O.i <- Omega.inv[[i]]%x%O.i
  }

  fix.beta = FALSE;
  if (is.null(max.iter.em)) {
    if (length(num.samp) > 1) {
      max.iter.em <- length(num.samp) - 1
    } else {
      max.iter.em <- 1
    }
  }
  if (length(num.samp) == 1) {
    num.samp <- rep(num.samp, max.iter.em + 1)
  }
  if (print.iter) {cat("Set Starting Value\n")}
  # Get initial values
  samples <- sampler(X = X, y = y, Omega.half = Omega.half, num.samp = num.samp[1], print.iter = FALSE,
                     max.iter = max.iter.slice, eps = eps.slice, diag.app = diag.app, burn.in = burn.in, prior = prior, c = c,
                     U = U, Psi.half = Psi.half, sig.sq = sig.sq, reg = reg, fix.beta = FALSE, rho = rho)
  post.mean <- c(colMeans(samples$gammas), colMeans(samples$Bs))
  post.median <- c(apply(samples$gammas, 2, median), apply(samples$Bs, 2, median))

  fix.beta = TRUE;
  beta.fix <- post.mean

  betas <- matrix(nrow = max.iter.em, ncol = prod(dim(X)[-1]) + ifelse(is.null(U), 0, ncol(U)))

  for (i in 1:max.iter.em) {

    if (print.iter) {cat("EM Iteration: ", i, "\n")}

    beta.fix[beta.fix == 0] <- rnorm(sum(beta.fix == 0))

    samples.diag <- sampler(X = X, y = y, Omega.half = Omega.half, num.samp = num.samp[i + 1], print.iter = FALSE,
                            max.iter = max.iter.slice, eps = eps.slice, diag.app = diag.app, burn.in = burn.in, prior = prior, c = c,
                            U = U, Psi.half = Psi.half, sig.sq = sig.sq, reg = reg, fix.beta = fix.beta,
                            beta.fix = beta.fix, rho = rho)

    if (prior != "spn") {
      inv.ss <- matrix(rowMeans(apply(samples.diag$Ss, 1, function(x) {tcrossprod(1/abs(x))})), nrow = prod(p), ncol = prod(p))
    } else {
      inv.ss <- matrix(rowMeans(apply(samples.diag$Ss, 1, function(x) {tcrossprod(1/x)})), nrow = prod(p), ncol = prod(p))
    }


    penC[2:nrow(penC), 2:ncol(penC)] <- (O.i*inv.ss)

    if (reg == "logit") {
      beta.fix <- coord.desc.logit(y = y, X = UW, Omega.inv = penC,
                                   print.iter = FALSE, max.iter = max.iter.slice, eps = eps.slice)$beta
    } else if (reg == "linear") {
      beta.fix <- coord.desc.lin(y = y, X = UW, Omega.inv = penC,
                                 print.iter = FALSE, max.iter = max.iter.slice, eps = eps.slice, sig.sq = sig.sq)$beta
    }

    if (is.null(U)) {
      betas[i, ] <- beta.fix[-1]
    } else {
      betas[i, ] <- beta.fix
    }

    if (i > 1) {
      if (mean(abs(betas[i - 1, ] - betas[i, ])) < eps.em) {
        break
      }
    }
  }

  if (is.null(U)) {
    post.mean <- post.mean[-1]
    post.median <- post.median[-1]
  }

  return(list("post.mean" = post.mean, "post.med" = post.median, betas = betas[1:i, ]))

}

