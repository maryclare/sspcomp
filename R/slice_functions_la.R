#' @export
get.beta.blocks <- function(X, U = NULL, min.block.size = 25, no.eig = TRUE) {

  p <- dim(X)[-1]
  if (is.null(U)) {
    q <- 0
  } else {
    q <- dim(U)[-1]
  }
  if (no.eig) {
    ends <- seq(1, prod(p) + q, by = min.block.size)
    joint.beta <- vector("list", length = length(ends) - 1)
    ends[-1] <- ends[-1] + (prod(p) + q - max(ends[-1]))

    for (i in 2:(length(ends))) {
      if (i == 2) {
        joint.beta[[i - 1]] <- ends[i - 1]:(ends[i])
      } else {
        joint.beta[[i - 1]] <- (ends[i - 1] + 1):(ends[i])
      }
    }
  } else {
    # Fix this later
    nc <- prod(p) + q
    if (!is.null(U)) {
      C <- crossprod(cbind(t(apply(X, 1, "c")), U))
    } else {
      C <- crossprod(cbind(t(apply(X, 1, "c"))))
    }
    C.ei <- eigen(C)
    # Get number of covariates per block
    num.per.block <- round(nc*C.ei$values/sum(C.ei$values), 0)
    # Decide how many groups to make
    max.ei <- max(which((num.per.block) > 0))
    num.block <- max.ei + (nc - sum(num.per.block[1:max.ei]))
    new.num.per.block <- numeric(num.block)
    new.num.per.block[1:max.ei] <- num.per.block[1:max.ei]
    new.num.per.block[(max.ei + 1):length(new.num.per.block)] <- 1
    num.per.block <- new.num.per.block
    joint.beta <- vector("list", length = num.block)

    remain <- 1:(nc)
    for (i in 1:max.ei) {

        order.vec <- order(abs(C.ei$vectors[i, ]), decreasing = FALSE)
        order.vec <- order.vec[order.vec %in% remain]

        block.beta <- order.vec[1:num.per.block[i]]

        remain <- remain[!remain %in% block.beta]

        joint.beta[[i]] <- block.beta

    }
    for (i in 1:length(remain)) {
      joint.beta[[max.ei + i]] <- remain[i]
    }

    }
  return(joint.beta)
}

#' @export
samp.Omega.inv <- function(Beta, pr.V.inv = diag(ncol(Beta)),
                           pr.df = ncol(Beta) + 2, str = "uns") {
  p <- ncol(Beta)
  if (str == "uns") {
    V.inv <- crossprod(Beta) + pr.V.inv
    df <- nrow(Beta) + pr.df
    # print(V.inv)
    V.half <- sym.sq.root.inv((V.inv + t(V.inv))/2)

    return(tcrossprod(crossprod(V.half, matrix(rnorm(p*df), nrow = p, ncol = df))))

    # return(matrix(rWishart(1, df, solve(V.inv))[, , 1], nrow = p, ncol = p))
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
                                                diag(sqrt(ifelse(V.inv.eig$values > 0, 1/V.inv.eig$values, 0)), nrow = length(V.inv.eig$values), ncol = length(V.inv.eig$values))),
                             V.inv.eig$vectors)
  V <- Matrix::tcrossprod(Matrix::tcrossprod(V.inv.eig$vectors,
                                             diag(ifelse(V.inv.eig$values > 0, 1/V.inv.eig$values, 0), nrow = length(V.inv.eig$values), ncol = length(V.inv.eig$values))),
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
    ARMat <- diag(p, nrow = p, ncol = p)
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

g.delta <- function(theta, xi, c) {

  f.val <- f.deltas(deltas = theta, c = c)
  log(f.val) - xi*f.val
}

# Performs univariate slice sampling for an arbitrary RV, requires
slice <- function(x.tilde, # Previous value of theta
                  ll.fun,  # Log likelihood function
                  var.lim, # Domain limits of slice variable theta
                  ll.args  # Arguments required by the log likelihood function
) {

  d <- x.tilde
  x.tilde.vals <- unique(x.tilde)

  arg.list <- vector("list", length = length(ll.args) + 1)
  names(arg.list)[1] <- "theta"
  names(arg.list)[2:length(arg.list)] <- names(ll.args)
  arg.list[2:length(arg.list)] <- ll.args
  arg.list[["theta"]] <- d

  for (i in 1:length(x.tilde.vals)) {

    z <- do.call(ll.fun, arg.list) - rexp(1)
    a <- var.lim[1]
    b <- var.lim[2]
    d[x.tilde == x.tilde.vals[i]] <- runif(1, a, b)
    arg.list[["theta"]] <- d

    ll.d <- do.call(ll.fun, arg.list)
    while (z > ll.d) {
      if (unique(d[x.tilde == x.tilde.vals[i]]) < x.tilde.vals[i]) {
        a <- unique(d[x.tilde == x.tilde.vals[i]])
      } else {
        b <- unique(d[x.tilde == x.tilde.vals[i]])
      }
      d[x.tilde == x.tilde.vals[i]] <- runif(1, a, b)

      arg.list[["theta"]] <- d

      ll.d <- do.call(ll.fun, arg.list)
      # If we end up in a bad spot, reinitalize interval choice
      if (is.nan(ll.d)) {
        a <- var.lim[1]
        b <- var.lim[2]
        d[x.tilde == x.tilde.vals[i]] <- runif(1, a, b)
        arg.list[["theta"]] <- d
        ll.d <- do.call(ll.fun, arg.list)
      }

    }
    x.tilde.p <- d
  }
  return(x.tilde.p)
}

#' @export
cond.rho.log <- function(theta, B, pr.a, pr.b, j) {
  p <- dim(B)
  B.mat <- mat(B, j)

  c1 <- -(prod(p[-j])*(p[j] - 1)/2)*log(1 - theta^2)

  O.i <- make.ar.mat(p = p[j], rho = theta, inv = TRUE)

  c2 <- -sum(diag(crossprod(B.mat, crossprod(O.i, B.mat))))/2

  c3 <- dbeta((theta + 1)/2, pr.a, pr.b, log = TRUE)

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

h.log.r.sng <- function(theta, d0, d1, Omega.inv, beta, c, r.tilde, V.r.inv, nu = NULL) {

  p <- unlist(lapply(Omega.inv, function(x) {nrow(x)}))

  r <- d0*sin(theta) + d1*cos(theta) + r.tilde

  if (is.null(r.tilde) & is.null(V.r.inv)) {
    r.tilde <- rep(0, length(d0))
    V.r.inv <- rep(0, length(d0))
  }
  s <- abs(r)
  val <- -sum(c(atrans.mc(array(beta/s, dim = p), Omega.inv))*(beta/s))/2 -
    sum(log(s)) + (2*c - 1)*sum(log(s)) - sum(c*s^(2))
  if (is.null(nu)) {
    if (is.vector(V.r.inv)) {
      val <- val + sum((r - r.tilde)^2*V.r.inv/2)
    } else if (is.matrix(V.r.inv)) {
      val <- val + crossprod((r - r.tilde), crossprod(V.r.inv, (r - r.tilde)))/2
    }
  } else {
    if (is.vector(V.r.inv)) {
      val <- val + sum((nu + 1)*log(1 + (r - r.tilde)^2*V.r.inv/nu)/2)
    } else if (is.matrix(V.r.inv)) {
      val <- val + sum((nu + 1)*log(1 + crossprod((r - r.tilde), crossprod(V.r.inv, (r - r.tilde)))/nu)/2)
    }
  }
  return(val)

}

h.log.r.spb <- function(theta, d0, d1, Omega.inv, beta, c, r.tilde = NULL, V.r.inv = NULL, deltas, nu = NULL) {

  p <- unlist(lapply(Omega.inv, function(x) {nrow(x)}))
  if (is.null(r.tilde) & is.null(V.r.inv)) {
    r.tilde <- rep(0, length(d0))
    V.r.inv <- rep(0, length(d0))
  }
  r = d0*sin(theta) + d1*cos(theta) + r.tilde
  alpha <- c/2
  s <- abs(r)
  multip <- (((2*gamma(3/c))/gamma(1/c))^(alpha/(1 - alpha))*f.deltas(deltas = deltas, c = c))

  val <- -sum(c(atrans.mc(array(beta/s, dim = p), Omega.inv))*(beta/s))/2 - sum(log(s)) + ((1 + alpha)/(1 - alpha) - 1)*sum(log(s)) - sum(multip*s^(2*alpha/(1 - alpha)))
  if (is.null(nu)) {
    if (is.vector(V.r.inv)) {
      val <- val + sum((r - r.tilde)^2*V.r.inv/2)
    } else if (is.matrix(V.r.inv)) {
      val <- val + crossprod((r - r.tilde), crossprod(V.r.inv, (r - r.tilde)))/2
    }
  } else {
    if (is.vector(V.r.inv)) {
      val <- val + sum((nu + 1)*log(1 + (r - r.tilde)^2*V.r.inv/nu)/2)
    } else if (is.matrix(V.r.inv)) {
      val <- val + sum((nu + 1)*log(1 + crossprod((r - r.tilde), crossprod(V.r.inv, (r - r.tilde)))/nu)/2)
    }
  }
  return(val)



}

h.log.r.spn <- function(theta, d0, d1, Omega.inv, beta, Psi.inv, r.tilde, V.r.inv, nu = NULL) {
  s <- r <- d0*sin(theta) + d1*cos(theta) + r.tilde
  p <- unlist(lapply(Omega.inv, function(x) {nrow(x)}))
  if (is.null(r.tilde) & is.null(V.r.inv)) {
    r.tilde <- rep(0, length(d0))
    V.r.inv <- rep(0, length(d0))
  }
  val <- -sum(c(atrans.mc(array(beta/s, dim = p), Omega.inv))*(beta/s))/2 - sum(log(abs(s))) - sum(c(atrans.mc(array(s, dim = p), Psi.inv))*(s))/2

  if (is.null(nu)) {
    if (is.vector(V.r.inv)) {
      val <- val + sum((r - r.tilde)^2*V.r.inv/2)
    } else if (is.matrix(V.r.inv)) {
      val <- val + crossprod((r - r.tilde), crossprod(V.r.inv, (r - r.tilde)))/2
    }
  } else {
    if (is.vector(V.r.inv)) {
      val <- val + sum((nu + 1)*log(1 + (r - r.tilde)^2*V.r.inv/nu)/2)
    } else if (is.matrix(V.r.inv)) {
      val <- val + sum((nu + 1)*log(1 + crossprod((r - r.tilde), crossprod(V.r.inv, (r - r.tilde)))/nu)/2)
    }
  }
  return(val)

}

sample.d <- function(theta, delta, V.half, V.r.inv = NULL, nu = NULL) {

  if (is.null(nu)) {
    precisions <- 1
  } else {
    if (is.vector(V.half)) {
      precisions <- rgamma(length(delta), (nu + 1)/2, (nu + (delta/V.half)^2)/2)
    } else if (is.matrix(V.half)) {
      precisions <- rgamma(1, (nu + 1)/2, (nu + crossprod(delta, crossprod(V.r.inv, delta)))/2)
    }

  }
  V.half <- V.half/sqrt(precisions)

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
                         q = NULL, deltas = NULL, nu = NULL) {
  delta <- r - r.tilde
  d <- sample.d(V.half = V.r.half, V.r.inv = V.r.inv, theta = eta, delta = delta, nu = nu)
  if (prior == "sng") {
    ll.args <- list("d0" = d$d0,
                    "d1" = d$d1,
                    "Omega.inv" = Omega.inv,
                    "beta" = beta,
                    "c" = c,
                    "r.tilde" = r.tilde,
                    "V.r.inv" = V.r.inv,
                    "nu" = nu)
    ll.fun <- "h.log.r.sng"
  } else if (prior == "spn") {
    ll.args <- list("d0" = d$d0,
                    "d1" = d$d1,
                    "Omega.inv" = Omega.inv,
                    "beta" = beta,
                    "Psi.inv" = Psi.inv,
                    "r.tilde" = r.tilde,
                    "V.r.inv" = V.r.inv, "nu" = nu)
    ll.fun <- "h.log.r.spn"
  } else if (prior == "spb") {
    ll.args <- list("d0" = d$d0,
                    "d1" = d$d1,
                    "Omega.inv" = Omega.inv,
                    "beta" = beta,
                    "c" = c,
                    "deltas" = deltas,
                    "r.tilde" = r.tilde,
                    "V.r.inv" = V.r.inv, "nu" = nu)
    ll.fun <- "h.log.r.spb"
  }
  eta <- slice(x.tilde = eta, ll.fun = ll.fun,
               var.lim = c(0, 2*pi),
               ll.args = ll.args)
  delta <- d$d0*sin(eta) + d$d1*cos(eta)
  return(list("eta" = eta,
              "r" = r.tilde + delta))
}

sample.beta.theta <- function(X, U, y, V.half, beta.prev, theta, beta.tilde, Omega.inv, V.inv,
                              sig.sq, reg) {
  delta <- beta.prev - beta.tilde
  # cat("Sample d\n")
  d <- sample.d(V.half = V.half, V.inv = V.inv, theta = theta, delta = delta)
  # cat("Sample theta\n")
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
  # cat("Assemble beta\n")
  delta <- d$d0*sin(theta) + d$d1*cos(theta)
  return(list("theta" = theta,
              "beta" = beta.tilde + delta))
}

#' @export
sampler <- function(
  ### Data and regression type
  X, # Array of penalized covariates, covariance along second dimension is AR-1
  y, # Outcome
  reg = "linear", # Regression model for data
  U = NULL, # Matrix of unpenalized covariates
  ### Prior Choice for beta
  prior = "sno",
  c = 1,
  ### Prior Parameters and Likelihood Parameters
  Omega.half = NULL, # A dim(X) - 1 list of symmetric square roots of covariance matrices
  Psi.half = NULL,   # A dim(X) - 1 list of symmetric square roots of covariance matrices
  sig.sq = NULL,
  ### MCMC Parameters
  fix.beta = NULL, # Null if beta should not be fixed, a q + p vector otherwise
  num.samp = 100, # Number of samples to return
  burn.in = 0, # Number of burn-in samples to discard
  thin = 1, # Number of samples to thin by
  print.iter = TRUE, # Indicator for whether or not iteration counter should be printed
  max.iter = 1000, # Maximum number of outer iterations in coordinate descent for beta (and r if not max.iter.r not specified)
  max.iter.r = max.iter, # Maximum number of outer iterations in coordinate descent for r
  eps = 10^(-12), # Convergence threshold for coordinate descent beta (and r if eps.r not specified)
  eps.r = eps, # Convergence threshold for coordinate descent r
  diag.app = FALSE, # Whether or not a diagonal approximation to the covariance
  diag.app.r = diag.app, #
  from.prior = FALSE,
  do.svd = TRUE,
  slice.beta = TRUE,
  joint.beta = list(1:(prod(dim(X)[-1]) + ifelse(is.null(U), 1, ncol(U)))),
  use.previous = FALSE,
  use.previous.r = use.previous,
  max.inner = 100,
  max.inner.r = max.inner,
  sep.theta = list(1:(prod(dim(X)[-1]) + ifelse(is.null(U), 1, ncol(U)))),
  sep.eta = list(1:(prod(dim(X)[-1]))),
  V.inv = NULL,
  V.r.inv = NULL,
  z.tilde = NULL,
  r.tilde = NULL,
  r.start = NULL, # Starting value for MCMC for r
  z.start = NULL,
  gamma.start = NULL,
  nu = NULL, # t-distribution parameter for slice proposals for beta (and r if not specified)
  nu.r = nu, # t-distribution parameter for slice proposals for r
  ### Hyperparameters (if prior/likelihood parameters not specified)
  pr.rho.a = 1,
  pr.rho.b = 1,
  str = "uns", # Variance-covariance matrix type
  pr.Omega.V.inv = lapply(dim(X)[-1], diag),
  pr.Psi.V.inv = lapply(dim(X)[-1], diag),
  pr.Omega.df = lapply(dim(X)[-1], function(x) {x + 2}),
  pr.Psi.df = lapply(dim(X)[-1], function(x) {x + 2}),
  pr.sig.sq.shape = 3/2,
  pr.sig.sq.rate = 1/2) {

  # Just to be safe, "initalize" all variables whose entries may depend on other objects
  max.iter.r = max.iter.r
  eps.r = eps.r
  diag.app.r = diag.app.r
  joint.beta = joint.beta
  max.inner.r = max.inner.r
  sep.theta = sep.theta
  sep.eta = sep.eta
  nu.r = nu.r
  pr.Omega.V.inv = pr.Omega.V.inv
  pr.Psi.V.inv = pr.Psi.V.inv
  pr.Omega.df = pr.Omega.df
  pr.Psi.df = pr.Psi.df

  null.V.inv <- is.null(V.inv)
  null.V.r.inv <- is.null(V.r.inv)
  null.z.tilde <- is.null(z.tilde)
  null.r.tilde <- is.null(r.tilde)

  if (!diag.app) {
    # Can only use separate thetas if we are using independent proposals
    if (print.iter) {cat("Separate values of the slice variable for beta are only possible for diagonal covariance matrices\n")}
    sep.theta <- list(1:(prod(dim(X)[-1]) + ifelse(is.null(U), 1, ncol(U))))
  }
  if (!diag.app.r) {
    # Can only use separate etas if we are using independent proposals
    if (print.iter) {cat("Separate values of the slice variable for r are only possible for diagonal covariance matrices\n")}
    sep.eta <- list(1:(prod(dim(X)[-1])))
  }

  n <- length(y)
  p <- dim(X)[-1]
  X.arr <- X
  X <- t(apply(X, 1, "c"))

  # Set up indicators for null arguments, will be used to decide whether or not to resample
  if (is.null(Omega.half[[1]]) & str == "con" | str == "het") {
    Omega.half[[1]] <- diag(p[1], nrow = p[1], ncol = p[1])
    if (prior == "spn") {
      Psi.half[[1]] <- diag(p[1], nrow = p[1], ncol = p[1])
    }
  }
  null.Omega.half <- unlist(lapply(Omega.half, function(x) {is.null(x)}))

  # Set starting values
  null.sig.sq <- is.null(sig.sq) & reg == "linear"
  if (null.sig.sq | reg == "logit") {
    sig.sq <- 1
  } else {
    sig.sq <- sig.sq
  }

  if (null.Omega.half[1] & dim(X.arr)[2] > 1) {
    rho <- 0
    rho.psi <- 0
  }

  null.U <- is.null(U)
  if (max(null.Omega.half) == 1) {
    res.Omega <- vector("list", length(p))
    res.Sigma <- vector("list", length(p))
    if (prior == "spn") {
      res.Psi <- vector("list", length(p))
    }
    for (i in which(null.Omega.half)) {
      if (i != 1) {
        Omega.half[[i]] <- diag(p[i], nrow = p[i], ncol = p[i])
      } else {
        Omega.half[[i]] <- sym.sq.root(make.ar.mat(p = p[i], rho = rho, inv = FALSE))
      }
      if (i > 1) {
        res.Omega[[i]] <- array(dim = c(num.samp, p[i], p[i]))
        res.Sigma[[i]] <- array(dim = c(num.samp, p[i], p[i]))
      }
      if (prior == "spn") {
        if (i != 1) {
          Psi.half[[i]] <- diag(p[i])
        } else {
          Psi.half[[i]] <- sym.sq.root(make.ar.mat(p = p[i], rho = rho.psi, inv = FALSE))
        }
        if (i > 1) {
          res.Psi[[i]] <- array(dim = c(num.samp, p[i], p[i]))
        }
      }
    }
  }
  Omega <- lapply(Omega.half, function(x) {crossprod(x)})
  Omega.inv <- lapply(Omega, function(x) {ei.inv(x)})
  Omega.half.inv <- lapply(Omega.half, function(x) {ei.inv(x)})

  if (prior == "spn") {
    Psi <- lapply(Psi.half, function(x) {crossprod(x)})
    Psi.inv <- lapply(Psi, function(x) {ei.inv(x)})
    Psi.half.inv <- lapply(Psi.half, function(x) {ei.inv(x)})
  }

  # No intercept if U is null
  if (null.U) {
    U <- matrix(0, nrow = n, ncol = 1)
  }
  q <- ncol(U)

  # Results objects
  res.rho <- array(dim = c(num.samp, 1))
  res.rho.psi <- array(dim = c(num.samp, 1))
  res.sig.sq <- array(dim = c(num.samp, 1))
  res.B <- array(dim = c(num.samp, prod(p)))
  res.Z <- array(dim = c(num.samp, prod(p) + q))
  res.R <-res.S <- array(dim = c(num.samp, prod(p)))
  res.theta <- array(dim = c(num.samp, prod(p) + q))
  res.ome <- array(dim = c(num.samp, n))
  res.eta <- array(dim = c(num.samp, prod(p)))
  res.gamma <- array(dim = c(num.samp, q))
  res.D <- array(dim = c(num.samp, prod(p)))

  S <- array(1, dim = p)
  if (!is.null(fix.beta)) {
    gamma <- fix.beta[1:q]
    beta <- fix.beta[(q + 1):length(fix.beta)]
  }
  if (is.null(gamma.start)) {
    gamma <- rep(0, q)
  } else {
    gamma <- gamma.start
  }
  if (is.null(z.start)) {
    Z <- array(0, dim = p)
  }  else {
    Z <- array(z.start[(q + 1):length(z.start)], dim = p)
  }
  if (is.null(r.start)) {
    R <- array(1, dim = p)
  } else {
    R <- array(r.start, dim = p)
  }
  if (slice.beta) {
    theta <- numeric(prod(p) + q)
    for (i in 1:length(sep.theta)) {
      theta[sep.theta[[i]]] <- runif(1, 0, 2*pi)
    }
  } else {
    if (reg == "logit" | reg == "nb") {
      ome <- rep(1, n)
      omeD <- Matrix::Matrix(0, nrow = n, ncol = n, sparse = TRUE)
      diag(omeD) <- ome
      if (reg == "logit") {
        offset <- rep(1/2, n)
      }
    } else {
      offset <- rep(0, n)
    }
  }
  eta <- numeric(prod(p))
  for (i in 1:length(sep.eta)) {
    eta[sep.eta[[i]]] <- runif(1, 0, 2*pi)
  }

  X.arr.s <- X.arr
  X.arr.s <- array(c(crossprod(apply(X.arr, 1, "c"), diag(c(S)))), dim = c(n, p))
  for (l in 1:length(p)) {
    X.arr.s <- amprod.mc(X.arr.s, Omega.half[[l]], l + 1)
  }
  W <- t(apply(X.arr.s, 1, "c"))
  UW <- cbind(U, W)

  penC <- c(rep(0, ncol(U)), rep(1, ncol(W)))

  B <- atrans.mc(Z, Omega.half)

  for (i in 1:(burn.in + thin*num.samp)) {
    if (print.iter) {cat("i=", i, "\n")}

    if (is.null(fix.beta)) {
      sample.beta <- c(gamma, c(Z))
      if (print.iter) {cat("Sample Beta\n")}
      if (slice.beta) {


        if (null.z.tilde & null.V.inv) {
          # Set mean for proposal distribution
          if (from.prior) {
            z.tilde <- rep(0, ncol(UW))
          } else if ((i == 1 & prior == "sno" & max(null.Omega.half[-1]) == 0 & (!null.Omega.half[1]) & reg == "linear" & !null.sig.sq) |
                     (i >= 1 & (prior != "sno" | max(null.Omega.half[-1]) != 0 | (null.Omega.half[1]) | reg != "linear" | null.sig.sq))) {

            if (use.previous & i > 1) {
              z.start <- c(z.tilde[1:q], atrans.mc(B.tilde/S, Omega.half.inv))
              # This can work poorly if elements of S are super small...
            } else {
              z.start <- rep(0, ncol(UW))
            }

            if (reg == "logit") {
              if (print.iter) {cat("Get Mode\n")}

              z.tilde <- coord.desc.logit(y = y, X = UW, Omega.inv = penC,
                                          print.iter = FALSE, max.iter = max.iter, eps = eps,
                                          start.beta = z.start,
                                          joint.beta = joint.beta, max.inner = max.inner)$beta
            } else if (reg == "linear") {
              if (print.iter) {cat("Get Mode\n")}
              z.tilde <- coord.desc.lin(y = y, X = UW, sig.sq = sig.sq, Omega.inv = penC,
                                        print.iter = FALSE, max.iter = max.iter, eps = eps,
                                        start.beta = z.start)$beta
            }
            B.tilde <- atrans.mc(array(z.tilde[(q + 1):length(z.tilde)], dim = p), Omega.half)*S
          }
          if (from.prior) {
            V.half <- c(rep(10^(12), q), rep(1, length(z.tilde) - q))
            V.inv <- 1/V.half^2
          } else {
            if  (print.iter) {cat("Get Pieces for Covariance Matrix\n")}
            UWz.tilde <- crossprod(t(UW), z.tilde)[, 1]
            if (reg == "logit") {

              if (ncol(U) != 1 | !diag.app) {
                AA <- diag(exp(UWz.tilde)/(1 + exp(UWz.tilde))^2)
                AAU <- crossprod(sqrt(AA), U)
              } else {
                AAU <- matrix(sqrt(exp(UWz.tilde)/(1 + exp(UWz.tilde))^2)*U, nrow = nrow(U), ncol = 1)
              }

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
                if (ncol(U) != 1) {
                  V.half <- sqrt(c(diag(solve(crossprod(AAU)/sig.sq)), rep(1, prod(p))))
                } else {
                  V.half <- sqrt(c(sum(AAU^2/sig.sq), rep(1, prod(p))))
                }

              }

              V.inv <- 1/V.half^2

            } else {
              if (print.iter) {cat("Get Covariance Matrix\n")}
              UWtBB <- crossprod(UW, BB)/sig.sq
              V.inv <- UWtBB + diag(penC)
              V.half <- sym.sq.root.inv(V.inv)
            }

          }
        } else {
          z.tilde <- z.tilde
          if (is.matrix(V.inv)) {
            V.half <- sym.sq.root.inv(V.inv)
          } else if (is.vector(V.inv)) {
            V.half <- sqrt(1/V.inv)
          }

        }
        if  (print.iter) {cat("Sample Beta and Theta\n")}
        sample <- sample.beta.theta(X = W, U = U, y = y, V.half = V.half, beta.prev = c(gamma, c(Z)),
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

            } else if ((i == 1 & reg == "linear" & max(null.Omega.half) == 0 & prior == "sno" & length(joint.beta) == 1) |
                       !(reg == "linear") | !(max(null.Omega.half) == 0) | !(prior == "sno") | length(joint.beta) > 1) {
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



    if (i == 1) {
      deltas <- runif(prod(p), 0, pi)
    } else {
      if (prior == "spb") {
        if (print.iter) {cat("Sample Delta\n")}
        for (ii in 1:length(deltas)) {
          alpha <- c/2
          xi <- (2*gamma(3/(2*alpha))*c(S)[ii]^2/gamma(1/(2*alpha)))^(alpha/(1 - alpha))
          deltas[ii] <- slice(x.tilde = deltas[ii], ll.fun = "g.delta",
                              var.lim = c(0, pi), ll.args = list("c" = c, "xi" = xi))
        }
      }
    }

    if (prior == "sng" | prior == "spn" | prior == "spb") {

      if (null.r.tilde & null.V.r.inv) {
        if (i == 1 | !use.previous.r) {
          start.r <- rep(1, prod(p))
        } else {
          start.r <- r.tilde
        }

        Omega.inv.diag <- rep(diag(Omega.inv[[1]]), times = prod(p[-1]))
        for (k in 2:length(p)) {
          Omega.inv.diag <- rep(diag(Omega.inv[[k]]), each = prod(p[1:(k - 1)]))*Omega.inv.diag
        }

        # r.tilde <- rep(0, prod(p))
        # V.r.inv <- c(B)^2*Omega.inv.diag

        if (prior == "sng") {


          kappa1 <- 1 - 2*c
          kappa2 <- -c
          kappa3 <- -1

          # if (i == 1 | (is.null(fix.beta) | max(null.Omega.half) == 1)) {
            if (print.iter) {cat("Set Sampling Values for R\n")}
            r.tilde <- coord.desc.r(Omega.inv = Omega.inv,
                                    beta = c(B), c = c, eps = eps.r, max.iter = max.iter.r,
                                    print.iter = FALSE, max.inner = max.inner.r,
                                    start.r = start.r,
                                    prior = prior)$r

          # }


          if (diag.app.r) {
            r.tilde <- abs(r.tilde) # Sign doesn't matter
            V.r.inv <- numeric(length(r.tilde))
            for (jj in 1:length(V.r.inv)) {

              Omega.inv.jj <- get.kron.row(jj, Omega = Omega.inv)
              alpha1 <- -c(B)[jj]^2*Omega.inv.jj[jj]/2
              alpha2 <- -sum(c(B)[jj]*Omega.inv.jj[-jj]*abs(r.tilde[-jj])*c(B)[-jj])

              if (r.tilde[jj] != 0) {
                hess <- sn.ll.dd(rj = r.tilde[jj], alpha1 = alpha1, alpha2 = alpha2, kappa1 = kappa1, kappa2 = kappa2, kappa3 = kappa3)
              } else {
                hess <- 2*alpha1
              }
              V.r.inv[jj] <- -1*hess
            }

            V.r.inv[V.r.inv < 10^(-10) | is.nan(V.r.inv)] <- 10^(-10) # Make sure we don't have problems with infinity/0 (should be careful about this)
            V.r.inv[is.infinite(V.r.inv)] <- 10^(10)
            V.r.half <- sqrt(1/V.r.inv)
          } else {
            # print(range(abs(r.tilde)))
            r.tilde <- abs(r.tilde + ifelse(r.tilde == 0, 10^(-12), r.tilde)) # Sign doesn't matter
            if (i == 1 | max(null.Omega.half) == 1) {
              O.i <- do.call("%x%", Omega.inv[length(Omega.inv):1])
            }
            V.r.inv <- -1*((-1/2)*(O.i)*tcrossprod(c(B)/r.tilde^2)*(matrix(1,
                                                                           nrow = length(r.tilde),
                                                                           ncol = length(r.tilde)) +
                                                                      5*diag(length(r.tilde))) +
                             diag(1/(r.tilde)^2 + -(2*c - 1)*1/(r.tilde)^2 - 2*c))

            V.r.half <- sym.sq.root.inv(V.r.inv)
            V.r.inv <- ei.inv(crossprod(V.r.half)) # Make sure it's pos-semi-def
          }


          if (print.iter) {cat("Sample R\n")}

          sample <- sample.r.eta(r = c(R), Omega.inv = Omega.inv, beta = c(B), c = c,
                                 eta = eta,
                                 r.tilde = r.tilde, V.r.inv = V.r.inv,
                                 V.r.half = V.r.half,
                                 prior = prior, nu = nu.r)

        } else if (prior == "spn") {

          if (i == 1 | max(null.Omega.half)==1) {
            Psi.inv.diag <- rep(diag(Psi.inv[[1]]), times = prod(p[-1]))
            for (k in 2:length(p)) {
              Psi.inv.diag <- rep(diag(Psi.inv[[k]]), each = prod(p[1:(k - 1)]))*Psi.inv.diag
            }
          }

          # if (i == 1 | (is.null(fix.beta) | max(null.Omega.half) == 1)) {
            if (print.iter) {cat("Set Sampling Values for R\n")}
            r.tilde <- coord.desc.r(Omega.inv = Omega.inv,
                                    beta = c(B), c = c, eps = eps.r, max.iter = max.iter.r,
                                    print.iter = FALSE, max.inner = max.inner.r,
                                  start.r = start.r,
                                  prior = prior, Psi.inv = Psi.inv)$r
            # r.tilde[abs(r.tilde) > 10^(14)] <- sign(r.tilde[abs(r.tilde) > 10^(14)])*10^(14)
            # print(range(abs(r.tilde)))
          # }

          if (diag.app.r) {
          V.r.inv <- numeric(length(r.tilde))
          for (jj in 1:length(V.r.inv)) {

            Omega.inv.jj <- get.kron.row(jj, Omega = Omega.inv)
            Psi.inv.jj <- get.kron.row(jj, Omega = Psi.inv)
            alpha1 <- -c(B)[jj]^2*Omega.inv.jj[jj]/2
            alpha2 <- -sum(c(B)[jj]*Omega.inv.jj[-jj]*abs(r.tilde[-jj])*c(B)[-jj])

            kappa1 <- 1/2
            kappa2 <- -Psi.inv.jj[jj]/2
            kappa3 <- -sum(Psi.inv.jj[-jj]*r.tilde[-jj]/2)

            if (r.tilde[jj] != 0) {
              hess <- sn.ll.dd(rj = r.tilde[jj], alpha1 = alpha1, alpha2 = alpha2, kappa1 = kappa1, kappa2 = kappa2, kappa3 = kappa3)
            } else {
              hess <- 2*alpha1
            }
            V.r.inv[jj] <- -1*hess
          }


          V.r.inv[V.r.inv < 10^(-10)] <- 10^(-10) # Make sure we don't have problems with infinity/0 (should be careful about this)
          V.r.inv[is.infinite(V.r.inv) | is.nan(V.r.inv)] <- 10^(10)
          V.r.half <- sqrt(1/V.r.inv)
          } else {
            # print(range(abs(r.tilde)))
            r.tilde <- abs(r.tilde + ifelse(r.tilde == 0, 10^(-12), r.tilde)) # Sign doesn't matter
            if (i == 1 | max(null.Omega.half) == 1) {
              O.i <- do.call("%x%", Omega.inv[length(Omega.inv):1])
              P.i <- do.call("%x%", Psi.inv[length(Psi.inv):1])
            }
            V.r.inv <- -1*((-1/2)*(O.i)*tcrossprod(c(B)/r.tilde^2)*(matrix(1,
                                                                           nrow = length(r.tilde),
                                                                           ncol = length(r.tilde)) +
                                                                      5*diag(length(r.tilde))) +
                             - diag(-1/abs(r.tilde)^2) - P.i/2)

            V.r.half <- sym.sq.root.inv(V.r.inv)
            V.r.inv <- ei.inv(crossprod(V.r.half)) # Make sure it's pos-semi-def
          }


          if (print.iter) {cat("Sample R\n")}
          sample <- sample.r.eta(r = c(R), Omega.inv = Omega.inv, beta = c(B), Psi.inv = Psi.inv, eta = eta,
                                 r.tilde = r.tilde, V.r.inv = V.r.inv, V.r.half = V.r.half, c = c, prior = prior,
                                 nu = nu.r)
        } else if (prior == "spb") {


          if (print.iter) {cat("Set Sampling Values for R\n")} # Have to do every iteration

          kappa1 <- 2*(c/2)/(c/2 - 1)
          kappa2 <- -(((2*gamma(3/c))/gamma(1/c))^(c/2/(1 - c/2))*f.deltas(deltas = deltas, c = c))
          kappa3 <- c/2/(c/2 - 1)

          start.r <- rep(1, prod(p))
          r.tilde <- coord.desc.r(Omega.inv = Omega.inv,
                                  beta = c(B), c = c, eps = eps.r, max.iter = max.iter.r, print.iter = FALSE,
                                  max.inner = max.inner.r, start.r = start.r,
                                  deltas = deltas, prior = prior)$r
          r.tilde <- abs(r.tilde) # Sign doesn't matter

          if (diag.app.r) {
          V.r.inv <- numeric(length(r.tilde))
          for (jj in 1:length(V.r.inv)) {

            Omega.inv.jj <- get.kron.row(jj, Omega = Omega.inv)
            alpha1 <- -c(B)[jj]^2*Omega.inv.jj[jj]/2
            alpha2 <- -sum(c(B)[jj]*Omega.inv.jj[-jj]*abs(r.tilde[-jj])*c(B)[-jj])

            if (r.tilde[jj] != 0) {
              hess <- sn.ll.dd(rj = r.tilde[jj], alpha1 = alpha1, alpha2 = alpha2, kappa1 = kappa1, kappa2 = kappa2[jj], kappa3 = kappa3)
            } else {
              hess <- 2*alpha1
            }
            V.r.inv[jj] <- -1*hess
          }

          V.r.inv[V.r.inv < 10^(-10)] <- 10^(-10) # Make sure we don't have problems with infinity/0 (should be careful about this)
          V.r.inv[is.infinite(V.r.inv)] <- 10^(10)
          V.r.half <- sqrt(1/V.r.inv)
          } else {

            r.tilde <- abs(r.tilde + ifelse(r.tilde == 0, 10^(-12), r.tilde)) # Sign doesn't matter
            if (i == 1 | max(null.Omega.half) == 1) {
              O.i <- do.call("%x%", Omega.inv[length(Omega.inv):1])
            }
            multip <- (((2*gamma(3/c))/gamma(1/c))^((c/2)/(1 - (c/2)))*f.deltas(deltas = deltas, c = c))
            V.r.inv <- -1*((-1/2)*(O.i)*tcrossprod(c(B)/r.tilde^2)*(matrix(1,
                                                                           nrow = length(r.tilde),
                                                                           ncol = length(r.tilde)) +
                                                                      5*diag(length(r.tilde))) +
                             - diag(-1/abs(r.tilde)^2 + multip*((2*(c/2)/(1 - (c/2))))*((2*(c/2)/(1 - (c/2))) - 1)*r.tilde^(2*(c/2)/(1 - (c/2)) - 2)) + diag(-((1 + (c/2))/(1 - (c/2)) - 1)/r.tilde^2))

            V.r.half <- sym.sq.root.inv(V.r.inv)
            V.r.inv <- ei.inv(crossprod(V.r.half)) # Make sure it's pos-semi-def

          }

          if (print.iter) {cat("Sample R\n")}

          sample <- sample.r.eta(r = c(R), Omega.inv = Omega.inv, beta = c(B), deltas = deltas, eta = eta,
                                 r.tilde = r.tilde, V.r.inv = V.r.inv, V.r.half =
                                   V.r.half, c = c, prior = prior, nu = nu.r)
        }
      } else {
        if (is.matrix(V.r.inv)) {
          V.r.half <- sym.sq.root.inv(V.r.inv)
        } else {
          V.r.half <- sqrt(1/V.r.inv)
        }
        sample <- sample.r.eta(r = c(R), Omega.inv = Omega.inv, beta = c(B), deltas = deltas, eta = eta,
                               r.tilde = r.tilde, V.r.inv = V.r.inv, V.r.half =
                                 V.r.half, c = c, prior = prior, Psi.inv = Psi.inv, nu = nu.r)
      }
      r <- sample$r
      # print(r)
      R <- array(r, p)
      eta <- sample$eta
      if (prior %in% c("sng", "spb")) {
        S <- array(abs(r), dim = p)
      } else if (prior  == "spn") {
        S <- array(r, dim = p)
      } else {
        S <- array(1, dim = p)
      }
      if (is.null(fix.beta)) {
        Z <- atrans.mc(B/S, Omega.half.inv)
        B <- S*atrans.mc(Z, Omega.half)
      }

    }

    for (k in which(null.Omega.half)) {
      Covar <- B/S
      for (l in 1:length(p)) {
        if (l != k) {
          Covar <- amprod.mc(Covar, Omega.half.inv[[l]], l)
        }
      }
      if (k == 1 & null.Omega.half[1]) {
        if (print.iter) {cat("Sample rho\n")}
        rho <- slice(x.tilde = rho, ll.fun = "cond.rho.log", var.lim = c(-1, 1)*(1 - 10^(-7)),
                     ll.args = list("B" = Covar,
                                    "pr.a" = pr.rho.a,
                                    "pr.b" = pr.rho.b,
                                    "j" = k))
        # if (print.iter) {
        #   cat("rho = ", rho, "\n")
        # }
        Omega.inv[[k]] <- make.ar.mat(p = p[k], rho = rho, inv = TRUE)

      } else {
        if (print.iter) {cat("Sample Omega.inv ", k, " \n")}
        Covar <- apply(Covar, k, "c")
        Omega.inv[[k]] <- samp.Omega.inv(Beta = Covar, str = str,
                                         pr.V.inv = pr.Omega.V.inv[[k]],
                                         pr.df = pr.Omega.df[[k]])
      }
      Omega[[k]] <- ei.inv(Omega.inv[[k]])
      if (i > burn.in & (i - burn.in)%%thin == 0 & k != 1) {
        res.Omega[[k]][(i - burn.in)/thin, , ] <- Omega[[k]]
      }

      Omega.half[[k]] <- sym.sq.root(Omega[[k]])
      Omega.half.inv[[k]] <- sym.sq.root(Omega.inv[[k]])

      if (prior == "spn") {

        Covar <- S
        for (l in 1:length(p)) {
          if (l != k) {
            Covar <- amprod.mc(Covar, Psi.half.inv[[l]], l)
          }
        }

        if (k == 1 & null.Omega.half[1]) {
          if (print.iter) {cat("Sample rho psi\n")}
          rho.psi <- slice(x.tilde = rho.psi, ll.fun = "cond.rho.log", var.lim = c(-1, 1),
                           ll.args = list("B" = Covar,
                                          "pr.a" = pr.rho.a,
                                          "pr.b" = pr.rho.b,
                                          "j" = k))
          Psi.inv[[k]] <- make.ar.mat(p = p[k], rho = rho.psi, inv = TRUE)
        } else {
          Covar <- apply(Covar, k, "c")
          Psi.inv[[k]] <- samp.Omega.inv(Beta = Covar, str = str,
                                         pr.V.inv = pr.Psi.V.inv[[k]],
                                         pr.df = pr.Psi.df[[k]])
        }
        Psi[[k]] <- ei.inv(Psi.inv[[k]])
        Psi.half[[k]] <- sym.sq.root(Psi[[k]])
        Psi.half.inv[[k]] <- sym.sq.root(Psi.inv[[k]])

        if (i > burn.in & (i - burn.in)%%thin == 0 & k != 1) {
          res.Psi[[k]][(i - burn.in)/thin, , ] <- Psi[[k]]
        }
      }

      if (reg == "linear" & null.sig.sq) {
        resid <- y - crossprod(t(U), gamma) - crossprod(t(X), c(B))
        sig.sq <- 1/rgamma(1, shape = pr.sig.sq.shape + length(resid)/2, rate = pr.sig.sq.rate + sum(resid^2)/2)
      }

      if (i > burn.in & (i - burn.in)%%thin == 0 & k != 1) {
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


    if (is.null(fix.beta) & max(null.Omega.half) == 1) {
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

    if (i > burn.in & (i - burn.in)%%thin == 0) {

      res.B[(i - burn.in)/thin, ] <- c(B)
      res.Z[(i - burn.in)/thin, ] <- c(gamma, c(B))
      res.gamma[(i - burn.in)/thin, ] <- gamma
      if (slice.beta) {
        res.theta[(i - burn.in)/thin, ] <- theta
      }
      res.S[(i - burn.in)/thin, ] <- c(S)
      res.R[(i - burn.in)/thin, ] <- c(R)
      if (prior %in% c("sng", "spb", "spn")) {
        res.eta[(i - burn.in)/thin, ] <- eta
      }
      if (prior == "spb") {
        res.D[(i - burn.in)/thin, ] <- c(deltas)
      }
      if (null.Omega.half[1]) {
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
                   "Ss" = res.S, "Zs" = res.Z, "Rs" = res.R)
  if (slice.beta) {
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
em.est <- function(print.iter.em = TRUE,
                   max.iter.em = NULL,
                   eps.em = 10^(-3),
                   ### Data and regression type
                   X, # Array of penalized covariates, covariance along second dimension is AR-1
                   y, # Outcome
                   reg = "linear", # Regression model for data
                   U = NULL, # Matrix of unpenalized covariates
                   ### Prior Choice for beta
                   prior = "sno",
                   c = 1,
                   ### Prior Parameters and Likelihood Parameters
                   Omega.half = NULL, # A dim(X) - 1 list of symmetric square roots of covariance matrices
                   Psi.half = NULL,   # A dim(X) - 1 list of symmetric square roots of covariance matrices
                   sig.sq = NULL,
                   ### MCMC Parameters
                   num.samp = 100, # Number of samples to return
                   burn.in = 0, # Number of burn-in samples to discard
                   thin = 1, # Number of samples to thin by
                   print.iter = TRUE, # Indicator for whether or not iteration counter should be printed
                   max.iter = 1000, # Maximum number of outer iterations in coordinate descent for beta (and r if not max.iter.r not specified)
                   max.iter.r = max.iter, # Maximum number of outer iterations in coordinate descent for r
                   eps = 10^(-12), # Convergence threshold for coordinate descent beta (and r if eps.r not specified)
                   eps.r = eps, # Convergence threshold for coordinate descent r
                   diag.app = FALSE, # Whether or not a diagonal approximation to the covariance
                   diag.app.r = diag.app, #
                   from.prior = FALSE,
                   do.svd = TRUE,
                   slice.beta = TRUE,
                   joint.beta = list(1:(prod(dim(X)[-1]) + ifelse(is.null(U), 1, ncol(U)))),
                   use.previous = FALSE,
                   use.previous.r = use.previous,
                   max.inner = 100,
                   max.inner.r = max.inner,
                   sep.theta = list(1:(prod(dim(X)[-1]) + ifelse(is.null(U), 1, ncol(U)))),
                   sep.eta = list(1:(prod(dim(X)[-1]))),
                   V.inv = NULL,
                   V.r.inv = NULL,
                   z.tilde = NULL,
                   r.tilde = NULL,
                   r.start = NULL, # Starting value for MCMC for r
                   z.start = NULL,
                   gamma.start = NULL,
                   nu = NULL, # t-distribution parameter for slice proposals for beta (and r if not specified)
                   nu.r = nu, # t-distribution parameter for slice proposals for r
                   ### Hyperparameters (if prior/likelihood parameters not specified)
                   pr.rho.a = 1,
                   pr.rho.b = 1,
                   str = "uns", # Variance-covariance matrix type
                   pr.Omega.V.inv = lapply(dim(X)[-1], diag),
                   pr.Psi.V.inv = lapply(dim(X)[-1], diag),
                   pr.Omega.df = lapply(dim(X)[-1], function(x) {x + 2}),
                   pr.Psi.df = lapply(dim(X)[-1], function(x) {x + 2}),
                   pr.sig.sq.shape = 3/2,
                   pr.sig.sq.rate = 1/2) {

  W <- t(apply(X, 1, "c"))

  if (is.null(U)) {
    UW <- cbind(rep(0, nrow(W)), W)
  } else {
    UW <- cbind(U, W)
  }

  penC <- matrix(0, nrow = ncol(UW), ncol = ncol(UW))

  Omega.inv <- lapply(Omega.half, function(x) {ei.inv(crossprod(x))})
  O.i <- do.call("%x%", Omega.inv[length(Omega.inv):1])

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
  if (length(burn.in) == 1) {
    burn.in <- rep(burn.in, max.iter.em + 1)
  }
  if (length(thin) == 1) {
    thin <- rep(thin, max.iter.em + 1)
  }
  if (print.iter.em) {cat("Set Starting Value\n")}
  # Get initial values
  samples <- sampler(### Data and regression type
    X = X, # Array of penalized covariates, covariance along second dimension is AR-1
    y = y, # Outcome
    reg = reg, # Regression model for data
    U = U, # Matrix of unpenalized covariates
    ### Prior Choice for beta
    prior = "sno",
    c = c,
    ### Prior Parameters and Likelihood Parameters
    Omega.half = Omega.half, # A dim(X) - 1 list of symmetric square roots of covariance matrices
    Psi.half = Psi.half,   # A dim(X) - 1 list of symmetric square roots of covariance matrices
    sig.sq = sig.sq,
    ### MCMC Parameters
    fix.beta = NULL, # Null if beta should not be fixed, a q + p vector otherwise
    num.samp = num.samp[1], # Number of samples to return
    burn.in = burn.in[1], # Number of burn-in samples to discard
    thin = thin[1], # Number of samples to thin by
    print.iter = print.iter, # Indicator for whether or not iteration counter should be printed
    max.iter = max.iter, # Maximum number of outer iterations in coordinate descent for beta (and r if not max.iter.r not specified)
    max.iter.r = max.iter.r, # Maximum number of outer iterations in coordinate descent for r
    eps = eps, # Convergence threshold for coordinate descent beta (and r if eps.r not specified)
    eps.r = eps.r, # Convergence threshold for coordinate descent r
    diag.app = diag.app, # Whether or not a diagonal approximation to the covariance
    diag.app.r = diag.app.r, #
    from.prior = from.prior,
    do.svd = do.svd,
    slice.beta = slice.beta,
    joint.beta = joint.beta,
    use.previous = use.previous,
    use.previous.r = use.previous.r,
    max.inner = max.inner,
    max.inner.r = max.inner.r,
    sep.theta = sep.theta,
    sep.eta = sep.eta,
    V.inv = V.inv,
    V.r.inv = V.r.inv,
    z.tilde = z.tilde,
    r.tilde = r.tilde,
    r.start = r.start, # Starting value for MCMC for r
    z.start = z.start,
    gamma.start = gamma.start,
    nu = nu, # t-distribution parameter for slice proposals for beta (and r if not specified)
    nu.r = nu.r, # t-distribution parameter for slice proposals for r
    ### Hyperparameters (if prior/likelihood parameters not specified)
    pr.rho.a = pr.rho.a,
    pr.rho.b = pr.rho.b,
    str = "uns", # Variance-covariance matrix type
    pr.Omega.V.inv = pr.Omega.V.inv,
    pr.Psi.V.inv = pr.Psi.V.inv,
    pr.Omega.df = pr.Omega.df,
    pr.Psi.df = pr.Psi.df,
    pr.sig.sq.shape = pr.sig.sq.shape,
    pr.sig.sq.rate = pr.sig.sq.rate)
  post.mean <- c(colMeans(samples$gammas), colMeans(samples$Bs))
  post.median <- c(apply(samples$gammas, 2, median), apply(samples$Bs, 2, median))

  beta.fix <- post.median

  betas <- matrix(nrow = max.iter.em, ncol = prod(dim(X)[-1]) + ifelse(is.null(U), 0, ncol(U)))
  ess <- matrix(nrow = max.iter.em, ncol = prod(dim(X)[-1]))

  for (i in 1:max.iter.em) {

    if (print.iter.em) {cat("EM Iteration: ", i, "\n")}

    beta.fix[beta.fix == 0] <- rnorm(sum(beta.fix == 0))

    samples <- sampler(### Data and regression type
      X = X, # Array of penalized covariates, covariance along second dimension is AR-1
      y = y, # Outcome
      reg = reg, # Regression model for data
      U = U, # Matrix of unpenalized covariates
      ### Prior Choice for beta
      prior = prior,
      c = c,
      ### Prior Parameters and Likelihood Parameters
      Omega.half = Omega.half, # A dim(X) - 1 list of symmetric square roots of covariance matrices
      Psi.half = Psi.half,   # A dim(X) - 1 list of symmetric square roots of covariance matrices
      sig.sq = sig.sq,
      ### MCMC Parameters
      fix.beta = beta.fix, # Null if beta should not be fixed, a q + p vector otherwise
      num.samp = num.samp[i + 1], # Number of samples to return
      burn.in = burn.in[i + 1], # Number of burn-in samples to discard
      thin = thin[i + 1], # Number of samples to thin by
      print.iter = print.iter, # Indicator for whether or not iteration counter should be printed
      max.iter = max.iter, # Maximum number of outer iterations in coordinate descent for beta (and r if not max.iter.r not specified)
      max.iter.r = max.iter.r, # Maximum number of outer iterations in coordinate descent for r
      eps = eps, # Convergence threshold for coordinate descent beta (and r if eps.r not specified)
      eps.r = eps.r, # Convergence threshold for coordinate descent r
      diag.app = diag.app, # Whether or not a diagonal approximation to the covariance
      diag.app.r = diag.app.r, #
      from.prior = from.prior,
      do.svd = do.svd,
      slice.beta = slice.beta,
      joint.beta = joint.beta,
      use.previous = use.previous,
      use.previous.r = use.previous.r,
      max.inner = max.inner,
      max.inner.r = max.inner.r,
      sep.theta = sep.theta,
      sep.eta = sep.eta,
      V.inv = V.inv,
      V.r.inv = V.r.inv,
      z.tilde = z.tilde,
      r.tilde = r.tilde,
      r.start = r.start, # Starting value for MCMC for r
      z.start = z.start,
      gamma.start = gamma.start,
      nu = nu, # t-distribution parameter for slice proposals for beta (and r if not specified)
      nu.r = nu.r, # t-distribution parameter for slice proposals for r
      ### Hyperparameters (if prior/likelihood parameters not specified)
      pr.rho.a = pr.rho.a,
      pr.rho.b = pr.rho.b,
      str = "uns", # Variance-covariance matrix type
      pr.Omega.V.inv = pr.Omega.V.inv,
      pr.Psi.V.inv = pr.Psi.V.inv,
      pr.Omega.df = pr.Omega.df,
      pr.Psi.df = pr.Psi.df,
      pr.sig.sq.shape = pr.sig.sq.shape,
      pr.sig.sq.rate = pr.sig.sq.rate)
    if (prior != "spn") {
      inv.ss <- matrix(rowMeans(apply(samples$Ss, 1, function(x) {tcrossprod(1/abs(x))})), nrow = prod(p), ncol = prod(p))
    } else {
      inv.ss <- matrix(rowMeans(apply(samples$Ss, 1, function(x) {tcrossprod(1/x)})), nrow = prod(p), ncol = prod(p))
    }


    penC[2:nrow(penC), 2:ncol(penC)] <- (O.i*inv.ss)

    if (reg == "logit") {
      beta.fix <- coord.desc.logit(y = y, X = UW, Omega.inv = penC,
                                   print.iter = FALSE, max.iter = max.iter, eps = eps,
                                   max.inner = max.inner, joint.beta = joint.beta)$beta
    } else if (reg == "linear") {
      beta.fix <- coord.desc.lin(y = y, X = UW, Omega.inv = penC,
                                 print.iter = FALSE, max.iter = max.iter, eps = eps, sig.sq = sig.sq)$beta
    }

    if (is.null(U)) {
      betas[i, ] <- beta.fix[-1]
    } else {
      betas[i, ] <- beta.fix
    }
    ess[i, ] <- coda::effectiveSize(samples$Ss)

    if (i > 1) {
      if (max(abs(betas[i - 1, ] - betas[i, ])) < eps.em) {
        break
      }
    }
  }

  if (is.null(U)) {
    post.mean <- post.mean[-1]
    post.median <- post.median[-1]
  }

  return(list("post.mean" = post.mean, "post.med" = post.median, betas = betas[1:i, ],
              esss = ess[1:i, ]))

}

