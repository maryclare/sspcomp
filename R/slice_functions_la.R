#' @export
get.beta.blocks <- function(X, U = NULL, min.block.size = 25, no.eig = TRUE) {

  p <- dim(X)[-1]
  if (is.null(U)) {
    q <- 0
  } else {
    q <- dim(U)[-1]
  }

  if (min.block.size > prod(p) + q) {
    min.block.size <- prod(p) + q - 1
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
samp.Omega.inv <- function(Beta, pr.V.inv = diag(1, nrow = ncol(Beta), ncol = ncol(Beta)),
                           pr.df = ncol(Beta) + 2, str = "uns") {
  p <- ncol(Beta)
  if (str == "uns") {
    V.inv <- crossprod(Beta) + pr.V.inv
    df <- nrow(Beta) + pr.df
    V.half <- as.matrix(sym.sq.root.inv((V.inv + t(V.inv))/2))
    return(tcrossprod(crossprod(V.half, matrix(rnorm(p*df), nrow = p, ncol = df))))
    # return(matrix(rWishart(1, df, solve(V.inv))[, , 1], nrow = p, ncol = p))
  } else if (str == "het") {
    b <- apply(Beta, 2, function(x) {sum(x^2)})/(2) + diag(pr.V.inv)/2
    a <- rep(nrow(Beta), ncol(Beta))/2 + pr.df/2
    return(diag(rgamma(p, shape = a, rate = b), nrow = p, ncol = p))
  } else if (str == "con") {
    b <- sum(apply(Beta, 2, function(x) {sum(x^2)}))/(2) + sum(diag(pr.V.inv))/2
    # I'm a little worried about the code below if 'Beta' is a matrix wtih more than 1 column,
    # Should check. I think it works!
    a <- sum(rep(nrow(Beta), ncol(Beta)))/2 + p*pr.df/2
    return(rgamma(1, shape = a, rate = b)*diag(p, nrow = p, ncol = p))
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
    ARMat <- diag(1, nrow = p, ncol = p)
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
    ARMat <- diag(1, nrow = p, ncol = p)
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

cond.xi.log <- function(theta, B, pr.a, pr.b, j, W.j, lower.xi, upper.xi, W.ei,
                        sum.B.mat.sq = NULL, sum.BWB = NULL) {

  p <- dim(B)

  if (is.null(sum.B.mat.sq)) {
    sum.B.mat.sq <- sum(B^2)
  }
  if (is.null(sum.BWB)) {
    B.mat <- mat(B, j)
    sum.BWB <- sum(Matrix::diag(Matrix::crossprod(B.mat, Matrix::crossprod(W.j, B.mat))))
  }

  c1 <- prod(p[-j])*(1/2)*sum(log(1 - theta*W.ei$values))

  c2 <- -(sum.B.mat.sq - theta*sum.BWB)/2

  c3 <- dbeta((theta - lower.xi)/(upper.xi - lower.xi), pr.a, pr.b, log = TRUE)

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

  X.tilde <- tcrossprod(diag(sqrt(1 + d^2), nrow = length(d), ncol = length(d)), R)
  V.half[[length(V.half)]] <- sqrt(rep(1, prod(p)) - apply(R, 1, function(x) {sum(x^2*(d^2/(1 + d^2)))}))
  V.inv[[length(V.inv)]] <- sqrt(diag(crossprod(X.tilde)))
  X.tilde.arr <- array(c(tcrossprod(X.tilde, diag(V.half[[length(V.half)]], nrow = length(V.half[[length(V.half)]]), ncol = length(V.half[[length(V.half)]])))),
                       c(nrow(X.tilde), p))

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
        V.mat <- V.mat + (abs(min(V.mat.val)) + 10^(-8))*diag(1, nrow = p[k], ncol = p[k])
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
  A.eig <- eigen((A + Matrix::t(A))/2)
  Matrix::crossprod(Matrix::t(A.eig$vectors), Matrix::tcrossprod(Matrix(diag(sqrt(ifelse(A.eig$values > 0, A.eig$values, 0)),
                                              nrow = nrow(A), ncol = ncol(A))), A.eig$vectors))
}

#' @export
sym.sq.root.inv <- function(A) {
  A.eig <- eigen(A)
  Matrix::crossprod(Matrix::t(A.eig$vectors), Matrix::tcrossprod(Matrix(diag(sqrt(ifelse(A.eig$values > 0, 1/A.eig$values, 0)),
                                              nrow = nrow(A), ncol = ncol(A))), A.eig$vectors))
}

#' @export
ei.inv <- function(A) {
  A.eig <- eigen(A)
  Matrix::crossprod(Matrix::t(A.eig$vectors), Matrix::tcrossprod(Matrix(diag(ifelse(A.eig$values > 0, 1/A.eig$values, 0),
                                              nrow = nrow(A), ncol = ncol(A))), A.eig$vectors))
}

#' @export
amprod.mc <- function (A, M, k) {
  K <- length(dim(A))
  AM <- Matrix::crossprod(Matrix::t(M), mat(A, k))
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
  d <- sample.d(V.half = V.r.half, V.r.inv = V.r.inv,
                theta = eta, delta = delta, nu = nu)
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

sample.beta.theta <- function(X, U, y, V.half, beta.prev, theta, beta.tilde, Omega.inv,
                              V.inv,
                              sig.sq, reg) {
  delta <- beta.prev - beta.tilde
  # cat("Sample d\n")
  d <- sample.d(V.half = V.half,
                theta = theta, delta = delta)
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
  Omega.half = vector("list", length = length(dim(X)[-1])), # A dim(X) - 1 list of symmetric square roots of covariance matrices
  Psi.half = vector("list", length = length(dim(X)[-1])),   # A dim(X) - 1 list of symmetric square roots of covariance matrices
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
  use.previous = TRUE,
  use.previous.r = use.previous,
  max.inner = 1000,
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
  pr.rho.omega.a = 1,
  pr.rho.omega.b = 1,
  pr.xi.omega.a = 1,
  pr.xi.omega.b = 1,
  pr.rho.psi.a = 1,
  pr.rho.psi.b = 1,
  pr.xi.psi.a = 1,
  pr.xi.psi.b = 1,
  Neighbs = NULL, # Neighbor matrix if needed
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
  Omega.half = Omega.half
  Psi.half = Psi.half
  pr.Omega.V.inv = pr.Omega.V.inv
  pr.Psi.V.inv = pr.Psi.V.inv
  pr.Omega.df = pr.Omega.df
  pr.Psi.df = pr.Psi.df

  # Set up indicators for which things are null
  null.V.inv <- is.null(V.inv)
  null.V.r.inv <- is.null(V.r.inv)
  null.z.tilde <- is.null(z.tilde)
  null.r.tilde <- is.null(r.tilde)

  if (!is.null(Neighbs)) {
    Neighbs.ei <- eigen(Neighbs)
    lower.xi <- 1/Neighbs.ei$values[length(Neighbs.ei$values)]
    upper.xi <- 1/Neighbs.ei$values[1]
  }

  n <- length(y)
  p <- dim(X)[-1]

  null.U <- is.null(U)
  # No intercept if U is null
  if (null.U) {
    U <- matrix(0, nrow = n, ncol = 1)
  }
  q <- ncol(U)

  if (!diag.app) {
    block.theta <- matrix(0, nrow = prod(p) + q, ncol = prod(p) + q)
    for (i in 1:length(sep.theta)) {
      block.theta[sep.theta[[i]], sep.theta[[i]]] <- 1
    }
    if (!is.null(V.inv) & min(block.theta) == 0) {
      V.inv <- block.theta*(V.inv)
    }
  }
  if (!diag.app.r) {
    # Can only use separate etas if we are using independent proposals
    block.eta <- matrix(0, nrow = prod(p), ncol = prod(p))
    for (i in 1:length(sep.eta)) {
      block.eta[sep.eta[[i]], sep.eta[[i]]] <- 1
    }
    if (!is.null(V.r.inv) & min(block.eta) == 0) {
      V.r.inv <- block.eta*V.r.inv
    }
  }

  X.arr <- X
  X <- t(apply(X, 1, "c"))

  if (is.null(Psi.half) & prior == "spn") {
    Psi.half <- vector("list", length = length(Omega.half))
  }

  # Set up indicators for null arguments, will be used to decide whether or not to resample
  if (is.null(Omega.half[[1]]) & str == "con" | str == "het") {
    Omega.half[[1]] <- diag(1, nrow = p[1], ncol = p[1])
    if (prior == "spn") {
      Psi.half[[1]] <- diag(1, nrow = p[1], ncol = p[1])
    }
  }
  null.Omega.half <- unlist(lapply(Omega.half, function(x) {is.null(x)}))
  if (prior == "spn") {
    null.Psi.half <- unlist(lapply(Psi.half, function(x) {is.null(x)}))
  }

  # Set starting values
  null.sig.sq <- is.null(sig.sq) & reg == "linear"
  if (null.sig.sq | reg == "logit") {
    sig.sq <- 1
  } else {
    sig.sq <- sig.sq
  }

  # Set starting values for rho if Omega.half[[1]] is NULL and is not 1x1
  if (null.Omega.half[1] & dim(X.arr)[2] > 1) {
    if (is.null(Neighbs)) {
      rho <- 0
    } else {
      rho <- (lower.xi + upper.xi)/2
    }
  }
  if (prior == "spn" & dim(X.arr)[2] > 1) {
    if (null.Psi.half[1]) {
      if (is.null(Neighbs)) {
        rho.psi <- 0
      } else {
        rho.psi <- (lower.xi + upper.xi)/2
      }
    }
  }

  if (prior == "spn") {
    record.Sigma <- max(null.Psi.half) == 1 | max(null.Omega.half) == 1
    record.which <- which(null.Psi.half | null.Omega.half)
  } else {
    record.Sigma <- max(null.Omega.half) == 1
    record.which <- which(null.Omega.half)
  }
  if (record.Sigma) {
    res.Sigma <- vector("list", length(p))
    for (i in record.which) {
      if (i > 1) {
        res.Sigma[[i]] <- array(dim = c(num.samp, p[i], p[i]))
      }
    }
  }
  if (max(null.Omega.half) == 1) {
    res.Omega <- vector("list", length(p))
    for (i in which(null.Omega.half)) {
      if (i != 1) {
        Omega.half[[i]] <- diag(1, nrow = p[i], ncol = p[i])
      } else {
        if (is.null(Neighbs)) {
          Omega.half[[i]] <- sym.sq.root(make.ar.mat(p = p[i], rho = rho, inv = FALSE))
        } else {
          Omega.half[[i]] <- Matrix::crossprod(Matrix::t(Neighbs.ei$vectors), Matrix::tcrossprod(diag(sqrt(ifelse(1 - rho*Neighbs.ei$values > 0, 1/(1 - rho*Neighbs.ei$values), 0)),
                                                                              nrow = nrow(Neighbs), ncol = ncol(Neighbs)), Neighbs.ei$vectors))
        }
      }
      if (i > 1) {
        res.Omega[[i]] <- array(dim = c(num.samp, p[i], p[i]))
      }
    }
  }
  if (prior == "spn") {
    if (max(null.Psi.half) == 1) {
      res.Psi <- vector("list", length(p))
      for (i in which(null.Psi.half)) {
        if (i != 1) {
          Psi.half[[i]] <- diag(1, nrow = p[i], ncol = p[i])
        } else {
          if (is.null(Neighbs)) {
            Psi.half[[i]] <- sym.sq.root(make.ar.mat(p = p[i], rho = rho.psi, inv = FALSE))
          } else {
            Psi.half[[i]] <- Matrix::crossprod(Matrix::t(Neighbs.ei$vectors), Matrix::tcrossprod(diag(sqrt(ifelse(1 - rho.psi*Neighbs.ei$values > 0, 1/(1 - rho.psi*Neighbs.ei$values), 0)),
                                                                                                      nrow = nrow(Neighbs), ncol = ncol(Neighbs)), Neighbs.ei$vectors))
          }
        }
        if (i > 1) {
          res.Psi[[i]] <- array(dim = c(num.samp, p[i], p[i]))
        }
      }
    }
  }
  Omega <- lapply(Omega.half, function(x) {Matrix::crossprod(x)})
  Omega.inv <- lapply(Omega, function(x) {ei.inv(x)})
  Omega.half.inv <- lapply(Omega.half, function(x) {ei.inv(x)})

  if (prior == "spn") {
    Psi <- lapply(Psi.half, function(x) {Matrix::crossprod(x)})
    Psi.inv <- lapply(Psi, function(x) {ei.inv(x)})
    Psi.half.inv <- lapply(Psi.half, function(x) {ei.inv(x)})
  }



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
  S <- R
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
  X.arr.s <- array(c(crossprod(apply(X.arr, 1, "c"),
                               diag(c(S), nrow = prod(p), ncol = prod(p)))), dim = c(n, p))
  for (l in 1:length(p)) {
    X.arr.s <- amprod.mc(X.arr.s, Omega.half[[l]], l + 1)
  }
  W <- t(apply(X.arr.s, 1, "c"))
  UW <- cbind(U, W)

  penC <- c(rep(0, ncol(U)), rep(1, ncol(W)))

  if (is.null(fix.beta)) {
    B <- atrans.mc(Z, Omega.half)
  } else {
    B <- array(c(fix.beta), dim = p)
  }

  for (i in 1:(burn.in + thin*num.samp)) {
    if (print.iter) {cat("i=", i, "\n")}

    sample.beta <- c(gamma, c(Z))

    if (is.null(fix.beta)) {

      if (prior == "spn") {
        samp.vars <- c("z", "s")
      } else {
        samp.vars <- c("z")
      }
      for (sv in samp.vars) {
        if (print.iter & sv == "z") {cat("Sample Z\n")}
        if (print.iter & sv == "s") {cat("Sample S\n")}


      if (sv == "z") {

            X.arr.s <- array(c(crossprod(apply(X.arr, 1, "c"),
                                         diag(c(S),
                                              nrow = prod(p),
                                              ncol = prod(p)))),
                             dim = c(n, p))
            for (l in 1:length(p)) {
              X.arr.s <- amprod.mc(X.arr.s, Omega.half[[l]], l + 1)
            }
            W <- t(apply(X.arr.s, 1, "c"))
            UW <- cbind(U, W)

          } else {

            X.arr.s <- X.arr
            X.arr.s <- array(c(crossprod(apply(X.arr, 1, "c"),
                                         diag(c(atrans.mc(Z, Omega.half)),
                                              nrow = prod(p), ncol = prod(p)))), dim = c(n, p))
            for (l in 1:length(p)) {
              X.arr.s <- amprod.mc(X.arr.s, Psi.half[[l]], l + 1)
            }
            W <- t(apply(X.arr.s, 1, "c"))
            UW <- cbind(U, W)

          }


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
                  AA <- diag(exp(UWz.tilde)/(1 + exp(UWz.tilde))^2, nrow = length(UWz.tilde), ncol = length(UWz.tilde))
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
                V.inv <- UWtBB + diag(penC, nrow = length(penC), ncol = length(penC))
                if (min(block.theta) == 0) {
                  V.inv <- solve(block.theta*solve(V.inv))
                }
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
                                      Omega.inv = lapply(p, function(x) {diag(1, nrow = x, ncol = x)}),
                                      V.inv = V.inv,
                                      sig.sq = sig.sq, reg = reg)
          sample.beta <- sample$beta
          theta <- sample$theta
        } else {
          if (reg == "logit") {
            if (prior != "spn" | (prior == "spn" & sv == "z")) {
              if (print.iter) {cat("Sample logit auxiliary variables\n")}
              ome <- BayesLogit::rpg(n, offset*2, crossprod(t(UW), c(gamma, c(Z))))
              diag(omeD) <- ome
            }
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
                                                             Omega.inv = diag(penC[(q + 1):nrow(UWty)],
                                                                              nrow = length(penC[(q + 1):nrow(UWty)]),
                                                                              ncol = length(penC[(q + 1):nrow(UWty)])), sig.sq = sig.sq)
              } else {
                sample.beta <- samp.beta(XtX = UWtUW, Xty = UWty,
                                         Omega.inv = diag(penC, nrow = length(penC), ncol = length(penC)), sig.sq = sig.sq)
              }


            } else {
              if (min(U) == 0 & max(U) == 0) {

                sample.beta[(q + 1):nrow(UWty)] <- samp.beta(XtX = UWtUW[(q + 1):nrow(UWtUW), (q + 1):ncol(UWtUW)],
                                                             Xty = UWty[(q + 1):nrow(UWty)],
                                                             Omega.inv = diag(penC[(q + 1):nrow(UWty)],
                                                                              nrow = length(penC[(q + 1):nrow(UWty)]),
                                                                              ncol = length(penC[(q + 1):nrow(UWty)])), sig.sq = 1)
              } else {
                sample.beta <- samp.beta(XtX = UWtUW, Xty = UWty,
                                         Omega.inv = diag(penC, nrow = length(penC), ncol = length(penC)), sig.sq = 1)
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
                                                Omega.inv = diag(penC[block], nrow = length(penC[block]), ncol = length(penC[block])), sig.sq = sig.sq)
              } else {
                sample.beta[block] <- samp.beta(XtX = UWtUW, Xty = UWty - UWtNUWZ,
                                                Omega.inv = diag(penC[block], nrow = length(penC[block]), ncol = length(penC[block])), sig.sq = 1)
              }
            }

          }

        }

        gamma <- sample.beta[1:q]
        if (prior != "spn" | (prior == "spn" & sv == "z")) {
          # if (prior != "spn") {
          Z <- array(sample.beta[(q + 1):length(sample.beta)], dim = p)
          # } else {
          #   Z <- array(1, prod(p))
          # }
        } else {
          S <- atrans.mc(array(sample.beta[(q + 1):length(sample.beta)], p),
                         Psi.half)
          # S <- array(1, prod(p))
        }

        S <- array(S, p)
        Z <- array(Z, p)
        B <- S*atrans.mc(Z, Omega.half)

      }
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

    if (prior == "sng" | (prior == "spn" & !is.null(fix.beta)) | prior == "spb") {

      if (null.r.tilde) {

        # Don't have to reset r.tilde every iteration if
        # - prior %in% c("sng", "spn)
        # - !is.null(fix.beta)
        # - max(null.Psi.half) == 0
        # - max(null.Omega.half) == 0
        # - !null.sig.sq
        once <- !is.null(fix.beta) & prior %in% c("sng", "spn") & max(null.Omega.half) == 0 & !null.sig.sq
        if (prior == "spn") {
          once <- once*(max(null.Psi.half) == 0)
        }

        if ((i == 1 & once) | !once) {

        if (print.iter) {cat("Set Sampling Values for R\n")}
        if (use.previous.r | i == 1) {
          start.r <- rep(1, prod(p))
        } else {
          start.r <- r.tilde
        }
        r.tilde <- coord.desc.r(Omega.inv = Omega.inv,
                                beta = c(B), c = c, eps = eps.r, max.iter = max.iter.r,
                                print.iter = FALSE, max.inner = max.inner.r,
                                start.r = start.r,
                                prior = prior, deltas = deltas, Psi.inv = Psi.inv)$r
        r.tilde[r.tilde == 0] <- 10^(-12)
        r.tilde[abs(r.tilde) > 10^(12)] <- sign(r.tilde[abs(r.tilde) > 10^(12)])*10^(12)
      }
      if (null.V.r.inv) {



        if (diag.app.r) {
          V.r.inv <- numeric(length(r.tilde))
          for (jj in 1:length(V.r.inv)) {

            Omega.inv.jj <- get.kron.row(jj, Omega = Omega.inv)
            alpha1 <- -c(B)[jj]^2*Omega.inv.jj[jj]/2

            if (r.tilde[jj] != 0) {
              hess <- kappa.ll.dd(s.j = r.tilde[jj], kappa = get.kappa(prior = prior, Omega.inv = Omega.inv, beta = c(B),
                                                                       r = r.tilde, c = c, deltas = deltas,
                                                                       Psi.inv = Psi.inv, jj))
            } else { # I don't think we ever use this because we reset values of r.tilde equal to zero
              hess <- 2*alpha1
            }
            V.r.inv[jj] <- -1*hess
          }

          V.r.inv[V.r.inv < 10^(-10) | is.nan(V.r.inv)] <- 10^(-10) # Make sure we don't have problems with infinity/0 (should be careful about this)
          V.r.inv[is.infinite(V.r.inv)] <- 10^(10)
          V.r.half <- sqrt(1/V.r.inv)
        } else {
          if (i == 1 | max(null.Omega.half) == 1) {
            O.i <- do.call("%x%", Omega.inv[length(Omega.inv):1])
            if (prior == "spn") {
              P.i <- do.call("%x%", Psi.inv[length(Psi.inv):1])
            }
          }
          V.r.inv.1 <- (-1/2)*(O.i)*Matrix::tcrossprod(c(B)/r.tilde^2)*(matrix(1,
                                                                       nrow = length(r.tilde),
                                                                       ncol = length(r.tilde)) +
                                                                  5*diag(1, nrow = length(r.tilde), ncol = length(r.tilde)))
          if (prior == "spn") {
            V.r.inv.2 <- -P.i/2
          } else {

            V.r.inv.2 <- diag(0, nrow = length(r.tilde), ncol = length(r.tilde))
            for (jj in 1:length(r.tilde)) {
              kappa <- get.kappa(prior = prior, Omega.inv = Omega.inv, beta = c(B),
                                 r = r.tilde, c = c, deltas = deltas,
                                 Psi.inv = Psi.inv, jj)
              kappa[1:2] <- 0
              V.r.inv.2[jj, jj] <- kappa.ll.dd(s.j = r.tilde[jj], kappa = kappa)
            }
          }
          V.r.inv <- -1*(V.r.inv.1 + V.r.inv.2)
          if (!is.null(V.r.inv) & min(block.eta) == 0) {
            V.r.inv <- block.eta*(V.r.inv)
          }
          V.r.half <- sym.sq.root.inv(V.r.inv)
          V.r.inv <- ei.inv(Matrix::crossprod(V.r.half)) # Make sure it's pos-semi-def
        }
      }
      }

      if (print.iter) {cat("Sample R\n")}

      sample <- sample.r.eta(r = c(R), Omega.inv = Omega.inv, beta = c(B), c = c,
                             eta = eta,
                             r.tilde = r.tilde, V.r.inv = V.r.inv,
                             V.r.half = V.r.half,
                             prior = prior, nu = nu.r, deltas = deltas, Psi.inv = Psi.inv)


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
        if (is.null(Neighbs)) {
        rho <- slice(x.tilde = rho, ll.fun = "cond.rho.log",
                     var.lim = c(-1, 1)*(1 - 10^(-7)), # Avoid boundaries
                     ll.args = list("B" = Covar,
                                    "pr.a" = pr.rho.omega.a,
                                    "pr.b" = pr.rho.omega.b,
                                    "j" = k))

        Omega.inv[[k]] <- make.ar.mat(p = p[k], rho = rho, inv = TRUE)

        Omega[[k]] <- make.ar.mat(p = p[k], rho = rho, inv = FALSE)
        Omega.half[[k]] <- sym.sq.root(Omega[[k]])
        Omega.half.inv[[k]] <- sym.sq.root(Omega.inv[[k]])

        } else {
          sum.B.mat.sq <- sum(c(Covar^2))
          B.mat <- mat(Covar, k)
          sum.BWB <- sum(Matrix::diag(Matrix::crossprod(B.mat, Matrix::crossprod(Neighbs, B.mat))))
          rho <- slice(x.tilde = rho, ll.fun = "cond.xi.log",
                       var.lim = c(lower.xi, upper.xi)*(1 - 10^(-7)), # Avoid boundaries
                       ll.args = list("B" = Covar,
                                      "pr.a" = pr.xi.omega.a,
                                      "pr.b" = pr.xi.omega.b,
                                      "j" = k,
                                      "lower.xi" = lower.xi,
                                      "upper.xi" = upper.xi,
                                      "W.j" = Neighbs,
                                      "W.ei" = Neighbs.ei
                                      ,
                                      "sum.B.mat.sq" = sum.B.mat.sq,
                                      "sum.BWB" = sum.BWB
                                      ))
          if (print.iter) {cat("Compute Omega Inv\n")}
          Omega.inv[[k]] <- diag(1, nrow = p[k], ncol = p[k]) - rho*Neighbs
          if (print.iter) {cat("Compute Omega Half\n")}
          Omega.half[[k]] <-  Matrix::crossprod(Matrix::tcrossprod(Matrix(diag(sqrt(sqrt(ifelse(1 - rho*Neighbs.ei$values > 0, 1/(1 - rho*Neighbs.ei$values), 0))),
                                                                               nrow = nrow(Neighbs), ncol = ncol(Neighbs))), Neighbs.ei$vectors))

          # I don't save this
          # if (print.iter) {cat("Compute Omega\n")}
          # Omega[[k]] <- Matrix::crossprod(Omega.half[[k]])

          if (print.iter) {cat("Compute Omega Half Inv\n")}
          Omega.half.inv[[k]] <- Matrix::crossprod(Matrix::tcrossprod(Matrix(diag(sqrt(sqrt(ifelse(1 - rho*Neighbs.ei$values > 0, 1 - rho*Neighbs.ei$values, 0))),
                                                                                                          nrow = nrow(Neighbs), ncol = ncol(Neighbs))), Neighbs.ei$vectors))
        }

      } else {

        if (print.iter) {cat("Sample Omega.inv ", k, " \n")}
        Covar <- matrix(apply(Covar, k, "c"), nrow = prod(p[-k]), ncol = p[k])
        Omega.inv[[k]] <- samp.Omega.inv(Beta = Covar, str = str,
                                         pr.V.inv = pr.Omega.V.inv[[k]],
                                         pr.df = pr.Omega.df[[k]])
        Omega[[k]] <- ei.inv(Omega.inv[[k]])
        Omega.half[[k]] <- sym.sq.root(Omega[[k]])
        Omega.half.inv[[k]] <- sym.sq.root(Omega.inv[[k]])
      }


      if (i > burn.in & (i - burn.in)%%thin == 0 & k != 1) {
        res.Omega[[k]][(i - burn.in)/thin, , ] <- as.matrix(Omega[[k]])
      }

      Z <- atrans.mc(B/S, Omega.half.inv)
      B <- S*atrans.mc(Z, Omega.half)

    }
    if (prior == "spn") {

      Psi.half.inv.old <- Psi.half.inv

      for (k in which(null.Psi.half)) {

        Covar <- S
        for (l in 1:length(p)) {
          if (l != k) {
            Covar <- amprod.mc(Covar, Psi.half.inv[[l]], l)
          }
        }
        if (k == 1 & null.Psi.half[1]) {
          if (print.iter) {cat("Sample rho.psi\n")}
          if (is.null(Neighbs)) {
          rho.psi <- slice(x.tilde = rho.psi, ll.fun = "cond.rho.log", var.lim = c(-1, 1)*(1 - 10^(-7)),
                           ll.args = list("B" = Covar,
                                          "pr.a" = pr.rho.psi.a,
                                          "pr.b" = pr.rho.psi.b,
                                          "j" = k))

          Psi.inv[[k]] <- make.ar.mat(p = p[k], rho = rho.psi, inv = TRUE)
          Psi[[k]] <- make.ar.mat(p = p[k], rho = rho.psi, inv = FALSE)

          Psi.half[[k]] <- sym.sq.root(Psi[[k]])
          Psi.half.inv[[k]] <- sym.sq.root(Psi.inv[[k]])

          } else {
            sum.B.mat.sq <- sum(c(Covar^2))
            B.mat <- mat(Covar, k)
            sum.BWB <- sum(Matrix::diag(Matrix::crossprod(B.mat, Matrix::crossprod(Neighbs, B.mat))))
            rho.psi <- slice(x.tilde = rho.psi, ll.fun = "cond.xi.log",
                             var.lim = c(lower.xi, upper.xi)*(1 - 10^(-7)),
                         ll.args = list("B" = Covar,
                                        "pr.a" = pr.xi.psi.a,
                                        "pr.b" = pr.xi.psi.b,
                                        "j" = k,
                                        "lower.xi" = lower.xi,
                                        "upper.xi" = upper.xi,
                                        "W.j" = Neighbs,
                                        "W.ei" = Neighbs.ei
                                        ,
                                        "sum.B.mat.sq" = sum.B.mat.sq,
                                        "sum.BWB" = sum.BWB
                                        ))
            Psi.inv[[k]] <- diag(1, nrow = p[k], ncol = p[k]) - rho.psi*Neighbs

            Psi.half[[k]] <-  Matrix::crossprod(Matrix::tcrossprod(Matrix(diag(sqrt(sqrt(ifelse(1 - rho.psi*Neighbs.ei$values > 0, 1/(1 - rho.psi*Neighbs.ei$values), 0))),
                                                                                                         nrow = nrow(Neighbs), ncol = ncol(Neighbs))), Neighbs.ei$vectors))

            Psi.half.inv[[k]] <- Matrix::crossprod(Matrix::tcrossprod(Matrix(diag(sqrt(sqrt(ifelse(1 - rho.psi*Neighbs.ei$values > 0, 1 - rho.psi*Neighbs.ei$values, 0))),
                                                                                                            nrow = nrow(Neighbs), ncol = ncol(Neighbs))), Neighbs.ei$vectors))
          }
        } else {
          if (print.iter) {cat("Sample Psi.inv ", k, " \n")}
          Covar <- matrix(apply(Covar, k, "c"), nrow = prod(p[-k]), ncol = p[k])
          Psi.inv[[k]] <- samp.Omega.inv(Beta = Covar, str = str,
                                         pr.V.inv = pr.Psi.V.inv[[k]],
                                         pr.df = pr.Psi.df[[k]])

          Psi[[k]] <- ei.inv(Psi.inv[[k]])
          Psi.half[[k]] <- sym.sq.root(Psi[[k]])
          Psi.half.inv[[k]] <- sym.sq.root(Psi.inv[[k]])
        }
        if (i > burn.in & (i - burn.in)%%thin == 0 & k != 1) {
          res.Psi[[k]][(i - burn.in)/thin, , ] <- as.matrix(Psi[[k]])
        }


        S <- atrans.mc(atrans.mc(S, Psi.half.inv.old), Psi.half)
        Z <- atrans.mc(B/S, Omega.half.inv)
        B <- S*atrans.mc(Z, Omega.half)

      }
  }

    for (k in record.which) {
      if (i > burn.in & (i - burn.in)%%thin == 0 & k != 1) {
        if (prior == "sno") {
          res.Sigma[[k]][(i - burn.in)/thin, , ] <- as.matrix(Omega[[k]])
        } else if (prior == "sng") {
          es <- c^(-1/2)*gamma(c + 1/2)/gamma(c)
          els <- (1 - es^2)*diag(1, nrow = p[k], ncol = p[k]) + es^2*matrix(1, nrow = p[k], ncol = p[k])
          res.Sigma[[k]][(i - burn.in)/thin, , ] <- as.matrix(Omega[[k]])*els
        } else if (prior == "spb") {
          es <- sqrt(pi/2)*gamma(2/c)/(sqrt(gamma(1/c)*gamma(3/c)))
          els <- (1 - es^2)*diag(1, nrow = p[k], ncol = p[k]) + es^2*matrix(1, nrow = p[k], ncol = p[k])
          res.Sigma[[k]][(i - burn.in)/thin, , ] <- as.matrix(Omega[[k]])*els
        } else if (prior == "spn") {
          res.Sigma[[k]][(i - burn.in)/thin, , ] <- as.matrix(Omega[[k]])*as.matrix(Psi[[k]])
        }
      }
    }

    if (reg == "linear" & null.sig.sq) {
      resid <- y - Matrix::crossprod(t(U), gamma) - Matrix::crossprod(t(X), c(B))
      sig.sq <- 1/rgamma(1, shape = pr.sig.sq.shape + length(resid)/2, rate = pr.sig.sq.rate + sum(resid^2)/2)
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
      if (prior == "logit" & !slice.beta) {
        res.ome[(i - burn.in)/thin, ] <- ome
      }
      if (prior %in% c("sng", "spb", "spn")) {
        res.eta[(i - burn.in)/thin, ] <- eta
      }
      if (prior == "spb") {
        res.D[(i - burn.in)/thin, ] <- c(deltas)
      }
      if (null.Omega.half[1]) {
        res.rho[(i - burn.in)/thin] <- rho
      }
      if (prior == "spn") {
        if (null.Psi.half[1]) {
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
    res.list[["rhos"]] <- res.rho
  }
  if (record.Sigma) {
    res.list[["Sigmas"]] <- res.Sigma
  }
  if (prior == "spn") {
    if (max(null.Psi.half) == 1) {
      res.list[["Psis"]] <- res.Psi
      res.list[["rho.psis"]] <- res.rho.psi
    }
  }
  if (reg == "linear" & null.sig.sq) {
    res.list[["sig.sqs"]] <- res.sig.sq
  }

  return(res.list)
}



#' @export
em.est <- function(max.iter.em = NULL,
                   print.iter.em = TRUE,
                   eps.em = NULL,
                   beta.start = NULL,
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
                   Psi.half = NULL,
                   sig.sq = NULL,
                   ### MCMC Parameters
                   num.samp = 100, # Number of samples to return
                   burn.in = 0, # Number of burn-in samples to discard
                   thin = 1, # Number of samples to thin by
                   print.iter = TRUE,
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
                   use.previous = TRUE,
                   use.previous.r = use.previous,
                   max.inner = 1000,
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
                   pr.rho.omega.a = 1,
                   pr.rho.omega.b = 1,
                   pr.xi.omega.a = 1,
                   pr.xi.omega.b = 1,
                   pr.rho.psi.a = 1,
                   pr.rho.psi.b = 1,
                   pr.xi.psi.a = 1,
                   pr.xi.psi.b = 1,
                   str = "uns", # Variance-covariance matrix type
                   pr.Omega.V.inv = lapply(dim(X)[-1], diag),
                   pr.Psi.V.inv = lapply(dim(X)[-1], diag),
                   pr.Omega.df = lapply(dim(X)[-1], function(x) {x + 2}),
                   pr.Psi.df = lapply(dim(X)[-1], function(x) {x + 2}),
                   pr.sig.sq.shape = 3/2,
                   pr.sig.sq.rate = 1/2) {

  W <- t(apply(X, 1, "c"))

  U.orig <- U
  if (is.null(U)) {
    U <- matrix(0, nrow = nrow(W), ncol = 1)
    UW <- cbind(U, W)
  } else {
    UW <- cbind(U, W)
  }

  q <- dim(U)[2]
  p <- dim(X)[-1]
  n <- dim(X)[1]

  Omega.inv <- lapply(Omega.half, function(x) {ei.inv(Matrix::crossprod(x))})
  O.i <- do.call("%x%", Omega.inv[length(Omega.inv):1])

  penC <- matrix(0, nrow = ncol(UW), ncol = ncol(UW))
  betas.em <- matrix(NA, nrow = max.iter.em + 1, ncol = ncol(UW))
  es.i <- Matrix::tcrossprod(rep(1, prod(p)))
  ess <- rep(NA, max.iter.em)

  for (k in 1:max.iter.em) {
    if (print.iter.em) {cat(toupper(prior), " ", c, " EM Iteration=", k, "\n")}
    penC[(q + 1):nrow(penC), (q + 1):ncol(penC)] <- as.matrix((O.i*es.i))
    if ((k > 1 & !is.null(beta.start)) | (is.null(beta.start))) {
    if (reg == "linear") {
      betas.em[k, ] <- coord.desc.lin(y = y, X = UW, sig.sq = sig.sq, Omega.inv = penC,
                                          print.iter = FALSE,
                                      max.iter = max.iter, eps = eps)$beta
    } else if (reg == "logit") {
      betas.em[k, ] <- coord.desc.logit(y = y, X = UW, Omega.inv = penC,
                                            print.iter = FALSE,
                                        max.iter = max.iter, eps = eps,
                                        max.inner = max.inner)$beta
    }
    if (k > 1) {
      diff <- mean((betas.em[k, ] - betas.em[k - 1, ])^2)
      if (!is.null(U.orig)) {
        zero.beta <- min(abs(betas.em[k, ])) == 0
      } else {
        zero.beta <- min(abs(betas.em[k, -1])) == 0
      }
      if (print.iter.em) {cat("Diff=", diff, "\n")}
      if (diff < eps.em | zero.beta) {
        break
      }
    }
    } else {
      if (is.null(U.orig)) {
        betas.em[k, ] <- c(0, beta.start)
      } else {
        betas.em[k, ] <- beta.start
      }
    }
    beta.fix <- betas.em[k, ]
    if (is.null(U.orig)) {
      beta.fix <- beta.fix[-1]
    }
    s.s <- sampler(### Data and regression type
      X = X, # Array of penalized covariates, covariance along second dimension is AR-1
      y = y, # Outcome
      reg = "linear", # Regression model for data
      U = U.orig, # Matrix of unpenalized covariates
      ### Prior Choice for beta
      prior = prior,
      c = c,
      ### Prior Parameters and Likelihood Parameters
      Omega.half = Omega.half, # A dim(X) - 1 list of symmetric square roots of covariance matrices
      Psi.half = Psi.half,
      sig.sq = sig.sq,
      ### MCMC Parameters
      fix.beta = beta.fix, # Null if beta should not be fixed, a q + p vector otherwise
      num.samp = num.samp, # Number of samples to return
      burn.in = burn.in, # Number of burn-in samples to discard
      thin = thin, # Number of samples to thin by
      print.iter = print.iter,
      max.iter.r = max.iter.r, # Maximum number of outer iterations in coordinate descent for r
      eps.r = eps.r, # Convergence threshold for coordinate descent r
      diag.app.r = diag.app.r, #
      use.previous.r = use.previous.r,
      max.inner.r = max.inner.r,
      sep.eta = sep.eta,
      V.r.inv = V.r.inv,
      r.tilde = r.start,
      r.start = r.start, # Starting value for MCMC for r
      nu.r = nu.r, # t-distribution parameter for slice proposals for r
      pr.rho.omega.a = pr.rho.omega.a,
      pr.rho.omega.b = pr.rho.omega.b,
      pr.xi.omega.a = pr.xi.omega.a,
      pr.xi.omega.b = pr.xi.omega.b,
      pr.rho.psi.a = pr.rho.psi.a,
      pr.rho.psi.b = pr.rho.psi.b,
      pr.xi.psi.a = pr.xi.psi.a,
      pr.xi.psi.b = pr.xi.psi.b,
      str = str,
      pr.Omega.V.inv = pr.Omega.V.inv,
      pr.Psi.V.inv = pr.Psi.V.inv,
      pr.Omega.df = pr.Omega.df,
      pr.Psi.df = pr.Psi.df,
      pr.sig.sq.shape = pr.sig.sq.shape,
      pr.sig.sq.rate = pr.sig.sq.rate)$Ss

    es.is <- t(apply(1/s.s, 1, tcrossprod))
    es.i <- matrix(colMeans(es.is, na.rm = TRUE), nrow = prod(p), ncol = prod(p))
    ess[k] <- min(coda::effectiveSize(s.s))
    if (print.iter.em) {cat("ESS=", ess[k], "\n")}

  }

  if (k == max.iter.em) {
  if (reg == "linear") {
    betas.em[k + 1, ] <- coord.desc.lin(y = y, X = UW, sig.sq = sig.sq, Omega.inv = penC,
                                        print.iter = FALSE,
                                        max.iter = max.iter, eps = eps)$beta
  } else if (reg == "logit") {
    betas.em[k + 1, ] <- coord.desc.logit(y = y, X = UW, Omega.inv = penC,
                                        print.iter = FALSE,
                                        max.iter = max.iter, eps = eps,
                                        max.inner = max.inner)$beta
  }
  } else {
    betas.em <- betas.em[1:k, ]
    ess <- ess[1:k]
  }

  if (is.null(U.orig)) {
    betas.em <- betas.em[, -1]
  }

  return(list("betas" = betas.em, "esss" = ess))

}

