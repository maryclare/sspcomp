# Function for setting kappa values
get.kappa <- function(prior, Omega.inv, beta, r, c = NULL, deltas = NULL,
                      Psi.inv = NULL, j) {

  Omega.inv.j <- get.kron.row(j, Omega = Omega.inv)
  kappa <- numeric(6)
  kappa[1] <- -beta[j]^2*Omega.inv.j[j]/2
  if (prior == "spn") {
    Psi.inv.j <- get.kron.row(j, Omega = Psi.inv)
    kappa[2] <- -beta[j]*sum(Omega.inv.j[-j]*beta[-j]/r[-j])
    kappa[3] <- -1/2
    kappa[4] <- -sum(Psi.inv.j[-j]*r[-j])
    kappa[5] <- -Psi.inv.j[j]/2
    kappa[6] <- 2
  } else if (prior == "sng") {
    kappa[2] <- -beta[j]*sum(Omega.inv.j[-j]*beta[-j]/abs(r[-j]))
    kappa[3] <- c - 1
    kappa[4] <- 0
    kappa[5] <- -c
    kappa[6] <- 2
  } else if (prior == "spb") {
    alpha <- c/2
    kappa[2] <- -beta[j]*sum(Omega.inv.j[-j]*beta[-j]/abs(r[-j]))
    kappa[3] <- (1 + alpha)/(2*(1 - alpha)) - 1
    kappa[4] <- 0
    kappa[5] <- -(((2*gamma(3/c))/gamma(1/c))^(alpha/(1 - alpha))*f.deltas(deltas = deltas[j], c = c))
    kappa[6] <- 2*(alpha/(1 - alpha))
  }

  return(kappa)

}

kappa.ll <- function(s.j, kappa) {
  kappa[1]*s.j^(-2) + kappa[2]*s.j^(-1) + kappa[3]*log(s.j^2) + kappa[4]*s.j + kappa[5]*s.j^kappa[6]
}

kappa.ll.dd <- function(s.j, kappa) {
  6*kappa[1]*s.j^(-4) + 2*kappa[2]*s.j^(-3) - 2*kappa[3]*s.j^(-2) + kappa[5]*kappa[6]*(kappa[6] - 1)*s.j^(kappa[6] - 2)
}

# Function for relating kappa values to polynomial coefficients
get.poly.coef <- function(kappa) {

  max.pow <- max(3, kappa[6] + 2)
  coef <- numeric(max.pow + 1)
  coef[1] <- -2*kappa[1]
  coef[2] <- -kappa[2]
  coef[3] <- 2*kappa[3]
  coef[4] <- kappa[4]
  coef[kappa[6] + 2 + 1] <- kappa[5]*kappa[6]

  return(coef)

}

solve.kappa <- function(kappa, prior) {
  poly <- polynom::polynomial(get.poly.coef(kappa))
  sol <- solve(poly)
  sol <- Re(sol[Im(sol) == 0])
  if (prior == "sng" | prior == "spb") {
    sol <- sol[sol >= 0]
  }
  # cat("sol=", sol, "\n ")
  if (length(sol) > 1) {
    sol.val <- numeric(length(sol))
    for (k in 1:length(sol.val)) {
      sol.val[k] <- kappa.ll(s.j = sol[k], kappa = kappa)
    }
    sol <- sol[which(sol.val == max(sol.val, na.rm = TRUE))]
  }
  return(sol[1])
}

#' @export
coord.desc.r <- function(Omega.inv, beta,
                         c = NULL, eps = 10^(-12),
                         max.iter = 1000, print.iter = TRUE,
                         start.r = rep(1, length(beta)),
                         max.inner = Inf, deltas = NULL, prior,
                         Psi.inv = NULL) {

  r <- start.r
  p <- length(r)
  objs <- rep(NA, max.iter)
  for (i in 1:max.iter) {
    for (j in 1:length(r)) {

      r.old <- r

      kappa <- get.kappa(prior = prior,
                         Omega.inv = Omega.inv,
                         beta = beta, r = r, c = c, deltas = deltas,
                         Psi.inv = Psi.inv, j)

      if (prior == "spn" | prior == "sng") {
        r[j] <- solve.kappa(kappa = kappa, prior = prior)
      } else if (prior == "spb") {

        ce <- kappa
        ce[6] <- ceiling(ce[6])
        low <- solve.kappa(kappa = ce, prior = prior)

        fl <- kappa
        fl[6] <- floor(fl[6])
        hig <- solve.kappa(kappa = fl, prior = prior)

        if (is.na(low) & is.na(hig)) {
          r[j] <- r[j]
        } else if (is.na(low)) {
          r[j] <- hig
        } else if (is.na(hig)) {
          r[j] <- low
        } else {

        g.low <- kappa.ll(s.j = low, kappa = kappa)
        g.hig <- kappa.ll(s.j = hig, kappa = kappa)

        inner <- 1
        mid <- (hig + low)/2
        while (abs(hig - low) > eps & inner <= max.inner) {
          mid <- (low + hig)/2
          g.mid <- kappa.ll(s.j = mid, kappa = kappa)
          if (g.mid < 0) {
            hig <- mid
            g.hig <- g.mid
          } else {
            low <- mid
            g.low <- g.mid
          }
          inner <- inner + 1
        }
        r[j] <- mid
        }
      }
    }

    if (prior == "spb") {
      objs[i] <- h.log.r.spb(theta = 0, d0 = rep(0, length(r)),
                             d1 = r, Omega.inv = Omega.inv, beta = beta, c = c,
                             deltas = deltas)/length(r)
    } else if (prior == "sng") {
      objs[i] <- h.log.r.sng(theta = 0, d0 = rep(0, length(r)),
                             d1 = r, Omega.inv = Omega.inv, beta = beta, c = c)/length(r)
    } else if (prior == "spn") {
      objs[i] <- h.log.r.spn(theta = 0, d0 = rep(0, length(r)),
                             d1 = r, Omega.inv = Omega.inv, beta = beta, Psi.inv = Psi.inv)/length(r)
    }
  if (print.iter) {
    cat("  i = ", i, "\n")
    cat("obj = ", objs[i], "\n")
  }
  if (i > 1) {
    if (sum(r == 0) == 0 & abs(objs[i] - objs[i - 1]) < eps) {
      break
    }
  }
  }
  return(list("r" = r, "objs" = objs[1:i]))
}
