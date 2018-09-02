# Can just use h.log.r.spb and h.log.r.sng functions for objective
sp.ll <- function(rj, alpha1, alpha2, kappa1, kappa2, kappa3) {
  alpha1*rj^(-2) -log(abs(rj)) + kappa2*(rj^2) + kappa3*(rj) + alpha2/rj
}
# sp.ll.d <- Deriv(sp.ll, "rj")
sp.ll.d <- function (rj, alpha1, alpha2, kappa1, kappa2, kappa3) {
  (kappa1 - 1) * sign(rj)/abs(rj) + 2 * (alpha1 * rj) + alpha2 -
    (2 * (kappa2/rj) + kappa3)/rj^2
}
# sp.ll.dd <- Deriv(sp.ll.d, "rj")
sp.ll.dd <- function (rj, alpha1, alpha2, kappa1, kappa2, kappa3)
{
  .e1 <- 2 * (kappa2/rj)
  ((2 * (.e1 + kappa3) + .e1)/rj - (kappa1 - 1) * sign(rj)^2)/rj^2 +
    2 * alpha1
}


sn.ll <- function(rj, alpha1, alpha2, kappa1, kappa2, kappa3) {
  alpha1*rj^(-2) - (kappa1 - 1)*log(abs(rj)) - 2*log(abs(rj)) + kappa2*abs(rj)^(-2*kappa3) + alpha2/abs(rj)
}
# sn.ll.d <- Deriv(sn.ll, "rj")
sn.ll.d <- function (rj, alpha1, alpha2, kappa1, kappa2, kappa3)
{
  .e1 <- abs(rj)
  -(((1 + 2 * (kappa2 * kappa3/.e1^(2 * kappa3)) + kappa1)/.e1 +
       alpha2/rj^2) * sign(rj) + 2 * (alpha1/rj^3))
}
# sn.ll.dd <- Deriv(sn.ll.d, "rj")
sn.ll.dd <- function (rj, alpha1, alpha2, kappa1, kappa2, kappa3)
{
  .e3 <- abs(rj)^(2 * kappa3)
  .e4 <- rj^2
  .e5 <- sign(rj)
  (((1 + kappa1 + kappa2 * kappa3 * (2/.e3 + 4 * (kappa3/.e3))) *
      .e5 + 2 * (alpha2/rj)) * .e5 + 6 * (alpha1/.e4))/.e4
}

#' @export
coord.desc.r <- function(Omega.inv, beta, c = NULL, eps = 10^(-12), max.iter = 1000, print.iter = TRUE,
                             start.r = rep(1, length(beta)), max.inner = Inf, deltas = NULL, prior,
                         Psi.inv = NULL) {

  r <- start.r
  p <- length(r)
  objs <- rep(NA, max.iter)
  for (i in 1:max.iter) {
    for (j in 1:length(r)) {

      r.old <- r

      # Newton Raphson
      diff <- Inf
      inner <- 1
      reset <- 1
      while(abs(diff) > eps & inner <= max.inner) {

        Omega.inv.j <- get.kron.row(j, Omega = Omega.inv)
        alpha1 <- -beta[j]^2*Omega.inv.j[j]/2
        alpha2 <- -sum(beta[j]*Omega.inv.j[-j]*abs(r[-j])*beta[-j])

        if (prior %in% c("spb", "sng")) {
          kappa1 <- ifelse(prior == "spb", 2*(c/2)/(c/2 - 1), 1 - 2*c)
          kappa2 <- ifelse(prior == "spb", -(((2*gamma(3/c))/gamma(1/c))^(c/2/(1 - c/2))*f.deltas(deltas = deltas[j], c = c)) , -c)
          kappa3 <- ifelse(prior == "spb", c/2/(c/2 - 1), -1)

          hess <- sn.ll.dd(rj = r[j], alpha1 = alpha1, alpha2 = alpha2, kappa1 = kappa1, kappa2 = kappa2, kappa3 = kappa3)
          grad <- sn.ll.d(rj = r[j], alpha1 = alpha1, alpha2 = alpha2, kappa1 = kappa1, kappa2 = kappa2, kappa3 = kappa3)
        } else if (prior == "spn") {

          Psi.inv.j <- get.kron.row(j, Omega = Psi.inv)

          kappa1 <- 1/2
          kappa2 <- -Psi.inv.j[j]/2
          kappa3 <- -sum(Psi.inv.j[-j]/r[-j])

          hess <- sp.ll.dd(rj = r[j], alpha1 = alpha1, alpha2 = alpha2, kappa1 = kappa1, kappa2 = kappa2, kappa3 = kappa3)
          grad <- sp.ll.d(rj = r[j], alpha1 = alpha1, alpha2 = alpha2, kappa1 = kappa1, kappa2 = kappa2, kappa3 = kappa3)
        }
        r[j] <- r[j] - grad/hess
        diff <- mean(abs(r[j] - r.old[j]))

        inner <- inner + 1

        if (is.na(diff) | abs(r[j]) > 100) {
          r[j] <- abs(runif(1, 0, 1))
          diff <- Inf
          inner <- 1
          reset <- reset + 1
          if (reset > 1) {
            r[j] <- 0
            diff <- 0
          }
        }

        r.old <- r

      }
    }

    if (prior == "spb") {
      objs[i] <- h.log.r.spb(theta = 0, d0 = rep(0, length(r)),
                             d1 = r, Omega.inv = Omega.inv, beta = beta, c = c, r.tilde = rep(0, length(r)),
                             V.r.inv = rep(0, length(r)), deltas = deltas)/length(r)
    } else if (prior == "sng") {
      objs[i] <- h.log.r.sng(theta = 0, d0 = rep(0, length(r)),
                             d1 = r, Omega.inv = Omega.inv, beta = beta, c = c, r.tilde = rep(0, length(r)),
                             V.r.inv = rep(0, length(r)))/length(r)
    } else if (prior == "spn") {
      objs[i] <- h.log.r.spn(theta = 0, d0 = rep(0, length(r)),
                             d1 = r, Omega.inv = Omega.inv, beta = beta, Psi.inv = Psi.inv,
                             r.tilde = rep(0, length(r)),
                             V.r.inv = rep(0, length(r)))/length(r)
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
