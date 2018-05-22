obj.r.spn <- function(r, Omega.inv, Psi.inv, beta) {
  p <- unlist(lapply(Omega.inv, function(x) {nrow(x)}))
  c1 <- -c(crossprod(Omega.inv[[1]], tcrossprod(array(r*beta, dim = p), Omega.inv[[2]])))%*%(r*beta)/2
  c2 <- -c(crossprod(Psi.inv[[1]], tcrossprod(array(1/r, dim = p), Psi.inv[[2]])))%*%(1/r)/2
  c3 <- -sum(log(abs(r)))
  return(c1 + c2 + c3)
}

f.r.spn <- function(r.j, r.nj, beta, Omega.tilde.inv.j, Psi.inv.j, j) {
  # Omega.tilde.inv <- Omega.inv*tcrossprod(beta)
  (-1/2)*(r.j^2*Omega.tilde.inv.j[j] + r.j*sum(2*Omega.tilde.inv.j[-j]%*%r.nj) + Psi.inv.j[j]/r.j^2 + (1/r.j)*sum(2*Psi.inv.j[-j]%*%(1/r.nj))) - log(abs(r.j))
}

update.r.spn <- function(r, beta, Omega.tilde.inv.j, Psi.inv.j, j) {
  
  # Polynomial coefficients for solution
  # Omega.tilde.inv <- Omega.inv*tcrossprod(beta)
  poly.coef <- c(-2*Psi.inv.j[j], -sum(2*Psi.inv.j[-j]/r[-j]), 2, sum(2*Omega.tilde.inv.j[-j]*r[-j]), 2*Omega.tilde.inv.j[j])
  pol <- polynomial(poly.coef)
  sol <- solve(pol)
  sol <- Re(sol[which(Im(sol) == 0 & Re(sol) > 0)])
  if (length(sol) > 0) {
    val <- f.r.spn(sol, r.nj = r[-j], beta = beta, Omega.tilde.inv.j = Omega.tilde.inv.j, Psi.inv.j = Psi.inv.j, j = j)
    sol <- sol[which(val == max(val))]
  } else {
    # Use a new random value if no solution
    sol <- rnorm(1)
  }
  return(sol)
  
}

coord.desc.r.spn <- function(beta, Omega.inv, Psi.inv, eps = 10^(-12), max.iter = 1000, print.iter = TRUE,
                             start.r.spn = rep(1, length(beta))) {
  
  r <- start.r.spn
  objs <- rep(NA, max.iter)
  for (i in 1:max.iter) {
    for (j in 1:length(beta)) {
      if (print.iter) {
        # cat("  j = ", j, "\n")
      }
      p <- unlist(lapply(Omega.inv, function(x) {nrow(x)}))
      Omega.tilde.inv.j <- rep(Omega.inv[[2]][ifelse(j%%p[1] == 0, j/p[1], floor(j/p[1]) + 1), ], each = p[1])*rep(Omega.inv[[1]][ifelse(j%%p[1] != 0, j%%p[1], p[1]), ], times = p[2])*beta[j]*beta
      Psi.inv.j <- rep(Psi.inv[[2]][ifelse(j%%p[1] == 0, j/p[1], floor(j/p[1]) + 1), ], each = p[1])*rep(Psi.inv[[1]][ifelse(j%%p[1] != 0, j%%p[1], p[1]), ], times = p[2])
      
      r[j] <- update.r.spn(r = r, beta = beta, Omega.tilde.inv.j = Omega.tilde.inv.j, Psi.inv.j = Psi.inv.j, j = j)
    }
    objs[i] <- obj.r.spn(r = r, Omega.inv = Omega.inv, beta = beta, Psi.inv = Psi.inv)/length(r)
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
  return(list("r" = r, "objs" = objs[1:i]))
}