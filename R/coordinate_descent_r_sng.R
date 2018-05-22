obj.r.sng <- function(r, Omega.inv, beta, c) {
  p <- unlist(lapply(Omega.inv, function(x) {nrow(x)}))
  c1 <- -sum(c(crossprod(Omega.inv[[1]], tcrossprod(array(abs(r)*beta, dim = p), Omega.inv[[2]])))*(abs(r)*beta)/2)
  c2 <- -c*sum(abs(r))
  c3 <- -sum(c/r^2)
  return(c1 + c2 + c3)
}

f.r.sng <- function(r.j, r.nj, beta, Omega.inv.j, c, j) {
  (-1/2)*(r.j^2*Omega.inv.j[j]*beta[j]^2 + beta[j]*abs(r.j)*(sum(2*Omega.inv.j[-j]*beta[-j]*r.nj))) - c*log(abs(r.j)) - c/r.j^2
}

update.r.sng <- function(r, beta, Omega.inv.j, c, j) {
  # Polynomial coefficients for solution
  poly.coef <- c(2*c, 0, -c, -beta[j]*(2*sum(Omega.inv.j[-j]*(beta[-j]*r[-j])))/2,
                 -beta[j]^2*Omega.inv.j[j])
  pol <- polynomial(poly.coef)
  sol <- solve(pol)
  sol <- Re(sol[which(Im(sol) == 0 & Re(sol) >= 0)])
  val <- f.r.sng(sol, r.nj = r[-j], beta = beta, Omega.inv.j = Omega.inv.j, c = c, j = j)
  return(sol[which(val == max(val))])
  # return(sol[sample(1:length(sol), 1)])
  
}

coord.desc.r.sng <- function(beta, Omega.inv, c, eps = 10^(-12), max.iter = 1000, print.iter = TRUE,
                         start.r.sng = rep(1, length(beta))) {
  
  r <- start.r.sng
  p <- unlist(lapply(Omega.inv, function(x) {nrow(x)}))
  
  objs <- rep(NA, max.iter)
  for (i in 1:max.iter) {
    for (j in 1:length(beta)) {
      if (print.iter) {
        # cat("  j = ", j, "\n")
      }
      
      Omega.inv.j <- rep(Omega.inv[[2]][ifelse(j%%p[1] == 0, j/p[1], floor(j/p[1]) + 1), ], each = p[1])*rep(Omega.inv[[1]][ifelse(j%%p[1] != 0, j%%p[1], p[1]), ], times = p[2])
      
      r[j] <- update.r.sng(r = r, beta = beta, Omega.inv.j = Omega.inv.j, c = c, j = j)
      }
    objs[i] <- obj.r.sng(r = r, Omega.inv = Omega.inv, beta = beta, c = c)
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
