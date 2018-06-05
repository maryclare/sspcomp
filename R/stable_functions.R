# library(Deriv)

# I think this is correct and more stable, it displays the expected behavior
# as alpha approaches 1

#' @export
g <- function(alpha, u) {
  (sin(alpha*u)/sin(u))^(1/(1 - alpha))*(sin((1 - alpha)*u)/sin(alpha*u))
}
gg <- function(alpha, u) {
  ((sin(alpha*u)/sin(u))^(1/(1 - alpha))*(sin((1 - alpha)*u)/sin(alpha*u)))^((alpha - 1)/(2*alpha))
}
# gg.d <- Deriv(gg, "u")
gg.d <- function (alpha, u) {
  .e1 <- 1 - alpha
  .e2 <- alpha * u
  .e3 <- sin(.e2)
  .e4 <- sin(u)
  .e5 <- u * .e1
  .e6 <- 1/.e1
  .e7 <- .e3/.e4
  .e8 <- sin(.e5)
  .e9 <- .e7^.e6
  .e10 <- alpha - 1
  .e11 <- alpha * cos(.e2)
  ((.e1 * cos(.e5) - .e11 * .e8/.e3) * .e9 + (.e11 - cos(u) *
                                                .e3/.e4) * .e7^(.e6 - 1) * .e8/(.e1 * .e4)) * (.e9 *
                                                                                                 .e8/.e3)^(.e10/(2 * alpha) - 1) * .e10/(2 * (alpha *
                                                                                                                                                .e3))
}
# gg.dd <- Deriv(gg.d, "u")
gg.dd <- function (alpha, u) {
  .e1 <- alpha * u
  .e2 <- 1 - alpha
  .e3 <- sin(.e1)
  .e4 <- sin(u)
  .e5 <- u * .e2
  .e6 <- 1/.e2
  .e7 <- .e3/.e4
  .e8 <- cos(.e1)
  .e9 <- sin(.e5)
  .e10 <- alpha * .e8
  .e11 <- cos(u)
  .e12 <- .e7^.e6
  .e13 <- .e6 - 1
  .e14 <- .e7^.e13
  .e15 <- .e2 * cos(.e5)
  .e16 <- .e10 - .e11 * .e3/.e4
  .e17 <- .e2 * .e4
  .e18 <- alpha - 1
  .e19 <- .e15 - .e10 * .e9/.e3
  .e20 <- .e18/(2 * alpha)
  .e21 <- .e20 - 1
  .e22 <- .e16 * .e14
  .e24 <- .e12 * .e9/.e3
  .e25 <- alpha * .e3
  .e27 <- .e19 * .e12 + .e22 * .e9/.e17
  .e28 <- .e24^.e21
  .e29 <- 2 * .e25
  .e30 <- alpha^2
  (((((.e19 * .e14 + .e15 * .e14) * .e16 + (.e13 * .e16^2 *
                                              .e7^(.e6 - 2)/.e4 - ((.e10 * .e11 - (.e11^2/.e4 + .e4) *
                                                                      .e3)/.e4 + .e30 * .e3) * .e14) * .e9)/.e17 - ((.e2^2 *
                                                                                                                       .e9 + alpha * (.e19 * .e8 - .e25 * .e9)/.e3) * .e12 +
                                                                                                                      .e2 * .e16 * .e11 * .e14 * .e9/.e17^2)) * .e28 + .e27 *
      ((.e22/.e17 - .e10 * .e12/.e3) * .e9 + .e15 * .e12) *
      .e21 * .e24^(.e20 - 2)/.e3)/.e29 - 2 * (.e30 * .e27 *
                                                .e28 * .e8/.e29^2)) * .e18
}

B <- function(x, alpha) {
  v <- numeric(length(x))
  for (i in 1:length(x)) {
    if (x[i] != 0) {
      v[i] <- sin(x[i])/(sin(alpha*x[i])^alpha*sin((1 - alpha)*x[i])^(1 - alpha))
    } else {
      v[i] <- (1 - alpha)^(-1 + alpha)*alpha^(-alpha)
    }
  }
  return(v)
}

dev.acc <- function(phi, alpha) {
  sig <- 1/sqrt((1/(2*alpha))*alpha*(1 - alpha))
  if (sig >= sqrt(2*pi)) {
    runif(1)*B(x = 0, alpha = alpha)^(1/(2*alpha)) <= B(x = phi, alpha = alpha)^(1/(2*alpha))
  } else {
    phi <= pi & runif(1)*B(x = 0, alpha = alpha)^(1/(2*alpha))*exp(-(phi/sig)^2/(2)) <= B(x = phi, alpha = alpha)^(1/(2*alpha))
  }
}

dev.phi <- function(n, alpha) {
  sig <- 1/sqrt((1/(2*alpha))*alpha*(1 - alpha))
  phi <- numeric(n)
  for (i in 1:n) {

    acc <- FALSE
    while (!acc) {
      if (sig >= sqrt(2*pi)) {
        phi[i] <- runif(1, 0, pi)
      } else {
        phi[i] <- sig*abs(rnorm(1))
      }
      acc <- dev.acc(phi = phi[i], alpha = alpha)
    }
  }

  return(phi)
}

dev.check.acc <- function(alpha, n = 5000) {
  sig <- 1/sqrt((1/(2*alpha))*alpha*(1 - alpha))
  phi <- numeric(n)
  accs <- numeric(n)
  for (i in 1:n) {

    if (sig >= sqrt(2*pi)) {
      phi[i] <- runif(1, 0, pi)
    } else {
      phi[i] <- sig*abs(rnorm(1))
    }
    accs[i] <- dev.acc(phi = phi[i], alpha = alpha)
  }

  return(mean(accs))
}

naive.phi <- function(n, alpha) {
  M <- (1 - alpha)^((alpha - 1)/(2*alpha))/sqrt(alpha)/(1/pi)
  phi <- numeric(n)
  for (i in 1:n) {
    phi[i] <- runif(1, 0, pi)
    while(runif(1) > gg(alpha = alpha, u = phi[i])/(M/pi)) {
      phi[i] <- runif(1, 0, pi)
    }
  }
  return(phi)
}
smart.acc <- function(phi, alpha, b = NULL, c = NULL) {
  if (is.null(b)) {
    b <- gg.d(u = pi/2, alpha = alpha)
  }
  if (is.null(c)) {
    c <- gg(u = pi/2, alpha = alpha) - b*pi/2
    if (c + b*pi < 0) {
      c <- -b*pi
    }
  }
  return(runif(1) < gg(u = phi, alpha = alpha)/(c + b*phi))
}

#' @export
smart.phi <- function(n, phi, alpha) {
  b <- gg.d(u = pi/2, alpha = alpha)
  c <- gg(u = pi/2, alpha = alpha) - b*pi/2
  if (c + b*pi < 0) {
    c <- -b*pi
  }
  phi <- numeric(n)
  # Generate proposal
  for (i in 1:n) {
    acc <- FALSE
    while (!acc) {
      u <- runif(1)
      d <- (c*pi + b*pi^2/2)
      aa <- (b/2)/d
      bb <- c/d
      cc <- -u
      p <- (-bb + c(-1, 1)*sqrt(bb^2 - 4*aa*cc))/(2*aa)
      phi[i] <- p[0 < p & p < pi]
      acc <- smart.acc(phi = phi[i], alpha = alpha, c = c, b = b)
    }
  }
  return(phi)
}
check.acc <- function(alpha, n = 10000) {
  b <- gg.d(u = pi/2, alpha = alpha)
  c <- gg(u = pi/2, alpha = alpha) - b*pi/2
  if (c + b*pi < 0) {
    c <- -b*pi
  }
  phi <- accs <- numeric(n)
  # Generate proposal
  for (i in 1:n) {
    u <- runif(1)
    d <- (c*pi + b*pi^2/2)
    aa <- (b/2)/d
    bb <- c/d
    cc <- -u
    p <- (-bb + c(-1, 1)*sqrt(bb^2 - 4*aa*cc))/(2*aa)
    phi[i] <- p[0 < p & p < pi]
    accs[i] <- smart.acc(phi = phi[i], alpha = alpha, c = c, b = b)
  }
  return(mean(accs))
}

int <- function(x, alpha) {
  d <- numeric(length(x))
  for (i in 1:length(x)) {
    d[i] <- (alpha/((1 - alpha)*pi))*x[i]^(1/(alpha - 1))*integrate(function(z) {exp(-x[i]^(alpha/(alpha - 1))*g(alpha = alpha, u = z))*g(alpha = alpha, u = z)},
                                                                    0, pi)$value
  }
  return(d)
}

p.tau.sq <- function(x, alpha) {
  d <- numeric(length(x))
  for (i in 1:length(x)) {
    d[i] <- (alpha/((1 - alpha)*pi))*x[i]^((3 - alpha)/(2*(alpha - 1)))*integrate(function(z) {exp(-x[i]^(alpha/(alpha - 1))*g(alpha = alpha, u = z))*g(alpha = alpha, u = z)},
                                                                                  0, pi)$value
  }
  return(d)
}

cfun <- function(x, alpha, t) {
  d <- numeric(length(x))
  for (i in 1:length(x)) {
    d[i] <- exp(-t*x)*(alpha/((1 - alpha)*pi))*x[i]^(1/(alpha - 1))*integrate(function(z) {exp(-x[i]^(alpha/(alpha - 1))*g(alpha = alpha, u = z))*g(alpha = alpha, u = z)},
                                                                              0, pi)$value
  }
  return(d)
}

g.s <- function(s, alpha) {
  d <- numeric(length(s))
  for (i in 1:length(s)) {
    d[i] <- (alpha/((1 - alpha)*pi))*s[i]^((2*alpha - 1)/(1 - alpha))*integrate(function(z) {exp(-s[i]^(alpha/(1 - alpha))*g(alpha = alpha, u = z))*g(alpha = alpha, u = z)},
                                                                                0, pi)$value
  }
  return(d)
}

p.s <- function(s, alpha) {
  d <- numeric(length(s))
  for (i in 1:length(s)) {
    d[i] <- s[i]^((3*alpha - 1)/(2*(1 - alpha)))*integrate(function(z) {exp(-s[i]^(alpha/(1 - alpha))*g(alpha = alpha, u = z))*g(alpha = alpha, u = z)},
                                                           0, pi)$value
  }
  return(d)
}

pi.s <- function(s) {
  exp(-s/4)/(2*sqrt(pi))
}

goal.dens <- function(x, alpha, q, beta) {
  t <- beta^2*gamma(3/q)/gamma(1/q)
  d <- numeric(length(x))
  for (i in 1:length(x)) {
    d[i] <- sqrt(x)*exp(-t*x)*(alpha/((1 - alpha)*pi))*x[i]^(1/(alpha - 1))*integrate(function(z) {exp(-x[i]^(alpha/(alpha - 1))*g(alpha = alpha, u = z))*g(alpha = alpha, u = z)},
                                                                                      0, pi)$value/sqrt(x)
  }
  return(d)
}
