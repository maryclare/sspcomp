// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; indent-tabs-mode: nil; -*-

// we only include RcppArmadillo.h which pulls Rcpp.h in for us
#include "RcppArmadillo.h"

// via the depends attribute we tell Rcpp to create hooks for
// RcppArmadillo so that the build process will know what to do
//
// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;

// [[Rcpp::export]]
double objVec(arma::colvec beta, arma::colvec Xbeta, arma::colvec y,
              arma::colvec OmegaInv) {

  double loglik = sum(Xbeta) - as_scalar(y.t()*Xbeta) + sum(log(1.0 + exp(-1.0*Xbeta)));

  loglik += as_scalar(beta.t()*(beta%OmegaInv)/2.0);

  return loglik;
}

// [[Rcpp::export]]
double objMat(arma::colvec beta, arma::colvec Xbeta, arma::colvec y,
              arma::mat OmegaInv) {


  double loglik = sum(Xbeta) - as_scalar(y.t()*Xbeta) + sum(log(1.0 + exp(-1.0*Xbeta)));

  loglik += as_scalar(beta.t()*OmegaInv*(beta)/2.0);

  return loglik;
}

// [[Rcpp::export]]
double gradVec(arma::mat X, arma::colvec y, arma::colvec OmegaInv,
               arma::colvec beta, arma::colvec Xbeta, int j) {

  return as_scalar(beta.row(j - 1)*OmegaInv.row(j - 1)) + as_scalar(X.col(j - 1).t()*((-1*y) + 1.0/(1.0 + exp(-Xbeta))));
}

// [[Rcpp::export]]
double gradMat(arma::mat X, arma::colvec y, arma::mat OmegaInv,
               arma::colvec beta, arma::colvec Xbeta, int j) {

  return as_scalar(beta.t()*OmegaInv.col(j - 1)) + as_scalar(X.col(j - 1).t()*((-1*y) + 1.0/(1.0 + exp(-Xbeta))));
}


// [[Rcpp::export]]
double hessVec(arma::mat X, arma::colvec OmegaInv, arma::colvec Xbeta, int j) {

  return as_scalar(OmegaInv.row(j - 1)) + sum(X.col(j - 1)%X.col(j - 1)%exp(Xbeta)/((1.0 + exp(Xbeta))%(1.0 + exp(Xbeta))));
}

// [[Rcpp::export]]
double hessMat(arma::mat X, arma::mat OmegaInv, arma::colvec Xbeta, int j) {

  return OmegaInv(j - 1, j - 1) + sum(X.col(j - 1)%X.col(j - 1)%exp(Xbeta)/((1.0 + exp(Xbeta))%(1.0 + exp(Xbeta))));
}

// [[Rcpp::export]]
double solveLogitVec(arma::colvec y, arma::colvec beta, arma::mat X,
                     arma::colvec OmegaInv, arma::colvec Xbeta, int j,
                     double epsInner) {

  arma::colvec beta_old = beta;
  arma::colvec b = beta_old;
  double diff = 100.0;
  double gr = 0.0;
  double he = 0.0;

  if (min(X.col(j - 1)) == 0.0 & max(X.col(j - 1)) == 0.0) {

    b.row(j - 1) = 0.0;
    Xbeta = Xbeta + X.col(j - 1)*(b.row(j - 1) - beta_old.row(j - 1));

  } else {

      while (diff > epsInner) {
        gr = gradVec(X = X, y = y, OmegaInv = OmegaInv, beta = b, Xbeta = Xbeta, j);
        he = hessVec(X = X, OmegaInv = OmegaInv, Xbeta = Xbeta, j);

        b.row(j - 1) = b.row(j - 1) - gr/he;
        diff = fabs(as_scalar(b.row(j - 1) - beta_old.row(j - 1)));
        Xbeta = Xbeta + X.col(j - 1)*(b.row(j - 1) - beta_old.row(j - 1));
        beta_old = b;
      }


  }
  return as_scalar(b.row(j - 1));
}

// [[Rcpp::export]]
double solveLogitMat(arma::colvec y, arma::colvec beta, arma::mat X,
                     arma::mat OmegaInv, arma::colvec Xbeta, int j,
                     double epsInner) {

  arma::colvec beta_old = beta;
  arma::colvec b = beta_old;
  double diff = 100.0;
  double gr = 0.0;
  double he = 0.0;

  if (min(X.col(j - 1)) == 0.0 & max(X.col(j - 1)) == 0.0) {

    b.row(j - 1) = 0.0;
    Xbeta = Xbeta + X.col(j - 1)*(b.row(j - 1) - beta_old.row(j - 1));

  } else {

    while (diff > epsInner) {
      gr = gradMat(X = X, y = y, OmegaInv = OmegaInv, beta = b, Xbeta = Xbeta, j);
      he = hessMat(X = X, OmegaInv = OmegaInv, Xbeta = Xbeta, j);

      b.row(j - 1) = b.row(j - 1) - gr/he;
      diff = fabs(as_scalar(b.row(j - 1) - beta_old.row(j - 1)));
      Xbeta = Xbeta + X.col(j - 1)*(b.row(j - 1) - beta_old.row(j - 1));
      beta_old = b;
    }


  }
  return as_scalar(b.row(j - 1));
}

// [[Rcpp::export]]
arma::colvec coordDescLogitVec(arma::colvec y, arma::mat X,
                            arma::colvec OmegaInv, int maxit, arma::colvec betaStart,
                            double eps, double epsInner) {
  arma::colvec beta = betaStart;
  int p = beta.size();
  arma::colvec Xbeta = X*beta;
  arma::colvec objs(maxit);
  arma::colvec beta_old = beta;
  double diff = 100.0;

  for (int i = 0; i < maxit; i++) {
    for (int j = 0; j < p; j++) {
      beta_old = beta;
      beta(j) = solveLogitVec(y = y, beta = beta, X = X, OmegaInv = OmegaInv, Xbeta = Xbeta, j + 1, epsInner);
      Xbeta = Xbeta + X.col(j)*(beta.row(j) - beta_old.row(j));
    }
    objs.row(i) = objVec(beta = beta, Xbeta = Xbeta, y = y, OmegaInv = OmegaInv);
    if (i > 0) {
      diff = fabs(as_scalar(objs.row(i) - objs.row(i - 1)));
    }
    if (diff < eps) {
      break;
    }
  }
  return beta;
}

// [[Rcpp::export]]
arma::colvec coordDescLogitMat(arma::colvec y, arma::mat X,
                               arma::mat OmegaInv, int maxit, arma::colvec betaStart,
                               double eps, double epsInner) {
  arma::colvec beta = betaStart;
  int p = beta.size();
  arma::colvec Xbeta = X*beta;
  arma::colvec objs(maxit);
  arma::colvec beta_old = beta;
  double diff = 100.0;

  for (int i = 0; i < maxit; i++) {
    for (int j = 0; j < p; j++) {
      beta_old = beta;
      beta(j) = solveLogitMat(y = y, beta = beta, X = X, OmegaInv = OmegaInv, Xbeta = Xbeta, j + 1, epsInner);
      Xbeta = Xbeta + X.col(j)*(beta.row(j) - beta_old.row(j));
    }
    objs.row(i) = objMat(beta = beta, Xbeta = Xbeta, y = y, OmegaInv = OmegaInv);
    if (i > 0) {
      diff = fabs(as_scalar(objs.row(i) - objs.row(i - 1)));
    }
    if (diff < eps) {
      break;
    }
  }
  return beta;
}

