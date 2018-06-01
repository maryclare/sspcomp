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
