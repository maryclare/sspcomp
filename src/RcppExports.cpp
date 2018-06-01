// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

// objVec
double objVec(arma::colvec beta, arma::colvec Xbeta, arma::colvec y, arma::colvec OmegaInv);
RcppExport SEXP _sspcomp_objVec(SEXP betaSEXP, SEXP XbetaSEXP, SEXP ySEXP, SEXP OmegaInvSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::colvec >::type beta(betaSEXP);
    Rcpp::traits::input_parameter< arma::colvec >::type Xbeta(XbetaSEXP);
    Rcpp::traits::input_parameter< arma::colvec >::type y(ySEXP);
    Rcpp::traits::input_parameter< arma::colvec >::type OmegaInv(OmegaInvSEXP);
    rcpp_result_gen = Rcpp::wrap(objVec(beta, Xbeta, y, OmegaInv));
    return rcpp_result_gen;
END_RCPP
}
// objMat
double objMat(arma::colvec beta, arma::colvec Xbeta, arma::colvec y, arma::mat OmegaInv);
RcppExport SEXP _sspcomp_objMat(SEXP betaSEXP, SEXP XbetaSEXP, SEXP ySEXP, SEXP OmegaInvSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::colvec >::type beta(betaSEXP);
    Rcpp::traits::input_parameter< arma::colvec >::type Xbeta(XbetaSEXP);
    Rcpp::traits::input_parameter< arma::colvec >::type y(ySEXP);
    Rcpp::traits::input_parameter< arma::mat >::type OmegaInv(OmegaInvSEXP);
    rcpp_result_gen = Rcpp::wrap(objMat(beta, Xbeta, y, OmegaInv));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_sspcomp_objVec", (DL_FUNC) &_sspcomp_objVec, 4},
    {"_sspcomp_objMat", (DL_FUNC) &_sspcomp_objMat, 4},
    {NULL, NULL, 0}
};

RcppExport void R_init_sspcomp(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
