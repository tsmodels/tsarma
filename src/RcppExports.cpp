// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppEigen.h>
#include <Rcpp.h>

using namespace Rcpp;

#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// armafilter
NumericVector armafilter(NumericVector y, NumericVector epsilon, NumericVector x, NumericVector initstate, NumericVector mu, NumericVector phi, NumericVector theta, IntegerVector model);
RcppExport SEXP _tsarma_armafilter(SEXP ySEXP, SEXP epsilonSEXP, SEXP xSEXP, SEXP initstateSEXP, SEXP muSEXP, SEXP phiSEXP, SEXP thetaSEXP, SEXP modelSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type y(ySEXP);
    Rcpp::traits::input_parameter< NumericVector >::type epsilon(epsilonSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type x(xSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type initstate(initstateSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type mu(muSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type phi(phiSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type theta(thetaSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type model(modelSEXP);
    rcpp_result_gen = Rcpp::wrap(armafilter(y, epsilon, x, initstate, mu, phi, theta, model));
    return rcpp_result_gen;
END_RCPP
}
// armafilter2
NumericVector armafilter2(NumericVector y, NumericVector epsilon, NumericVector x, NumericVector mu, NumericVector phi, NumericVector theta, IntegerVector model);
RcppExport SEXP _tsarma_armafilter2(SEXP ySEXP, SEXP epsilonSEXP, SEXP xSEXP, SEXP muSEXP, SEXP phiSEXP, SEXP thetaSEXP, SEXP modelSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type y(ySEXP);
    Rcpp::traits::input_parameter< NumericVector >::type epsilon(epsilonSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type x(xSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type mu(muSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type phi(phiSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type theta(thetaSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type model(modelSEXP);
    rcpp_result_gen = Rcpp::wrap(armafilter2(y, epsilon, x, mu, phi, theta, model));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_tsarma_armafilter", (DL_FUNC) &_tsarma_armafilter, 8},
    {"_tsarma_armafilter2", (DL_FUNC) &_tsarma_armafilter2, 7},
    {NULL, NULL, 0}
};

RcppExport void R_init_tsarma(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
