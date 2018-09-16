// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppEigen.h>
#include <Rcpp.h>

using namespace Rcpp;

// pickCoordinates
NumericMatrix pickCoordinates(unsigned Dim, unsigned N, unsigned fe, ListOf<NumericMatrix> VT, NumericMatrix U);
RcppExport SEXP _gfimm_pickCoordinates(SEXP DimSEXP, SEXP NSEXP, SEXP feSEXP, SEXP VTSEXP, SEXP USEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< unsigned >::type Dim(DimSEXP);
    Rcpp::traits::input_parameter< unsigned >::type N(NSEXP);
    Rcpp::traits::input_parameter< unsigned >::type fe(feSEXP);
    Rcpp::traits::input_parameter< ListOf<NumericMatrix> >::type VT(VTSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type U(USEXP);
    rcpp_result_gen = Rcpp::wrap(pickCoordinates(Dim, N, fe, VT, U));
    return rcpp_result_gen;
END_RCPP
}
// tsolveAndMultiply
Eigen::MatrixXd tsolveAndMultiply(const Eigen::MatrixXd& A, const Eigen::MatrixXd& C);
RcppExport SEXP _gfimm_tsolveAndMultiply(SEXP ASEXP, SEXP CSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const Eigen::MatrixXd& >::type A(ASEXP);
    Rcpp::traits::input_parameter< const Eigen::MatrixXd& >::type C(CSEXP);
    rcpp_result_gen = Rcpp::wrap(tsolveAndMultiply(A, C));
    return rcpp_result_gen;
END_RCPP
}
// nullSpace
Rcpp::List nullSpace(const Eigen::MatrixXd M);
RcppExport SEXP _gfimm_nullSpace(SEXP MSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const Eigen::MatrixXd >::type M(MSEXP);
    rcpp_result_gen = Rcpp::wrap(nullSpace(M));
    return rcpp_result_gen;
END_RCPP
}
// QRdecomp
Rcpp::List QRdecomp(const Eigen::MatrixXd& M);
RcppExport SEXP _gfimm_QRdecomp(SEXP MSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const Eigen::MatrixXd& >::type M(MSEXP);
    rcpp_result_gen = Rcpp::wrap(QRdecomp(M));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_gfimm_pickCoordinates", (DL_FUNC) &_gfimm_pickCoordinates, 5},
    {"_gfimm_tsolveAndMultiply", (DL_FUNC) &_gfimm_tsolveAndMultiply, 2},
    {"_gfimm_nullSpace", (DL_FUNC) &_gfimm_nullSpace, 1},
    {"_gfimm_QRdecomp", (DL_FUNC) &_gfimm_QRdecomp, 1},
    {NULL, NULL, 0}
};

RcppExport void R_init_gfimm(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
