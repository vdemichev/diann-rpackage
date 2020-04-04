// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppEigen.h>
#include <Rcpp.h>

using namespace Rcpp;

// rcppeigen_hello_world
Eigen::MatrixXd rcppeigen_hello_world();
RcppExport SEXP _diann_rcppeigen_hello_world() {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    rcpp_result_gen = Rcpp::wrap(rcppeigen_hello_world());
    return rcpp_result_gen;
END_RCPP
}
// col_max
std::vector<double> col_max(std::vector<double>& quantities, int m, int n);
RcppExport SEXP _diann_col_max(SEXP quantitiesSEXP, SEXP mSEXP, SEXP nSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< std::vector<double>& >::type quantities(quantitiesSEXP);
    Rcpp::traits::input_parameter< int >::type m(mSEXP);
    Rcpp::traits::input_parameter< int >::type n(nSEXP);
    rcpp_result_gen = Rcpp::wrap(col_max(quantities, m, n));
    return rcpp_result_gen;
END_RCPP
}
// maxlfq_solve
std::vector<double> maxlfq_solve(std::vector<double>& quantities, int peptides, int samples, double margin);
RcppExport SEXP _diann_maxlfq_solve(SEXP quantitiesSEXP, SEXP peptidesSEXP, SEXP samplesSEXP, SEXP marginSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< std::vector<double>& >::type quantities(quantitiesSEXP);
    Rcpp::traits::input_parameter< int >::type peptides(peptidesSEXP);
    Rcpp::traits::input_parameter< int >::type samples(samplesSEXP);
    Rcpp::traits::input_parameter< double >::type margin(marginSEXP);
    rcpp_result_gen = Rcpp::wrap(maxlfq_solve(quantities, peptides, samples, margin));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_diann_rcppeigen_hello_world", (DL_FUNC) &_diann_rcppeigen_hello_world, 0},
    {"_diann_col_max", (DL_FUNC) &_diann_col_max, 3},
    {"_diann_maxlfq_solve", (DL_FUNC) &_diann_maxlfq_solve, 4},
    {NULL, NULL, 0}
};

RcppExport void R_init_diann(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}