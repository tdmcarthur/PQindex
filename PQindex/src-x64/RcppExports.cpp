// This file was generated by Rcpp::compileAttributes
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

// linspacetest
arma::mat linspacetest(arma::vec x);
RcppExport SEXP PQindex_linspacetest(SEXP xSEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< arma::vec >::type x(xSEXP);
    __result = Rcpp::wrap(linspacetest(x));
    return __result;
END_RCPP
}
// regspacetest
arma::mat regspacetest(arma::vec x);
RcppExport SEXP PQindex_regspacetest(SEXP xSEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< arma::vec >::type x(xSEXP);
    __result = Rcpp::wrap(regspacetest(x));
    return __result;
END_RCPP
}
// regspacetest2
arma::mat regspacetest2(arma::vec x);
RcppExport SEXP PQindex_regspacetest2(SEXP xSEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< arma::vec >::type x(xSEXP);
    __result = Rcpp::wrap(regspacetest2(x));
    return __result;
END_RCPP
}
// sugar_in
std::vector<int> sugar_in(IntegerVector x, IntegerVector y);
RcppExport SEXP PQindex_sugar_in(SEXP xSEXP, SEXP ySEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< IntegerVector >::type x(xSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type y(ySEXP);
    __result = Rcpp::wrap(sugar_in(x, y));
    return __result;
END_RCPP
}
// sugar_in2
std::vector<int> sugar_in2(IntegerVector x, IntegerVector y);
RcppExport SEXP PQindex_sugar_in2(SEXP xSEXP, SEXP ySEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< IntegerVector >::type x(xSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type y(ySEXP);
    __result = Rcpp::wrap(sugar_in2(x, y));
    return __result;
END_RCPP
}
// myInOperator
arma::uvec myInOperator(const arma::uvec myBigVec, const arma::uvec mySmallVec);
RcppExport SEXP PQindex_myInOperator(SEXP myBigVecSEXP, SEXP mySmallVecSEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< const arma::uvec >::type myBigVec(myBigVecSEXP);
    Rcpp::traits::input_parameter< const arma::uvec >::type mySmallVec(mySmallVecSEXP);
    __result = Rcpp::wrap(myInOperator(myBigVec, mySmallVec));
    return __result;
END_RCPP
}
// fisherIndfastestfurious
arma::mat fisherIndfastestfurious(const arma::sp_mat Q_consol, const arma::sp_mat P_consol, arma::mat Q_freq, arma::mat P_freq, arma::uvec Q_ind, arma::uvec P_ind);
RcppExport SEXP PQindex_fisherIndfastestfurious(SEXP Q_consolSEXP, SEXP P_consolSEXP, SEXP Q_freqSEXP, SEXP P_freqSEXP, SEXP Q_indSEXP, SEXP P_indSEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< const arma::sp_mat >::type Q_consol(Q_consolSEXP);
    Rcpp::traits::input_parameter< const arma::sp_mat >::type P_consol(P_consolSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type Q_freq(Q_freqSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type P_freq(P_freqSEXP);
    Rcpp::traits::input_parameter< arma::uvec >::type Q_ind(Q_indSEXP);
    Rcpp::traits::input_parameter< arma::uvec >::type P_ind(P_indSEXP);
    __result = Rcpp::wrap(fisherIndfastestfurious(Q_consol, P_consol, Q_freq, P_freq, Q_ind, P_ind));
    return __result;
END_RCPP
}
// fisherIndfastest
arma::mat fisherIndfastest(arma::sp_mat Q_consol, arma::sp_mat P_consol, arma::mat Q_freq, arma::mat P_freq, arma::uvec Q_ind, arma::uvec P_ind);
RcppExport SEXP PQindex_fisherIndfastest(SEXP Q_consolSEXP, SEXP P_consolSEXP, SEXP Q_freqSEXP, SEXP P_freqSEXP, SEXP Q_indSEXP, SEXP P_indSEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< arma::sp_mat >::type Q_consol(Q_consolSEXP);
    Rcpp::traits::input_parameter< arma::sp_mat >::type P_consol(P_consolSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type Q_freq(Q_freqSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type P_freq(P_freqSEXP);
    Rcpp::traits::input_parameter< arma::uvec >::type Q_ind(Q_indSEXP);
    Rcpp::traits::input_parameter< arma::uvec >::type P_ind(P_indSEXP);
    __result = Rcpp::wrap(fisherIndfastest(Q_consol, P_consol, Q_freq, P_freq, Q_ind, P_ind));
    return __result;
END_RCPP
}
// fisherInd
arma::mat fisherInd(arma::mat Q, arma::mat P, int base_period);
RcppExport SEXP PQindex_fisherInd(SEXP QSEXP, SEXP PSEXP, SEXP base_periodSEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< arma::mat >::type Q(QSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type P(PSEXP);
    Rcpp::traits::input_parameter< int >::type base_period(base_periodSEXP);
    __result = Rcpp::wrap(fisherInd(Q, P, base_period));
    return __result;
END_RCPP
}
// fisherIndfast
arma::mat fisherIndfast(arma::mat Q, arma::mat P, arma::mat Q_consol, arma::mat P_consol, arma::mat Q_freq, arma::mat P_freq);
RcppExport SEXP PQindex_fisherIndfast(SEXP QSEXP, SEXP PSEXP, SEXP Q_consolSEXP, SEXP P_consolSEXP, SEXP Q_freqSEXP, SEXP P_freqSEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< arma::mat >::type Q(QSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type P(PSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type Q_consol(Q_consolSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type P_consol(P_consolSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type Q_freq(Q_freqSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type P_freq(P_freqSEXP);
    __result = Rcpp::wrap(fisherIndfast(Q, P, Q_consol, P_consol, Q_freq, P_freq));
    return __result;
END_RCPP
}
// testfn
arma::mat testfn(arma::mat Q_consol, arma::uvec Q_ind);
RcppExport SEXP PQindex_testfn(SEXP Q_consolSEXP, SEXP Q_indSEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< arma::mat >::type Q_consol(Q_consolSEXP);
    Rcpp::traits::input_parameter< arma::uvec >::type Q_ind(Q_indSEXP);
    __result = Rcpp::wrap(testfn(Q_consol, Q_ind));
    return __result;
END_RCPP
}
// matmult_cpp
arma::sp_mat matmult_cpp(const arma::sp_mat X, const arma::sp_mat Y);
RcppExport SEXP PQindex_matmult_cpp(SEXP XSEXP, SEXP YSEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< const arma::sp_mat >::type X(XSEXP);
    Rcpp::traits::input_parameter< const arma::sp_mat >::type Y(YSEXP);
    __result = Rcpp::wrap(matmult_cpp(X, Y));
    return __result;
END_RCPP
}
// fisherIndfaster
arma::mat fisherIndfaster(arma::mat Q_consol, arma::mat P_consol, arma::mat Q_freq, arma::mat P_freq, arma::uvec Q_ind, arma::uvec P_ind);
RcppExport SEXP PQindex_fisherIndfaster(SEXP Q_consolSEXP, SEXP P_consolSEXP, SEXP Q_freqSEXP, SEXP P_freqSEXP, SEXP Q_indSEXP, SEXP P_indSEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< arma::mat >::type Q_consol(Q_consolSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type P_consol(P_consolSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type Q_freq(Q_freqSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type P_freq(P_freqSEXP);
    Rcpp::traits::input_parameter< arma::uvec >::type Q_ind(Q_indSEXP);
    Rcpp::traits::input_parameter< arma::uvec >::type P_ind(P_indSEXP);
    __result = Rcpp::wrap(fisherIndfaster(Q_consol, P_consol, Q_freq, P_freq, Q_ind, P_ind));
    return __result;
END_RCPP
}
// fisherIndfasterold
arma::mat fisherIndfasterold(arma::mat Q_consol, arma::mat P_consol, arma::mat Q_freq, arma::mat P_freq, arma::uvec Q_ind, arma::uvec P_ind);
RcppExport SEXP PQindex_fisherIndfasterold(SEXP Q_consolSEXP, SEXP P_consolSEXP, SEXP Q_freqSEXP, SEXP P_freqSEXP, SEXP Q_indSEXP, SEXP P_indSEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< arma::mat >::type Q_consol(Q_consolSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type P_consol(P_consolSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type Q_freq(Q_freqSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type P_freq(P_freqSEXP);
    Rcpp::traits::input_parameter< arma::uvec >::type Q_ind(Q_indSEXP);
    Rcpp::traits::input_parameter< arma::uvec >::type P_ind(P_indSEXP);
    __result = Rcpp::wrap(fisherIndfasterold(Q_consol, P_consol, Q_freq, P_freq, Q_ind, P_ind));
    return __result;
END_RCPP
}
// fisherIndfullmat
arma::mat fisherIndfullmat(arma::mat Q, arma::mat P);
RcppExport SEXP PQindex_fisherIndfullmat(SEXP QSEXP, SEXP PSEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< arma::mat >::type Q(QSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type P(PSEXP);
    __result = Rcpp::wrap(fisherIndfullmat(Q, P));
    return __result;
END_RCPP
}
// rcpp_hello
List rcpp_hello();
RcppExport SEXP PQindex_rcpp_hello() {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    __result = Rcpp::wrap(rcpp_hello());
    return __result;
END_RCPP
}
