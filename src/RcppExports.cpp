// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

// baseSVD
arma::vec baseSVD(const arma::mat& X);
RcppExport SEXP _SSGL615_baseSVD(SEXP XSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::mat& >::type X(XSEXP);
    rcpp_result_gen = Rcpp::wrap(baseSVD(X));
    return rcpp_result_gen;
END_RCPP
}
// c_which
NumericVector c_which(NumericVector groups, int g);
RcppExport SEXP _SSGL615_c_which(SEXP groupsSEXP, SEXP gSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type groups(groupsSEXP);
    Rcpp::traits::input_parameter< int >::type g(gSEXP);
    rcpp_result_gen = Rcpp::wrap(c_which(groups, g));
    return rcpp_result_gen;
END_RCPP
}
// c_crossprod
NumericMatrix c_crossprod(NumericMatrix X);
RcppExport SEXP _SSGL615_c_crossprod(SEXP XSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericMatrix >::type X(XSEXP);
    rcpp_result_gen = Rcpp::wrap(c_crossprod(X));
    return rcpp_result_gen;
END_RCPP
}
// matrix_subset_idx_rcpp
Rcpp::NumericMatrix matrix_subset_idx_rcpp(Rcpp::NumericMatrix x, Rcpp::NumericVector y);
RcppExport SEXP _SSGL615_matrix_subset_idx_rcpp(SEXP xSEXP, SEXP ySEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::NumericMatrix >::type x(xSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type y(ySEXP);
    rcpp_result_gen = Rcpp::wrap(matrix_subset_idx_rcpp(x, y));
    return rcpp_result_gen;
END_RCPP
}
// c_Xtilde
NumericMatrix c_Xtilde(NumericMatrix Xstar, NumericVector groups, int G, int n);
RcppExport SEXP _SSGL615_c_Xtilde(SEXP XstarSEXP, SEXP groupsSEXP, SEXP GSEXP, SEXP nSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericMatrix >::type Xstar(XstarSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type groups(groupsSEXP);
    Rcpp::traits::input_parameter< int >::type G(GSEXP);
    Rcpp::traits::input_parameter< int >::type n(nSEXP);
    rcpp_result_gen = Rcpp::wrap(c_Xtilde(Xstar, groups, G, n));
    return rcpp_result_gen;
END_RCPP
}
// c_Xstar
NumericMatrix c_Xstar(NumericMatrix Xs);
RcppExport SEXP _SSGL615_c_Xstar(SEXP XsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericMatrix >::type Xs(XsSEXP);
    rcpp_result_gen = Rcpp::wrap(c_Xstar(Xs));
    return rcpp_result_gen;
END_RCPP
}
// c_which2
NumericVector c_which2(LogicalVector v);
RcppExport SEXP _SSGL615_c_which2(SEXP vSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< LogicalVector >::type v(vSEXP);
    rcpp_result_gen = Rcpp::wrap(c_which2(v));
    return rcpp_result_gen;
END_RCPP
}
// matrix_subset_idx_rcpp2
NumericMatrix matrix_subset_idx_rcpp2(NumericMatrix x, NumericVector y);
RcppExport SEXP _SSGL615_matrix_subset_idx_rcpp2(SEXP xSEXP, SEXP ySEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericMatrix >::type x(xSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type y(ySEXP);
    rcpp_result_gen = Rcpp::wrap(matrix_subset_idx_rcpp2(x, y));
    return rcpp_result_gen;
END_RCPP
}
// matrix_assign_rcpp2
NumericMatrix matrix_assign_rcpp2(NumericMatrix X, NumericMatrix Y, NumericVector id);
RcppExport SEXP _SSGL615_matrix_assign_rcpp2(SEXP XSEXP, SEXP YSEXP, SEXP idSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericMatrix >::type X(XSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type Y(YSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type id(idSEXP);
    rcpp_result_gen = Rcpp::wrap(matrix_assign_rcpp2(X, Y, id));
    return rcpp_result_gen;
END_RCPP
}
// vector_assign_rcpp2
NumericVector vector_assign_rcpp2(NumericVector X, NumericVector Y, NumericVector id);
RcppExport SEXP _SSGL615_vector_assign_rcpp2(SEXP XSEXP, SEXP YSEXP, SEXP idSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type X(XSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type Y(YSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type id(idSEXP);
    rcpp_result_gen = Rcpp::wrap(vector_assign_rcpp2(X, Y, id));
    return rcpp_result_gen;
END_RCPP
}
// vector_subset_idx_rcpp2
NumericVector vector_subset_idx_rcpp2(NumericVector x, NumericVector y);
RcppExport SEXP _SSGL615_vector_subset_idx_rcpp2(SEXP xSEXP, SEXP ySEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type x(xSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type y(ySEXP);
    rcpp_result_gen = Rcpp::wrap(vector_subset_idx_rcpp2(x, y));
    return rcpp_result_gen;
END_RCPP
}
// c_betainact
LogicalVector c_betainact(NumericVector beta, NumericVector active);
RcppExport SEXP _SSGL615_c_betainact(SEXP betaSEXP, SEXP activeSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type beta(betaSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type active(activeSEXP);
    rcpp_result_gen = Rcpp::wrap(c_betainact(beta, active));
    return rcpp_result_gen;
END_RCPP
}
// psi
double psi(NumericVector beta, double lambda);
RcppExport SEXP _SSGL615_psi(SEXP betaSEXP, SEXP lambdaSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type beta(betaSEXP);
    Rcpp::traits::input_parameter< double >::type lambda(lambdaSEXP);
    rcpp_result_gen = Rcpp::wrap(psi(beta, lambda));
    return rcpp_result_gen;
END_RCPP
}
// pStar
double pStar(NumericVector beta, double lambda1, double lambda0, double theta);
RcppExport SEXP _SSGL615_pStar(SEXP betaSEXP, SEXP lambda1SEXP, SEXP lambda0SEXP, SEXP thetaSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type beta(betaSEXP);
    Rcpp::traits::input_parameter< double >::type lambda1(lambda1SEXP);
    Rcpp::traits::input_parameter< double >::type lambda0(lambda0SEXP);
    Rcpp::traits::input_parameter< double >::type theta(thetaSEXP);
    rcpp_result_gen = Rcpp::wrap(pStar(beta, lambda1, lambda0, theta));
    return rcpp_result_gen;
END_RCPP
}
// c_gFunc
double c_gFunc(NumericVector beta, double lambda1, double lambda0, double theta, double sigmasq, double n);
RcppExport SEXP _SSGL615_c_gFunc(SEXP betaSEXP, SEXP lambda1SEXP, SEXP lambda0SEXP, SEXP thetaSEXP, SEXP sigmasqSEXP, SEXP nSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type beta(betaSEXP);
    Rcpp::traits::input_parameter< double >::type lambda1(lambda1SEXP);
    Rcpp::traits::input_parameter< double >::type lambda0(lambda0SEXP);
    Rcpp::traits::input_parameter< double >::type theta(thetaSEXP);
    Rcpp::traits::input_parameter< double >::type sigmasq(sigmasqSEXP);
    Rcpp::traits::input_parameter< double >::type n(nSEXP);
    rcpp_result_gen = Rcpp::wrap(c_gFunc(beta, lambda1, lambda0, theta, sigmasq, n));
    return rcpp_result_gen;
END_RCPP
}
// update
List update(NumericVector Y, NumericMatrix Xtilde, NumericVector groups, LogicalVector updateSigma, double sigmasq, NumericVector beta, double intercept, double lambda0_base, double lambda1, double lambda0, NumericVector betaOld, double a, double b, int M, LogicalVector Z, double theta, int G, IntegerVector forceGroups, double n);
RcppExport SEXP _SSGL615_update(SEXP YSEXP, SEXP XtildeSEXP, SEXP groupsSEXP, SEXP updateSigmaSEXP, SEXP sigmasqSEXP, SEXP betaSEXP, SEXP interceptSEXP, SEXP lambda0_baseSEXP, SEXP lambda1SEXP, SEXP lambda0SEXP, SEXP betaOldSEXP, SEXP aSEXP, SEXP bSEXP, SEXP MSEXP, SEXP ZSEXP, SEXP thetaSEXP, SEXP GSEXP, SEXP forceGroupsSEXP, SEXP nSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type Y(YSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type Xtilde(XtildeSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type groups(groupsSEXP);
    Rcpp::traits::input_parameter< LogicalVector >::type updateSigma(updateSigmaSEXP);
    Rcpp::traits::input_parameter< double >::type sigmasq(sigmasqSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type beta(betaSEXP);
    Rcpp::traits::input_parameter< double >::type intercept(interceptSEXP);
    Rcpp::traits::input_parameter< double >::type lambda0_base(lambda0_baseSEXP);
    Rcpp::traits::input_parameter< double >::type lambda1(lambda1SEXP);
    Rcpp::traits::input_parameter< double >::type lambda0(lambda0SEXP);
    Rcpp::traits::input_parameter< NumericVector >::type betaOld(betaOldSEXP);
    Rcpp::traits::input_parameter< double >::type a(aSEXP);
    Rcpp::traits::input_parameter< double >::type b(bSEXP);
    Rcpp::traits::input_parameter< int >::type M(MSEXP);
    Rcpp::traits::input_parameter< LogicalVector >::type Z(ZSEXP);
    Rcpp::traits::input_parameter< double >::type theta(thetaSEXP);
    Rcpp::traits::input_parameter< int >::type G(GSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type forceGroups(forceGroupsSEXP);
    Rcpp::traits::input_parameter< double >::type n(nSEXP);
    rcpp_result_gen = Rcpp::wrap(update(Y, Xtilde, groups, updateSigma, sigmasq, beta, intercept, lambda0_base, lambda1, lambda0, betaOld, a, b, M, Z, theta, G, forceGroups, n));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_SSGL615_baseSVD", (DL_FUNC) &_SSGL615_baseSVD, 1},
    {"_SSGL615_c_which", (DL_FUNC) &_SSGL615_c_which, 2},
    {"_SSGL615_c_crossprod", (DL_FUNC) &_SSGL615_c_crossprod, 1},
    {"_SSGL615_matrix_subset_idx_rcpp", (DL_FUNC) &_SSGL615_matrix_subset_idx_rcpp, 2},
    {"_SSGL615_c_Xtilde", (DL_FUNC) &_SSGL615_c_Xtilde, 4},
    {"_SSGL615_c_Xstar", (DL_FUNC) &_SSGL615_c_Xstar, 1},
    {"_SSGL615_c_which2", (DL_FUNC) &_SSGL615_c_which2, 1},
    {"_SSGL615_matrix_subset_idx_rcpp2", (DL_FUNC) &_SSGL615_matrix_subset_idx_rcpp2, 2},
    {"_SSGL615_matrix_assign_rcpp2", (DL_FUNC) &_SSGL615_matrix_assign_rcpp2, 3},
    {"_SSGL615_vector_assign_rcpp2", (DL_FUNC) &_SSGL615_vector_assign_rcpp2, 3},
    {"_SSGL615_vector_subset_idx_rcpp2", (DL_FUNC) &_SSGL615_vector_subset_idx_rcpp2, 2},
    {"_SSGL615_c_betainact", (DL_FUNC) &_SSGL615_c_betainact, 2},
    {"_SSGL615_psi", (DL_FUNC) &_SSGL615_psi, 2},
    {"_SSGL615_pStar", (DL_FUNC) &_SSGL615_pStar, 4},
    {"_SSGL615_c_gFunc", (DL_FUNC) &_SSGL615_c_gFunc, 6},
    {"_SSGL615_update", (DL_FUNC) &_SSGL615_update, 19},
    {NULL, NULL, 0}
};

RcppExport void R_init_SSGL615(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}