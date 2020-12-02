#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
NumericVector c_which2(LogicalVector v){
  // calling rnorm()
  Function f("which");

  // Next code is interpreted as rnorm(n=5, mean=10, sd=2)
  return f(v);
}



// [[Rcpp::export]]
NumericMatrix matrix_subset_idx_rcpp2(
    NumericMatrix x, NumericVector y) {

  int n_cols_out = y.size();
  NumericMatrix out = no_init(x.nrow(), n_cols_out);
  for(unsigned int z = 0; z < n_cols_out; ++z) {
    out(_, z) = x(_, y[z]);
  }
  return out;
}


// [[Rcpp::export]]
NumericMatrix matrix_assign_rcpp2(
    NumericMatrix X, NumericMatrix Y, NumericVector id) {

  int l = id.length();
  for(unsigned int z = 0; z < l; z++) {
    int pos = id[z];
    NumericMatrix::Column col = X( _ , pos);
    col = Y(_,z);
  }
  return X;
}



// [[Rcpp::export]]
NumericVector vector_assign_rcpp2(
    NumericVector X, NumericVector Y, NumericVector id) {

  int l = id.length();

  for(unsigned int z = 0; z < l; z++) {
    int pos = id[z];
    X[pos] = Y[z];
  }
  return X;
}

// [[Rcpp::export]]
NumericVector vector_subset_idx_rcpp2(
    NumericVector x, NumericVector y) {

  int l = y.length();
  NumericVector out;
  for(unsigned int z = 0; z < l; z++) {
    int id = y[z];
    out.push_back(x[id]);
  }
  return out;
}
