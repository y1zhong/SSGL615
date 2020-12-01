#include <Rcpp.h>
using namespace Rcpp;

// This is a simple example of exporting a C++ function to R. You can
// source this function into an R session using the Rcpp::sourceCpp
// function (or via the Source button on the editor toolbar). Learn
// more about Rcpp at:
//
//   http://www.rcpp.org/
//   http://adv-r.had.co.nz/Rcpp.html
//   http://gallery.rcpp.org/
//
// [[Rcpp::export]]
NumericVector c_which2(NumericVector groups, int g){
  // calling rnorm()
  Function f("which");

  // Next code is interpreted as rnorm(n=5, mean=10, sd=2)
  return f(groups == g);
}
// [[Rcpp::export]]
NumericVector c_which3(NumericVector beta, NumericVector active){
  // calling rnorm()
  NumericVector ind {Range}
  Function f("which");

  return f((beta != 0) & (!1:beta.length() %in% active));
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
NumericVector vector_subset_idx_rcpp2(
    NumericVector x, NumericVector y) {

  int l = y.length();
  NumericVector out;
  for(unsigned int z = 0; z < l; z++) {
    out.push_back(x[z]);
  }
  return out;
}

// [[Rcpp::export]]
NumericVector c_in(
    NumericVector beta, NumericVector active) {

  int l = beta.length();
  LogicalVector b0 = (beta != 0);
  NumericVector blen(l);
  blen = Range(1, l);
  for(unsigned int z = 0; z < l; z++) {
    out.push_back(x[z]);
  }
  return out;
}

// [[Rcpp::export]]
List update(NumericMatrix Y, NumericMatrix Xtilde, NumericVector groups,
            double sigmasq, NumericVector beta, double intercept, double lambda0_base,
            NumericVector Z, double theta, int G, NumericVector forceGroups){
  for(int g = 0; g < G; g++){
    NumericVector active2 = c_which2(beta, 0);
    if (active2.length() == 0) {
      intercept = mean(Y);
    } else if (active2.length() == 1) {
      int act = active2[0];
      intercept = mean(Y - Xtilde(_,act) * beta[act]);
    } else {
      NumericVector tempXX = vector_subset_idx_rcpp2(beta, active2-1);
      NumericMatrix tempXtilde = matrix_subset_idx_rcpp2(Xtilde, active2-1);
      intercept = mean(Y - tempXtilde * tempXX);
    }

    // which parameters refer to this group
    NumericVector active = c_which2(groups, g);
    int m = active.length();
    double lambda0 = sqrt(m) * lambda0_base;

    if ( std::find(forceGroups.begin(), forceGroups.end(), g) != forceGroups.end() ) {
      active2 = c_which3(beta, active);
    } else {

    }
  }
}
