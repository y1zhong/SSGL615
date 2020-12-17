// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
using namespace Rcpp;
using namespace arma;

//' Center X matrix by column
//'
//' @param Xs input numeric matrix
// [[Rcpp::export]]
arma::mat c_Xstar(arma::mat Xs){
  int nr = Xs.n_rows;
  int nc = Xs.n_cols;
  arma::mat Xnew( nr , nc );
  for(int j=0; j < nc; j++){
    Xnew.col(j) = (Xs.col(j) - mean(Xs.col(j)));
  }
  return(Xnew);
}
