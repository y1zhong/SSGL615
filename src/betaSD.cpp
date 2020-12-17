#include <RcppArmadillo.h>
using namespace Rcpp;
using namespace arma;

// [[Rcpp::depends(RcppArmadillo)]]

//' re-standardize the centered X matrix
//' @param Xstar centered numeric matrix
//' @param group vector group index
//' @param beta estimated parameters
//' @param Qmat,Dvec SVD results
//' @param G number of groups
//' @param n length of subjects
//' @return A standardized matrix
// [[Rcpp::export]]
arma::mat c_betaSD(arma::mat Xstar, arma::vec groups, arma::vec beta,
             List Qmat, List Dvec, int G, int n){
  arma::uvec active;
  arma::vec betaSD(beta.n_elem);
  int nr = Xstar.n_rows;
  int act;

  for(int g = 0; g < G; g++){
    active = find(groups==(g+1));
    if(active.size() == 1) {
      act = active(0);
      vec XstarAct = pow(Xstar, 2);
      double fi = sum(XstarAct.col(act));
      betaSD(act) = beta(act) * (sqrt((double)nr) / sqrt(fi));
    } else {
      vec Dvecg = Dvec[g];
      mat Qmatg = Qmat[g];
      betaSD.elem(active) = (Qmatg * diagmat(1/sqrt(Dvecg)) * beta.elem(active));
    }
  }
  return(betaSD);
}
