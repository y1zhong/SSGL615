// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
using namespace Rcpp;
using namespace arma;

//' Parametrize X matrix
//'
//' @param Xstar input centered matrix
//' @param groups group index
//' @param G the number of unique groups
//' @param n length of outcome variable
//' @return A list with output matrix and SVD calculation
// [[Rcpp::export]]
List c_Xtilde(arma::mat Xstar, arma::vec groups, int G, int n){

  // store orthonormalized design matrix
  int nr = Xstar.n_rows;
  int nc = Xstar.n_cols;
  arma::mat Xnew( nr, nc );
  arma::uvec active;
  int act;
  List Qmat;
  List Dvec;
  for(int g = 0; g < G; g++){

    active = find(groups==(g+1));

    if(active.size() == 1) {
      act = active(0);
      vec nnew = Xstar.col(act);
      Xnew.col(act) = sqrt(nr * (nnew / sqrt(sum(pow(nnew,2)))));
      Qmat.insert(g, R_NilValue);
      Dvec.insert(g, R_NilValue);
    } else {
      mat tempXX = Xstar.cols(active);
      mat tempX = tempXX.t() * tempXX / n;
      // Output matrices
      arma::mat U, V;
      arma::vec S;
      arma::svd(U, S, V, tempX, "standard");
      Qmat.insert(g, U);
      Dvec.insert(g, S);

      arma::mat repl = tempXX * U * diagmat(1/sqrt(S));
      Xnew.cols(active) = repl;
    }
  }
  List L = List::create(Named("Xnew") = Xnew , _["Qmat"] = Qmat,
                        _["Dvec"] = Dvec);
  return(L);
}
