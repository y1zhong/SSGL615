#include <RcppArmadillo.h>
using namespace Rcpp;
using namespace arma;
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
List c_Xtilde(arma::mat Xstar, arma::vec groups, int G, int n){

  // store orthonormalized design matrix
  int nr = Xstar.n_rows;
  int nc = Xstar.n_cols;
  arma::mat Xnew( nr, nc );
  //arma::vec active;
  arma::uvec active;
  int act;
  //Rcout << "The value is " << nr << std::endl;
  // store relevant matrices from SVD within each group
  List Qmat;
  List Dvec;
  for(int g = 0; g < G; g++){
    //Rcout << "The g is " << g << std::endl;
    active = find(groups==(g+1));
    //active = as<vec>(activeuvec);
    //Rcout << "The value is " << active << std::endl;
    if(active.size() == 1) {
      //Xnew(_,1) = sum(Xstar(_,1)^2);
      //NumericVector nnew = sum(pow(Xstar(_,1), 2));
      act = active(0);
      vec nnew = Xstar.col(act);
      Xnew.col(act) = sqrt(nr * (nnew / sqrt(sum(pow(nnew,2)))));
      Qmat.insert(g, R_NilValue);
      Dvec.insert(g, R_NilValue);
    } else {
      mat tempXX = Xstar.cols(active);
      mat tempX = tempXX.t() * tempXX / n;
      //Rcout << "N is" << std::endl << n << std::endl;
      //Rcout << "tempXX is" << std::endl << tempXX << std::endl;
      //Rcout << "tempX is" << std::endl << tempX<< std::endl;
      //arma::mat mattempX = as<arma::mat>(tempX) ;
     // arma::mat mattempXX = as<arma::mat>(tempXX) ;
      //Rcout << "Armadillo matrix is" << std::endl << mattempX << std::endl;
      // Output matrices
      arma::mat U, V;
      arma::vec S;
      arma::svd(U, S, V, tempX, "standard");
      //arma::vec Sother = baseSVD(mattempX);
      Qmat.insert(g, U);
      Dvec.insert(g, S);

      //Rcout << "Armadillo matrix is" << std::endl << mattempX << std::endl;
      arma::mat repl = tempXX * U * diagmat(1/sqrt(S));
      //Rcout << "repl matrix is" << std::endl << repl << std::endl;
      //NumericMatrix replMtx = wrap(repl);
      Xnew.cols(active) = repl;
      //as<arma::mat>(tempX)
      //repl.diag(active.length(), 1/sqrt(S));
      //Qmat[[g]] = SVD$u
      //Dvec[[g]] = SVD$d
    }
  }
  List L = List::create(Named("Xnew") = Xnew , _["Qmat"] = Qmat,
                        _["Dvec"] = Dvec);
  return(L);
}
