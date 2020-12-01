//#include <Rcpp.h>
#include <RcppArmadillo.h>
using namespace Rcpp;
using namespace arma;
// [[Rcpp::depends(RcppArmadillo)]]

//' @export
// [[Rcpp::export]]
arma::vec baseSVD(const arma::mat & X) {
  arma::mat U, V;
  arma::vec S;
  arma::svd(U, S, V, X, "standard");
  return S;
}

// [[Rcpp::export]]
NumericVector c_which(NumericVector groups, int g){
  // calling rnorm()
  Function f("which");

  // Next code is interpreted as rnorm(n=5, mean=10, sd=2)
  return f(groups == g);
}
// [[Rcpp::export]]
NumericMatrix c_crossprod(NumericMatrix X){
  // calling rnorm()
  Function f("crossprod");

  // Next code is interpreted as rnorm(n=5, mean=10, sd=2)
  return f(X);
}
// https://stackoverflow.com/questions/62118084/rcpp-select-subset-numericmatrix-column-by-a-numericvector
// [[Rcpp::export]]
Rcpp::NumericMatrix matrix_subset_idx_rcpp(
    Rcpp::NumericMatrix x, Rcpp::NumericVector y) {

  // Determine the number of observations
  int n_cols_out = y.size();

  // Create an output matrix
  Rcpp::NumericMatrix out = Rcpp::no_init(x.nrow(), n_cols_out);

  // Loop through each column and copy the data.
  for(unsigned int z = 0; z < n_cols_out; ++z) {
    out(Rcpp::_, z) = x(Rcpp::_, y[z]);
  }

  return out;
}

NumericMatrix matrix_assign_rcpp(
    NumericMatrix X, NumericMatrix Y, NumericVector id) {

  int l = id.length();


  for(unsigned int z = 0; z < l; z++) {
    int pos = id[z];
    NumericMatrix::Column col = X( _ , pos);
    col = Y(_,z);
  }
  return X;
}



//' @export
// [[Rcpp::export]]
NumericMatrix c_Xtilde(NumericMatrix Xstar, NumericVector groups, int G, int n){

  // store orthonormalized design matrix
  int nr = Xstar.nrow();
  int nc = Xstar.ncol();
  NumericMatrix Xnew( nr , nc );
  //Rcout << "The value is " << nr << std::endl;
  // store relevant matrices from SVD within each group
  List Qmat;
  List Dvec;
  for(int g = 0; g < G; g++){
    //Rcout << "The g is " << g << std::endl;
    NumericVector active = c_which(groups, g+1);
    //Rcout << "The value is " << active << std::endl;
    if(active.length() == 1) {
      //Xnew(_,1) = sum(Xstar(_,1)^2);
      //NumericVector nnew = sum(pow(Xstar(_,1), 2));
      int act = active[0];
      NumericVector nnew = Xstar(_,act);
      Xnew(_,act) = sqrt(nr * (nnew / sqrt(sum(pow(nnew,2)))));
      Qmat.insert(g, R_NilValue);
      Dvec.insert(g, R_NilValue);
    } else {
      NumericMatrix tempXX = matrix_subset_idx_rcpp(Xstar, active-1);
      NumericMatrix tempX = c_crossprod(tempXX) / n;
      //Rcout << "N is" << std::endl << n << std::endl;
      //Rcout << "tempXX is" << std::endl << tempXX << std::endl;
      //Rcout << "tempX is" << std::endl << tempX<< std::endl;
      arma::mat mattempX = as<arma::mat>(tempX) ;
      arma::mat mattempXX = as<arma::mat>(tempXX) ;
      //Rcout << "Armadillo matrix is" << std::endl << mattempX << std::endl;
      // Output matrices
      arma::mat U, V;
      arma::vec S;
      arma::svd(U, S, V, mattempX, "standard");
      //arma::vec Sother = baseSVD(mattempX);
      Qmat.insert(g, U);
      Dvec.insert(g, S);
      //NumericMatrix repl(active.length());
      //NumericVector a(3,4);
      //repl.fill_diag(a);
      //Rcout << "Armadillo matrix is" << std::endl << U << std::endl;
      //Rcout << "Armadillo matrix is" << std::endl << mattempX << std::endl;
      arma::mat repl = mattempXX * U * diagmat(1/sqrt(S));
      //Rcout << "repl matrix is" << std::endl << repl << std::endl;
      NumericMatrix replMtx = wrap(repl);
      Xnew = matrix_assign_rcpp(Xnew, replMtx, active-1);
      //as<arma::mat>(tempX)
      //repl.diag(active.length(), 1/sqrt(S));
      //Qmat[[g]] = SVD$u
      //Dvec[[g]] = SVD$d
    }
  }
  return(Xnew);
}
