#include <Rcpp.h>
using namespace Rcpp;

//' @export
// [[Rcpp::export]]
NumericMatrix c_Xstar(NumericMatrix Xs){
  int nr = Xs.nrow();
  int nc = Xs.ncol();
  NumericMatrix Xnew( nr , nc );
  //for (j in 1 : dim(X)[2])
  for(int j=0; j < nc; j++){
    Xnew(_,j) = (Xs(_,j) - mean(Xs(_,j)));
  }
  return(Xnew);
}
