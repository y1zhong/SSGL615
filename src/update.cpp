// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
using namespace Rcpp;
using namespace arma;



//' find element from one vector in another vector
//' @return the target vector elements
// [[Rcpp::export]]
arma::uvec c_betainact(
    arma::vec beta, arma::uvec active) {

  int l = beta.size();
  uvec b0 = find(beta != 0)+1;
  uvec betalen = regspace<uvec>(1, l);
  NumericVector activevec = wrap(active);
  NumericVector betalenvec = wrap(betalen);
  NumericVector dif = setdiff(betalenvec, activevec);
  uvec different = as<uvec>(dif);
  uvec interFi = intersect(b0, different);
  return (interFi);
}


// [[Rcpp::export]]
double psi(NumericVector beta, double lambda){
  double m = beta.length();
  double pi = 3.141592653589793238463 ;
  double C = pow(2, (-m)) * pow(pi, (-(m-1)/2)) / (tgamma((m+1)/2));
  double logDens = log(C) + m*log(lambda) - lambda*sqrt(sum(pow(beta, 2)));
  double dens = exp(logDens);
  return dens;
}

// [[Rcpp::export]]
double pStar(NumericVector beta, double lambda1,
             double lambda0, double theta){
  double psi1 = psi(beta, lambda1);
  double psi0 = psi(beta, lambda0);
  double p;
  if ((theta*psi1) == 0 & (1 - theta)*psi0 == 0) {
    p = 1;
  } else {
    p = (theta*psi1) / (theta*psi1 + (1 - theta)*psi0);
  }
  return(p);
}

double lambdaStar(NumericVector beta, double lambda1,
                  double lambda0, double theta){
  double p = pStar(beta, lambda1, lambda0, theta);
  double l = lambda1*p + lambda0*(1 - p);
  return(l);
}
// [[Rcpp::export]]
double c_gFunc(NumericVector beta, double lambda1, double lambda0,
               double theta, double sigmasq, double n){
  double l = lambdaStar(beta, lambda1, lambda0, theta);
  double p = pStar(beta, lambda1, lambda0, theta);

  double g = pow((l - lambda1), 2) + (2*n/sigmasq)*log(p);
  return(g);

}

//' @export
// [[Rcpp::export]]
List update(arma::vec Y, arma::mat Xtilde, arma::vec groups, LogicalVector& updateSigma,
            double sigmasq, arma::vec betaa, double intercept, double lambda0_base,
            double lambda1, double lambda0, arma::vec betaOld, double a, double b, int M,
            LogicalVector Z, double theta, int G, IntegerVector forceGroups, double n){

  vec beta = betaa;
  vec yResid;

  int act2int;
  double diff;
  double delta;

  vec betaVecAct2;
  uvec active2uvec;
  uvec active2;
  uvec active;

  mat XtildeMatAct;
  mat XtildeMatAct2;
  mat zg;

  mat XtildeBeta;

  for(int g = 0; g < G; g++){

    active2 = find(beta!=0)+1;

    /////////////////////Intercept//////////////////////
    if (active2.size() == 0) {
      intercept = mean(Y);
    } else if (active2.size() == 1) {
      act2int = active2(0);
      intercept = mean(Y - Xtilde.col(act2int) * beta(act2int));
    } else {
      XtildeMatAct2 = Xtilde.cols(active2-1);
      betaVecAct2 = beta.elem(active2-1);
      XtildeBeta= XtildeMatAct2 * betaVecAct2;
      intercept = mean(Y - XtildeBeta.col(0));
    }

    // which parameters refer to this group
    active = find(groups==(g+1))+1;
    int m = active.size();
    lambda0 = sqrt(m) * lambda0_base;

    if ( std::find(forceGroups.begin(), forceGroups.end(), (g+1)) != forceGroups.end() ) {
      active2 = c_betainact(beta, active);
      if (active2.size() == 0) {
        yResid = Y - intercept;
      } else if (active2.size() == 1) {
        act2int = active2(0);
        yResid = Y - intercept - Xtilde.col(act2int) * beta(act2int);
      } else {
        XtildeMatAct2 = Xtilde.cols(active2-1);
        betaVecAct2 = beta.elem(active2-1);
        XtildeBeta= XtildeMatAct2 * betaVecAct2;
        yResid = Y - intercept - XtildeBeta.col(0);
      }

      XtildeMatAct = Xtilde.cols(active-1);
      mat XtildeMatActCrossInv = inv(XtildeMatAct.t() * XtildeMatAct);
      mat XildeMatActyResid = XtildeMatAct.t() * yResid;
      mat tempupdate = XtildeMatActCrossInv * XildeMatActyResid;
      vec betatobeupdate = tempupdate.col(0);
      beta.elem(active-1) = betatobeupdate;

    } else {

      NumericVector beta2(m);
      double gf = c_gFunc(beta2, lambda1, lambda0, theta, sigmasq, n);
      if (gf > 0) {
        double pstarr = pStar(beta2, lambda1, lambda0, theta);
        delta =  sqrt(2*n*sigmasq*log(1/pstarr)) + sigmasq*lambda1;
      } else {

        double lambdaStarr = lambdaStar(beta2, lambda1, lambda0, theta);
        delta = sigmasq*lambdaStarr;
      }

      //////////////////Calculate necessary quantities/////////////
      active2 = c_betainact(beta, active);

      vec tempyMinInter = Y - intercept;
      XtildeMatAct = Xtilde.cols(active-1);
      if (active2.size() == 0) {
       zg = XtildeMatAct.t() * tempyMinInter;
     } else if (active2.size()  == 1) {
        act2int = active2(0);
        zg = XtildeMatAct.t() * (tempyMinInter - XtildeMatAct2 * beta(act2int));
      } else {
        XtildeMatAct2 = Xtilde.cols(active2-1);
        betaVecAct2 = beta.elem(active2-1);
        zg = XtildeMatAct.t() * (tempyMinInter - XtildeMatAct2 * betaVecAct2);
      }

      mat zgsquare = pow(zg, 2);
      vec zgvec = zgsquare.col(0);

      double norm_zg = sqrt(sum(zgvec));
      vec betaVecAct = beta.elem(active-1);
      NumericVector betaactive = wrap(betaVecAct);
      double shrinkageLambda =  lambdaStar(betaactive, lambda1, lambda0, theta);
      double shrinkageTerm = (1/n) * (1 - sigmasq*shrinkageLambda/norm_zg);
      shrinkageTerm = shrinkageTerm*(1*(shrinkageTerm > 0));
      vec tobeupdate = shrinkageTerm*(zg.col(0))*(1*(norm_zg > delta));

      beta.elem(active-1) = (tobeupdate);
    }
    Z[g] = any(beta.elem(active-1) != 0);
    diff = sqrt(sum(pow((beta - betaOld),2)));

    if ((g+1) % M == 0) {

      if (forceGroups.length() == 0) {

        NumericVector tempZ = wrap(Z);
        theta = (a + sum(tempZ)) / (a + b + G);
      } else {
        IntegerVector zi = seq(0, Z.size() - 1);
        IntegerVector fi = forceGroups - 1;
        IntegerVector dif = setdiff(zi, fi);
        NumericVector tempZ = wrap(Z[dif]);
        theta = (a + sum(tempZ)) / (a + b + G - forceGroups.length());
      }

      if (updateSigma(0) == 1) {
        active2 = find(beta!=0)+1;
        if (active2.size() == 0) {
          sigmasq = sum(pow(Y - intercept,2)) / (n + 2);
        } else if (active2.size() == 1) {
          act2int = active2[0];
          sigmasq = sum(pow((Y - Xtilde.col(act2int)* beta(act2int) - intercept), 2)) / (n + 2);
        } else {
          XtildeMatAct2 = Xtilde.cols(active2-1);
          betaVecAct2 = beta.elem(active2-1);
          mat temp = XtildeMatAct2*betaVecAct2;
          sigmasq = sum(pow((Y - temp.col(0) - intercept), 2)) / (n + 2);
          }
        }
      }
  }
  List L = List::create(Named("sigmasq") = sigmasq , _["beta"] = beta,
                        _["intercept"] = intercept, _["diff"] = diff,
                        _["lambda1"] = lambda1, _["lambda0"] = lambda0,
                        _["theta"] = theta,_["Z"] = Z);
  return(L);
}
