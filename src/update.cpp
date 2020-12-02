#include <RcppArmadillo.h>
using namespace Rcpp;
using namespace arma;
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

// [[Rcpp::export]]
LogicalVector c_betainact(
    NumericVector beta, NumericVector active) {

  int l = beta.length();
  LogicalVector b0 = (beta != 0);
  NumericVector blen(l);
  blen = Range(1, l);
  LogicalVector out;
  for(unsigned int z = 0; z < l; z++) {
    if(std::find(active.begin(), active.end(), blen[z]) == active.end()){
      out.push_back(1);
    } else{
      out.push_back(0);
    }
  }
  return (b0 & out);
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


// [[Rcpp::export]]
List update(NumericVector Y, NumericMatrix Xtilde, NumericVector groups, LogicalVector updateSigma,
            double sigmasq, NumericVector beta, double intercept, double lambda0_base,
            double lambda1, double lambda0, NumericVector betaOld, double a, double b, int M,
            LogicalVector Z, double theta, int G, IntegerVector forceGroups, double n){
  //Rcout << "The begin " << std::endl;
  //NumericMatrix Ytemp( 2 , 3 , v.begin() );
  //mat Ymat = as<mat>(Y);
  //Rcout << "The Y is " << Ymat << std::endl;
  beta = clone(beta);
  NumericVector yResid;
  //mat Xtildemat = as<mat>(Xtilde);
  double diff;

  for(int g = 0; g < G; g++){
    //Rcout << "The g is " << g << std::endl;
    NumericVector active2 = c_which2(beta!=0);
    int act2 = active2[0];
    vec betaAct2Vec = as<vec>(vector_subset_idx_rcpp2(beta, active2-1));
    mat XtildeAct2Mat = as<mat>(matrix_subset_idx_rcpp2(Xtilde, active2-1));
    NumericMatrix tempinter= wrap(XtildeAct2Mat * betaAct2Vec);

    /////////////////////Intercept//////////////////////
    //Rcout << "The active2 is " << active2 << std::endl;
    if (active2.length() == 0) {
      intercept = mean(Y);
    } else if (active2.length() == 1) {
      intercept = mean(Y - Xtilde(_,act2) * beta[act2]);
    } else {
      //vec tempXX = as<vec>(vector_subset_idx_rcpp2(beta, active2-1));
      //mat tempXtilde = as<mat>(matrix_subset_idx_rcpp2(Xtilde, active2-1));
      //NumericVector tempXX = vector_subset_idx_rcpp2(beta, active2-1);
      //NumericMatrix tempXtilde = matrix_subset_idx_rcpp2(Xtilde, active2-1);
      //NumericMatrix tempinter= wrap(tempXtilde * tempXX);
      //NumericMatrix tempinter= wrap(XtildeAct2Mat * betaAct2Vec);
      intercept = mean(Y - tempinter(_,0));
    }

    // which parameters refer to this group
    NumericVector active = c_which2(groups==(g+1));
    //Rcout << "The active is " << active << std::endl;
    int m = active.length();
    lambda0 = sqrt(m) * lambda0_base;
    //Rcout << "The lambda0 is " << lambda0 << std::endl;



    if ( std::find(forceGroups.begin(), forceGroups.end(), (g+1)) != forceGroups.end() ) {
      LogicalVector v = c_betainact(beta, active);
      active2 = c_which2(v);
      if (active2.length() == 0) {
        yResid = Y - intercept;
      } else if (active2.length() == 1) {
        act2 = active2[0];
        yResid = Y - intercept - Xtilde(_,act2) * beta[act2];
      } else {
        //vec tempXX = as<vec>(vector_subset_idx_rcpp2(beta, active2-1));
       // mat tempXtilde = as<mat>(matrix_subset_idx_rcpp2(Xtilde, active2-1));
       // mat XXmul = tempXX * tempXtilde;
       // NumericMatrix XXmulmtx = wrap(XXmul);
        yResid = Y - intercept - tempinter(_,0);
      }
      //Rcout << "update Beta1 " << std::endl;
      mat tempXtilde = as<mat>(matrix_subset_idx_rcpp2(Xtilde, active-1));
      //Rcout << "tempXtildecrossp " << std::endl;
      mat tempXtildecrossp = tempXtilde.t() * tempXtilde;
     // Rcout << "tempXtildeyRes" << std::endl;
      mat tempXtildeyRes = tempXtilde.t() * as<vec>(yResid);
     //Rcout << "tempupdate" << std::endl;
      NumericMatrix tempupdate = wrap(tempXtildecrossp * tempXtildeyRes);
      NumericVector betatobeupdate = tempupdate(_,0);
      beta = vector_assign_rcpp2(beta, betatobeupdate, active-1);
    } else {
      //Calculate delta for this size of a group
      //Rcout << "call gFunc " << std::endl;
      NumericVector beta2(active.length());
      double gf = c_gFunc(beta2, lambda1, lambda0, theta, sigmasq, n);
      //Rcout << "The gf is " << gf << std::endl;
      double delta;
      NumericVector beta3(m);
      if (gf > 0) {
        double pstarr = pStar(beta3, lambda1, lambda0, theta);
        delta =  sqrt(2*n*sigmasq*log(1/pstarr)) + sigmasq*lambda1;
      } else {
       // Rcout << "The sigmasq is " << sigmasq<< std::endl;
        double lambdaStarr = lambdaStar(beta3, lambda1, lambda0, theta);
       // Rcout << "The lambdaStarr is " << lambdaStarr << std::endl;
        delta = sigmasq*lambdaStarr;
      }

      //////////////////Calculate necessary quantities/////////////
     // Rcout << "The delta is " << delta << std::endl;
      LogicalVector v = c_betainact(beta, active);
      active2 = c_which2(v);
      //Rcout << "The active2 is " << active2 << std::endl;
      NumericVector yMinInter = Y - intercept;
      //Rcout << "calculate yMinInter " << yMinInter <<std::endl;
      vec tempyMinInter = as<vec>(yMinInter);
      mat XtildeActMat = as<mat>(matrix_subset_idx_rcpp2(Xtilde, active-1));
      //mat tempXtilde2 = as<mat>(matrix_subset_idx_rcpp2(Xtilde, active2-1));
      //mat tempbeta = as<mat>(beta[active2]);

      mat zg;
      if (active2.length() == 0) {
       // Rcout << "calculate tempXtilde" << tempXtilde <<std::endl;

        zg = XtildeActMat.t() * tempyMinInter;
       // zg = t(Xtilde[,active]) %*% (Y - intercept)
      } else if (active2.length() == 1) {
        act2 = active2[0];
        zg = XtildeActMat.t() * (tempyMinInter - XtildeAct2Mat * beta[act2]);
       // zg = t(Xtilde[,active]) %*% (Y - intercept - Xtilde[,active2] * beta[active2])
      } else {
        //vec tempbeta = as<vec>(beta[active2-1]);
        // Rcout << "calculate tempyMinInter" << tempyMinInter <<std::endl;
       //  Rcout << "calculate XtildeAct2Mat" << XtildeAct2Mat<<std::endl;
       //  Rcout << "betaAct2Vec"  <<betaAct2Vec<<std::endl;
         betaAct2Vec = as<vec>(vector_subset_idx_rcpp2(beta, active2-1));
         XtildeAct2Mat = as<mat>(matrix_subset_idx_rcpp2(Xtilde, active2-1));
        zg = XtildeActMat.t() * (tempyMinInter - XtildeAct2Mat * betaAct2Vec);
       // zg = t(Xtilde[,active]) %*% (Y - intercept - Xtilde[,active2] %*% as.matrix(beta[active2]))
      }
     // Rcout << "calculate zg " << zg <<std::endl;
      NumericMatrix tempzg = wrap(zg);
      NumericVector zgvec = tempzg(_,0);
   //   Rcout << "calculate zg " << zgvec <<std::endl;
      double norm_zg = sqrt(sum(pow(zgvec, 2)));
     // Rcout << "calculate norm_zg " << norm_zg <<std::endl;
      NumericVector betaactive = vector_subset_idx_rcpp2(beta, active-1);
      //double tempLambda = lambdaStar(betaactive, lambda1, lambda0, theta);
      //Rcout << "update beta2 " << std::endl;
      double shrinkageLambda =  lambdaStar(betaactive, lambda1, lambda0, theta);
      double shrinkageTerm = (1/n) * (1 - sigmasq*shrinkageLambda/norm_zg);
      shrinkageTerm = shrinkageTerm*(1*(shrinkageTerm > 0));

      NumericVector tobeupdate = shrinkageTerm*zgvec*(1*(norm_zg > delta));
      //Rcout << " norm_zg " << norm_zg << std::endl;
     // Rcout << " delta " << delta << std::endl;
     // Rcout << "tobeupdate " << tobeupdate << std::endl;
      beta = vector_assign_rcpp2(beta, tobeupdate, active-1);
    }
   // Rcout << "update Z " << std::endl;
   //Rcout << "active beyond" << active << std::endl;
    NumericVector betaactive2 = vector_subset_idx_rcpp2(beta, active-1);
    //Rcout << "betaactive2 " << betaactive2 << std::endl;
    //Rcout << "bool is " << as<bool>(any(betaactive2 != 0)) << std::endl;
    Z[g] = as<bool>(any(betaactive2 != 0));
    //Z[g] = tempBool[0];
    //Rcout << "beta is " << beta << std::endl;
    //Rcout << "betaOld is " << betaOld << std::endl;
   // NumericVector bbb = beta - betaOld;
    diff = sqrt(sum(pow((beta - betaOld),2)));
   // Rcout << "beta is " << beta << std::endl;
    //Rcout << "betaOld is " << betaOld << std::endl;
  //  Rcout << "diff is " << diff << std::endl;

    if ((g+1) % M == 0) {
      //Rcout << "g%m " << ((g+1) % M == 0 )<< std::endl;
  // Update theta
 // Rcout << "Update theta " << std::endl;
      if (forceGroups.length() == 0) {

        NumericVector tempZ = wrap(Z);
        //Rcout << "tempZ " << tempZ << std::endl;
        theta = (a + sum(tempZ)) / (a + b + G);
      } else {
        IntegerVector zi = seq(0, Z.size() - 1);
        IntegerVector fi = forceGroups - 1;
        IntegerVector dif = setdiff(zi, fi);
        NumericVector tempZ = wrap(Z[dif]);
        theta = (a + sum(tempZ)) / (a + b + G - forceGroups.length());
      }

     // Rcout << "Update sigmasq " << std::endl;
      if (updateSigma) {
        active2 = c_which2(beta != 0);
        if (active2.length() == 0) {
          sigmasq = sum(pow(Y - intercept,2)) / (n + 2);
        } else if (active2.length() == 1) {
          act2 = active2[0];
          sigmasq = sum(pow((Y - Xtilde(_,act2) * beta[act2] - intercept), 2)) / (n + 2);
        } else {
          mat tempXtilde = as<mat>(matrix_subset_idx_rcpp2(Xtilde, active2-1));
          vec betaactive333 = as<vec>(vector_subset_idx_rcpp2(beta, active2-1));
          //vec tempbeta = as<vec>(betaactive333);
          mat tempMul = tempXtilde * betaactive333;
          NumericMatrix temp = wrap(tempMul);
          sigmasq = sum(pow((Y - temp(_,0) - intercept), 2)) / (n + 2);
         //  Rcout << "sigmasq " << sigmasq << std::endl;
        }

      }
    }
   // Rcout << "=========== " << std::endl;
  // if(g==0) break;
  }
  //List L = List::create(3, 4);
  List L = List::create(Named("sigmasq") = sigmasq , _["beta"] = beta,
                        _["intercept"] = intercept, _["diff"] = diff,
                        _["lambda1"] = lambda1, _["lambda0"] = lambda0,
                        _["theta"] = theta,_["Z"] = Z);
  return(L);
}
