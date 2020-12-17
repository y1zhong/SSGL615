#' ### C++ version with SVD adjusted
#' A function to pergorm fast ssgl algorithm
#'
#' @param Y a vector of response variable
#' @param X a matrix of predictors
#' @param lambda1 inital value of lambda1
#' @param lambda0 inital value of lambda0
#' @param lambda0seq sequence from lambda0 to lambda1
#' @param a,b default parameter setup
#' @param M default is 10
#' @param error numeric value for tolerance
#' @param forceGroups a vector indicating index to force into same group
#'
#' @importFrom Rcpp sourceCpp
#' @import RcppArmadillo
#' @export
#' @return A list of results from ssgl path algorithm
SSGLpath_onceSVD = function(Y, X, lambda1, lambda0,G,
                            lambda0seq = seq(lambda1, lambda0, length=20),groups,
                            a = 1, b = length(unique(groups)),
                            M = 10, error = 0.001,
                            forceGroups = c()) {

  ## Final model
  betaStart = rep(0, dim(X)[2])
  updateSigma=FALSE
  printWarnings = FALSE


  ## pre-standardize matrix
  Xstar <- c_Xstar(X)
  Xtilde_list <- c_Xtilde(Xstar, groups, G, n)
  Xtilde <- Xtilde_list$Xnew
  Qmat <- Xtilde_list$Qmat
  Dvec <- Xtilde_list$Dvec



  for (nl in 1 : length(lambda0seq)) {
    lambda0 = lambda0seq[nl]

    if (nl == length(lambda0seq)) printWarnings = TRUE

    # starting values for lambda0 = lambda1
    if ( nl == 1) {
      modSSGL = SSGL_noSVD(Y=Y, X=X, Xstar=Xstar, Xtilde=Xtilde,
                           Qmat=Qmat, Dvec=Dvec,
                           lambda1=lambda1, lambda0=lambda0,
                           groups = groups,
                           a = 1, b = G,
                           updateSigma = updateSigma,
                           M = 10, error = error,
                           betaStart = betaStart,
                           theta = 0.5,
                           printWarnings = printWarnings
      )

    } else {
      modSSGL = SSGL_noSVD(Y=Y, X=X, Xstar=Xstar, Xtilde=Xtilde,
                           Qmat=Qmat, Dvec=Dvec,
                           lambda1=lambda1, lambda0=lambda0,
                           groups = groups,
                           a = 1, b = G,
                           updateSigma = updateSigma,
                           M = 10, error = error,
                           betaStart = betaStart,
                           sigmasqStart = sigmasqStart,
                           printWarnings = printWarnings)

    }


    if (modSSGL$nIter < 100 & modSSGL$converged) {
      updateSigma = TRUE
    }

    betaStart = modSSGL$betaStart
    sigmasqStart = modSSGL$sigmasqStart

  }

  l = list(beta = modSSGL$beta, theta=modSSGL$theta,
           sigmasq=modSSGL$sigmasq, intercept=modSSGL$intercept)

  return(l)
}
