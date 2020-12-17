### C++ version with SVD adjusted
#' A function to pergorm fast ssgl algorithm
#'
#' @param Y a vector of response variable
#' @param X a matrix of predictors
#' @param lambda0 inital value of lambda
#' @param groups a vector indicting the group number
#' @param a,b default parameter setup
#' @param updateSigma T/F indicating to update sigma value
#' @param M default is 10
#' @param error numeric value for tolerance
#' @param betaStart vector default set to all 0
#' @param sigmasqStart updated value through loop
#' @param theta initial theta value
#' @param forceGroups a vector indicating index to force into same group
#'
#' @importFrom Rcpp sourceCpp
#' @import RcppArmadillo
#' @export
#' @return A list of results from ssgl algorithm
SSGL_noSVD <- function(Y, X, Xstar, Xtilde, Qmat, Dvec, lambda1, lambda0, groups,
                       a = 1, b = length(unique(groups)),
                       updateSigma = TRUE,
                       M = 10, error = 0.001,
                       betaStart = rep(0, dim(X)[2]),
                       sigmasqStart,
                       theta,
                       printWarnings = TRUE,
                       forceGroups = numeric()) {

  G = length(unique(groups))
  p = length(groups)
  n = length(Y)


  # get initial value for sigma
  df <- 3
  sigquant <- 0.9
  sigest <- sd(Y)
  qchi <- qchisq(1 - sigquant, df)
  ncp <- sigest^2 * qchi / df
  min_sigma2 <- sigest^2 / n

  if (missing(sigmasqStart)) {
    sigmasqStart <- sqrt(df * ncp / (df + 2))
  }


  converged=TRUE

  ## Initialize values
  sigmasq = sigmasqStart
  beta = betaStart
  intercept = mean(Y - Xtilde %*% beta)
  Z = 1*((1:G) %in% unique(groups[which(beta != 0)]))
  if (missing(theta)) {
    theta = (a + sum(Z))/(a + b + G)
  }
  counter = 0
  diff = 10*error

  counter = 0

  lambda0_base = lambda0
  #forceGroups = numeric()
  #sourceCpp('src/update.cpp')

  while(diff > error & counter < 300) {

    ## Store an old beta so we can check for convergence at the end
    # betaOld = beta
    betaOld = beta

    l_update <- update( Y, Xtilde, groups, updateSigma,
                        sigmasq, beta,  intercept, lambda0_base,
                        lambda1, lambda0, betaOld, a, b, M,
                        Z, theta, G, forceGroups, n)


    beta = l_update$beta
    intercept = l_update$intercept
    sigmasq = l_update$sigmasq
    diff = l_update$diff
    lambda0 = l_update$lambda0
    theta = l_update$theta
    Z = l_update$Z
    #cat('beta:', beta)
    ## Check to make sure algorithm doesn't explode for values of lambda0 too small
    active2 = which(beta != 0)
    if (length(active2) == 0) {
      tempSigSq = sum((Y - intercept)^2) / (n + 2)
    } else if (length(active2) == 1) {
      tempSigSq = sum((Y - Xtilde[,active2] * beta[active2] - intercept)^2) / (n + 2)
    } else {
      tempSigSq = sum((Y - Xtilde[,active2] %*% as.matrix(beta[active2]) - intercept)^2) / (n + 2)
    }


    if (updateSigma==FALSE & (tempSigSq < min_sigma2 | tempSigSq > 100*var(Y))) {
      if (printWarnings == TRUE) {
        print(paste("lambda0 = ", lambda0, ",", " Algorithm diverging. Increase lambda0 or lambda0seq", sep=""))
      }
      sigmasq = sigmasqStart
      betaStart = rep(0, dim(Xtilde)[2])
      converged = FALSE
      diff = 0
    }

    if (updateSigma==TRUE & (tempSigSq < min_sigma2 | tempSigSq > 100*var(Y))) {
      sigmasq = sigmasqStart
      beta = betaStart
      updateSigma = FALSE
      diff = 10*error
      counter = 0
    }

    if (sum(beta != 0) >= min(n, p)) {
      if (printWarnings == TRUE) {
        print("Beta is saturated. Increase lambda0 or lambda0seq")
      }
      sigmasq = sigmasqStart
      betaStart = rep(0, dim(Xtilde)[2])
      converged = F
      break
    }

    counter = counter + 1
  }


  betaSD = c_betaSD(Xstar, groups, beta, Qmat, Dvec, G, n)

  ## Update new intercept
  interceptSD = mean(Y - X %*% betaSD)

  ## update the starting value for next iteration only if model converged
  if (converged == TRUE) {
    betaStart = beta
    if (updateSigma) {
      sigmasqStart = sigmasq
    }
  }

  ## estimate sigma regardless of convergence
  active2 = which(betaSD != 0)
  if (length(active2) == 0) {
    sigmasq = sum((Y - interceptSD)^2) / (n + 2)
  } else if (length(active2) == 1) {
    sigmasq = sum((Y - X[,active2] * betaSD[active2] - interceptSD)^2) / (n + 2)
  } else {
    sigmasq = sum((Y - X[,active2] %*% as.matrix(betaSD[active2]) - interceptSD)^2) / (n + 2)
  }


  l = list(beta = as.numeric(betaSD), betaStart = as.numeric(betaStart), theta=theta, sigmasqStart = sigmasqStart,
           sigmasq=sigmasq, intercept=interceptSD, nIter = counter,
           converged = converged)

  return(l)
}
