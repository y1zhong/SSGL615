#' A function to pergorm ssgl algorithm
#' @importFrom Rcpp sourceCpp
#' @import RcppArmadillo
#' @useDynLib SSGL615
#' @export
SSGL <- function(Y, X, lambda1, lambda0, groups,
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

 # library(Rcpp)
  #sourceCpp('src/c_Xstar.cpp')
  Xstar <- c_Xstar(X)

  #sourceCpp('src/Xtilde.cpp')
  Xtilde_list <- c_Xtilde(Xstar, groups, G, n)
  Xtilde <- Xtilde_list$Xnew
  Qmat <- Xtilde_list$Qmat
  Dvec <- Xtilde_list$Dvec

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
  l = list(beta = betaSD, betaStart = betaStart, theta=theta, sigmasqStart = sigmasqStart,
           sigmasq=sigmasq, intercept=interceptSD, nIter = counter,
           converged = converged)

  return(l)
}


#' @importFrom Rcpp sourceCpp
#' @import RcppArmadillo
#' @export
SSGLpath = function(Y, X, lambda1, lambda0, G,
                    lambda0seq = seq(lambda1, lambda0, length=20),groups,
                    a = 1, b = length(unique(groups)),
                    M = 10, error = 0.001,
                    forceGroups = c()) {

  ## Final model
  betaStart = rep(0, dim(X)[2])
  updateSigma=FALSE
  printWarnings = FALSE
  for (nl in 1 : length(lambda0seq)) {
    lambda0 = lambda0seq[nl]

    if (nl == length(lambda0seq)) printWarnings = TRUE

    # starting values for lambda0 = lambda1
    if ( nl == 1) {
      modSSGL = SSGL(Y=Y, X=X, lambda1=lambda1, lambda0=lambda0,
                     groups = groups,
                     a = 1, b = G,
                     updateSigma = updateSigma,
                     M = 10, error = error,
                     betaStart = betaStart,
                     theta = 0.5,
                     printWarnings = printWarnings
      )

    } else {
      modSSGL = SSGL(Y=Y, X=X, lambda1=lambda1, lambda0=lambda0,
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


