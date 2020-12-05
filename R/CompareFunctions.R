#####Original Functions in R package
SSGL0 = function(Y, X, lambda1, lambda0, groups,
                 a = 1, b = length(unique(groups)),
                 updateSigma = TRUE,
                 M = 10, error = 0.001,
                 betaStart = rep(0, dim(X)[2]),
                 sigmasqStart,
                 theta,
                 printWarnings = TRUE,
                 forceGroups = c()) {

  ## Number of groups and covariates overall
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

  Xstar = matrix(NA, dim(X)[1], dim(X)[2])
  for (j in 1 : dim(X)[2]) {
    Xstar[,j] = (X[,j] - mean(X[,j]))
  }

  ## store orthonormalized design matrix
  Xtilde = matrix(NA, dim(X)[1], dim(X)[2])

  ## store relevant matrices from SVD within each group
  Qmat = list()
  Dvec = list()

  ## create orthonormal matrix
  for (g in 1 : G) {
    active = which(groups == g)

    if (length(active) == 1) {
      Xtilde[,active] = sqrt(dim(Xstar)[1]) * (Xstar[,active] / sqrt(sum(Xstar[,active]^2)))
      Qmat[[g]] = NULL
      Dvec[[g]] = NULL
    } else {

      tempX = Xstar[,active]
      SVD = svd((1/n) * t(tempX) %*% tempX)
      Qmat[[g]] = SVD$u
      Dvec[[g]] = SVD$d

      Xtilde[,active] = Xstar[,active] %*% Qmat[[g]] %*% diag(1/sqrt(Dvec[[g]]))
    }
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

  lambda0_base = lambda0

  while(diff > error & counter < 300) {

    ## Store an old beta so we can check for convergence at the end
    betaOld = beta

    ## Now update each group of parameters, beta_G
    for (g in 1 : G) {

      ## First update intercept
      active2 = which(beta != 0)
      if (length(active2) == 0) {
        intercept = mean(Y)
      } else if (length(active2) == 1) {
        intercept = mean(Y - Xtilde[,active2]*beta[active2])
      } else {
        intercept = mean(Y - Xtilde[,active2] %*% beta[active2])
      }

      ## which parameters refer to this group
      active = which(groups == g)
      m = length(active)
      lambda0 = sqrt(m) * lambda0_base

      if (g %in% forceGroups) {
        active2 = which(beta != 0 & !1:length(beta) %in% active)

        if (length(active2) == 0) {
          yResid = Y - intercept
        } else if (length(active2) == 1) {
          yResid = Y - intercept - Xtilde[,active2] * beta[active2]
        } else {
          yResid = Y - intercept - Xtilde[,active2] %*% as.matrix(beta[active2])
        }
        beta[active] = solve(t(Xtilde[,active]) %*% Xtilde[,active]) %*% t(Xtilde[,active]) %*% yResid
      } else {

        ## Calculate delta for this size of a group
        gf = gFunc(beta = rep(0, length(active)), lambda1 = lambda1,
                   lambda0 = lambda0, theta = theta, sigmasq = sigmasq, n=n)
        if (gf > 0) {
          delta =  sqrt(2*n*sigmasq*log(1/pStar(beta = rep(0,m), lambda1=lambda1,
                                                lambda0=lambda0, theta=theta))) +
            sigmasq*lambda1
        } else {
          delta = sigmasq*lambdaStar(beta = rep(0,m), lambda1=lambda1,
                                     lambda0=lambda0, theta=theta)
        }

        ## Calculate necessary quantities

        active2 = which(beta != 0 & !1:length(beta) %in% active)

        if (length(active2) == 0) {
          zg = t(Xtilde[,active]) %*% (Y - intercept)
        } else if (length(active2) == 1) {
          zg = t(Xtilde[,active]) %*% (Y - intercept - Xtilde[,active2] * beta[active2])
        } else {
          zg = t(Xtilde[,active]) %*% (Y - intercept - Xtilde[,active2] %*% as.matrix(beta[active2]))
        }

        norm_zg = sqrt(sum(zg^2))
        tempLambda = lambdaStar(beta = beta[active], lambda1 = lambda1,
                                lambda0 = lambda0, theta = theta)

        ## Update beta
        shrinkageTerm = (1/n) * (1 - sigmasq*lambdaStar(beta = beta[active], lambda1 = lambda1,
                                                        lambda0 = lambda0, theta = theta)/norm_zg)
        shrinkageTerm = shrinkageTerm*(1*(shrinkageTerm > 0))
        beta[active] = as.numeric(shrinkageTerm)*zg*as.numeric((1*(norm_zg > delta)))
      }


      ## Update Z
      Z[g] = 1*(any(beta[active] != 0))

      diff = sqrt(sum((beta - betaOld)^2))

      if (g %% M == 0) {
        ## Update theta
        if (length(forceGroups) == 0) {
          theta = (a + sum(Z)) / (a + b + G)
        } else {
          theta = (a + sum(Z[-forceGroups])) / (a + b + G - length(forceGroups))
        }

        ## Update sigmasq
        if (updateSigma) {
          active2 = which(beta != 0)
          if (length(active2) == 0) {
            sigmasq = sum((Y - intercept)^2) / (n + 2)
          } else if (length(active2) == 1) {
            sigmasq = sum((Y - Xtilde[,active2] * beta[active2] - intercept)^2) / (n + 2)
          } else {
            sigmasq = sum((Y - Xtilde[,active2] %*% as.matrix(beta[active2]) - intercept)^2) / (n + 2)
          }
        }
      }
    }

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

  ## need to re-standardize the variables
  betaSD = rep(NA, length(beta))
  for (g in 1 : G) {
    active = which(groups == g)
    if (length(active) == 1) {
      betaSD[active] = beta[active] * (sqrt(dim(Xstar)[1]) / sqrt(sum(Xstar[,active]^2)))
    } else {
      betaSD[active] = (Qmat[[g]] %*% diag(1/sqrt(Dvec[[g]])) %*% beta[active])
    }
  }

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

SSGLpath0 = function(Y, X, lambda1, lambda0,
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
      modSSGL = SSGL0(Y=Y, X=X, lambda1=lambda1, lambda0=lambda0,
                      groups = groups,
                      a = 1, b = G,
                      updateSigma = updateSigma,
                      M = 10, error = error,
                      betaStart = betaStart,
                      theta = 0.5,
                      printWarnings = printWarnings
      )

    } else {
      modSSGL = SSGL0(Y=Y, X=X, lambda1=lambda1, lambda0=lambda0,
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



#### C++ version
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
  Dvec <- lapply(Dvec,as.numeric)

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

SSGLpath = function(Y, X, lambda1, lambda0,
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


### C++ version with SVD adjusted
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

SSGLpath_onceSVD = function(Y, X, lambda1, lambda0,
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



