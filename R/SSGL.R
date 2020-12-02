#' @importFrom Rcpp sourceCpp
#' @import RcppArmadillo
SSGL <- function(Y, X, lambda1, lambda0, groups,
                 a = 1, b = length(unique(groups)),
                 updateSigma = TRUE,
                 M = 10, error = 0.001,
                 betaStart = rep(0, dim(X)[2]),
                 sigmasqStart,
                 theta,
                 printWarnings = TRUE,
                 forceGroups = c()) {

  n = 200
  G = 100
  set.seed(2020)
  x = mvtnorm::rmvnorm(n, sigma=diag(G))

  ## Now create matrix that has linear and squared functions## of the original covariates.
  X = matrix(NA, nrow=n, ncol=G*2)

  for (g in 1 : G) {
    X[,2*(g-1) + 1] = x[,g]
    X[,2*g] = x[,g]^2
  }
  set.seed(2020)
  Y = 200 + x[,1] + x[,2] + 0.6*x[,2]^2 + rnorm(n, sd=1)

  lambda1=.1
  lambda0=10
  groups = rep(1:G, each=2)

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


  sigmasqStart <- sqrt(df * ncp / (df + 2))


  library(Rcpp)
  sourceCpp('src/c_Xstar.cpp')
  system.time(Xstar <- c_Xstar(X))

  sourceCpp('src/Xtilde.cpp')
  system.time(Xtilde <- c_Xtilde(Xstar, groups, G, n))

  a = 1
  b = length(unique(groups))
  M = 10
  error = 0.001
  betaStart = rep(0, dim(X)[2])
  forceGroups = c()
  updateSigma = TRUE

  converged=TRUE

  ## Initialize values
  sigmasq = sigmasqStart
  beta = betaStart
  intercept = mean(Y - Xtilde %*% beta)
  Z = 1*((1:G) %in% unique(groups[which(beta != 0)]))
  #if (missing(theta)) {
    theta = (a + sum(Z))/(a + b + G)
  #}
  counter = 0
  diff = 10*error

  counter = 0

  lambda0_base = lambda0
  forceGroups = numeric()
  sourceCpp('src/update.cpp')

system.time(while(diff > error & counter < 300) {

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
 })
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

}
#require install gfortran
#sourceCpp('src/c_svd.cpp')
#sourceCpp('src/Xtilde.cpp')
#system.time(Xtil <- c_Xtilde(Xstar, groups, G, n))
#c_crossprod(tempX)
#my_fun(g, groups)
#a <- baseSVD((1/n) * t(tempX) %*% tempX)
#diag(1/sqrt(SVD$d))


src <-
  "NumericVector rcpp_Xstar(){
    NumericVector Xnew(7);
    Xnew = Range(0,6);
    return(Xnew);
  }
  "
Rcpp::cppFunction(src)

gFunc = function(beta, lambda1, lambda0, theta, sigmasq, n) {
  l = lambdaStar(beta = beta, lambda1 = lambda1,
                 lambda0 = lambda0, theta = theta)
  p = pStar(beta = beta, lambda1 = lambda1,
            lambda0 = lambda0, theta = theta)

  g = (l - lambda1)^2 + (2*n/sigmasq)*log(p)
  return(g)
}


psi = function(beta, lambda) {
  m = length(beta)
  C = 2^(-m) * pi^(-(m-1)/2) / (gamma((m+1)/2))
  logDens = log(C) + m*log(lambda) - lambda*sqrt(sum(beta^2))
  #dens = C * lambda^m * exp(-lambda*sqrt(sum(beta^2)))
  dens = exp(logDens)

  return(dens)
}

## pStar function
pStar = function(beta, lambda1, lambda0, theta) {
  psi1 = psi(beta=beta, lambda=lambda1)
  psi0 = psi(beta=beta, lambda=lambda0)

  ## if a coefficient is really large then both these will
  ## numerically be zero because R can't handle such small numbers
  if ((theta*psi1) == 0 & (1 - theta)*psi0 == 0) {
    p = 1
  } else {
    p = (theta*psi1) / (theta*psi1 + (1 - theta)*psi0)
  }

  return(p)
}

## Lambda star function
lambdaStar = function(beta, lambda1, lambda0, theta) {
  p = pStar(beta = beta, lambda1 = lambda1,
            lambda0 = lambda0, theta = theta)

  l = lambda1*p + lambda0*(1 - p)
  return(l)
}
