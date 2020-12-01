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
  set.seed(123)
  x = mvtnorm::rmvnorm(n, sigma=diag(G))

  ## Now create matrix that has linear and squared functions## of the original covariates.
  X = matrix(NA, nrow=n, ncol=G*2)

  for (g in 1 : G) {
    X[,2*(g-1) + 1] = x[,g]
    X[,2*g] = x[,g]^2
  }

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


  # Xstar_old = matrix(NA, dim(X)[1], dim(X)[2])
  # system.time(
  # for (j in 1 : dim(X)[2]) {
  #   Xstar_old[,j] = (X[,j] - mean(X[,j]))
  # })

  # src <-
  #   "NumericMatrix rcpp_Xstar(NumericMatrix Xs){
  #     int nr = Xs.nrow();
  #     int nc = Xs.ncol();
  #     NumericMatrix Xnew( nr , nc );
  #     for(int j=0; j < nc; j++){
  #       Xnew(_,j) = (Xs(_,j) - mean(Xs(_,j)));
  #     }
  #     return(Xnew);
  #   }
  #   "
  #Rcpp::cppFunction(src)
  library(Rcpp)
  sourceCpp('src/c_Xstar.cpp')
  system.time(Xstar <- c_Xstar(X))


  # ## store orthonormalized design matrix
  # Xtilde = matrix(NA, dim(X)[1], dim(X)[2])
  #
  # ## store relevant matrices from SVD within each group
  # Qmat = list()
  # Dvec = list()

# system.time(
#   for (g in 1 : G) {
#     active = which(groups == g)
#
#     if (length(active) == 1) {
#       Xtilde[,active] = sqrt(dim(Xstar)[1]) * (Xstar[,active] / sqrt(sum(Xstar[,active]^2)))
#       Qmat[[g]] = NULL
#       Dvec[[g]] = NULL
#     } else {
#
#       tempX = Xstar[,active]
#       SVD = svd((1/n) * t(tempX) %*% tempX)
#       Qmat[[g]] = SVD$u
#       Dvec[[g]] = SVD$d
#
#       Xtilde[,active] = Xstar[,active] %*% Qmat[[g]] %*% diag(1/sqrt(Dvec[[g]]))
#     }
#   }
# )
  sourceCpp('src/Xtilde.cpp')
  system.time(Xtil <- c_Xtilde(Xstar, groups, G, n))

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
