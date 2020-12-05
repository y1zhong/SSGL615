# n = 200
# G = 100
# set.seed(2020)
# x = mvtnorm::rmvnorm(n, sigma=diag(G))
#
# ## Now create matrix that has linear and squared functions## of the original covariates.
# X = matrix(NA, nrow=n, ncol=G*2)
#
# for (g in 1 : G) {
#   X[,2*(g-1) + 1] = x[,g]
#   X[,2*g] = x[,g]^2
# }
# set.seed(2020)
# Y = 200 + x[,1] + x[,2] + 0.6*x[,2]^2 + rnorm(n, sd=1)
#
# ## Now fit model for chosen lambda0 and lambda1 values
# modSSGL = SSGL(Y=Y, X=X, lambda1=.1, lambda0=10,
#                groups = rep(1:G, each=2), updateSigma = F)
#
#
