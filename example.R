source("ag1.R")

set.seed(10)
mu1 = c(0,0)
mu2 = c(2,0)
mu3 = c(0,5)
mu4 = c(2,5)

sigma0 = 0.5

n1 = 1000
n2 = 1000
n3 = 200
n4 = 200

X1 = cbind(rnorm(n1, mean=mu1[1], sd=sigma0),rnorm(n1, mean=mu1[2], sd=sigma0))
X2 = cbind(rnorm(n2, mean=mu2[1], sd=sigma0),rnorm(n2, mean=mu2[2], sd=sigma0))
X3 = cbind(rnorm(n3, mean=mu3[1], sd=sigma0),rnorm(n3, mean=mu3[2], sd=sigma0))
X4 = cbind(rnorm(n4, mean=mu4[1], sd=sigma0),rnorm(n4, mean=mu4[2], sd=sigma0))

X = rbind(X1,X2,X3,X4)

h = h.dD(X)

x_test = GD_flow2(X, X, h)
