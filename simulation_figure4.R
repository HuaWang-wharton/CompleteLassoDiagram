#########################################################################################
## This R code is used for plotting Figure 4 (left) in paper "The Complete Lasso Diagram"
##---------------------------------------------------------------------------------------
## Copyright @ Hua Wang, Yachong Yang, Zhiqi Bu, and Weijie Su, 2020
##---------------------------------------------------------------------------------------
#########################################################################################

library(MASS)
library(glmnet)
library(ggplot2)
source('q_lower.R')
library(RColorBrewer)

#########################################################################################
n=1000 ## can be changed to 700 and 500 to obtain other plots.
p=1000
s=300
total = 100
#sim_name = "VarySigma_n1000p1000s300"
##Set X's baselines
X_raw = 0
## Set total number of different configurations as N
N = 8
##Set signals
beta0=rep(0,p)
beta0[1:s]= c(0.01, 0.1, 1, 10, 100)
betas = list(beta0, beta0, beta0, beta0, beta0, beta0, beta0, beta0)
lseq = function(from, to, length.out) {
  exp(seq(log(from), log(to), length.out = length.out))
}
##Set lambdas
lambda1 <- c(lseq(0.00005, 0.01, 50), lseq(0.01, 10, length = 100), lseq(1, 200, length = 50))
lambda2 <- c(lseq(0.0005, 0.01, 30), lseq(0.01, 10, length = 100), lseq(1, 200, length = 70))

## Set different noise levels
sigmas = list(0, 0.001, 0.01, 0.1, 1, 10, 100, 1000)
lambdas = list()
for (i in 1:5) {
  lambdas[[i]] = lambda1 * (1 + (sigmas[[i]]))
}
for (i in 6:N) {
  lambdas[[i]] = lambda2 * (1 + (sigmas[[i]]))
}

########################################################################################
fdpps = list()
tppps = list()
for (i in 1:N){
  fdpps[[i]] = matrix(0, total,length(lambda1))
  tppps[[i]] = matrix(0, total,length(lambda1))
}

##essential func
fdp_tpp_lasso <- function(X, beta, eps, lambda){
  ##X is matrix, beta, eps, lambda are vectors
  Y=X%*%beta + eps
  n = length(Y)
  p = dim(X)[2]
  fit<-glmnet(X,Y,intercept=F, lambda = lambda, thresh = 1e-8)
  betahat <- coef(fit, exact = T)[-1, ]
  discov <- apply(betahat!=0, 2, sum) ## This is wrong, since glmnet is naive, and can give more than n discoveries.
  discov <- pmin(n, discov) ## Cap with n
  if (n >= p) {
    falsediscov <- apply(betahat[-(1:s), ]!=0, 2, sum)
  } else {
    falsediscov = rep(0, dim(betahat)[2])
    for (i in 1:dim(betahat)[2]) {
      thre = sort(abs(betahat[, i]))[n+1] # This is >= 0
      falsediscov[i] = sum(abs(betahat[-(1:s),i]) > thre)
    }
  }
  fdp_vec=falsediscov/discov
  tpp_vec=(discov-falsediscov)/s
  return(list(fdp_vec, tpp_vec))
}


##Simulation
for (k in 1:total){
  print(k)
  X=mvrnorm(n,rep(0,p),diag(1/n,p))
  X <- scale(X_raw + X)/sqrt(n-1)
  for (i in 1:N){
    eps = rnorm(n, 0, sigmas[[i]]) # with i-th level of noise.
    ft = fdp_tpp_lasso(X, betas[[i]], eps, lambdas[[i]])
    fdpps[[i]][k,] = ft[[1]]
    tppps[[i]][k,] = ft[[2]]
    print(paste(i, "th config is done."))
  }
}

##Get fdp and tpp
fdp = list()
tpp = list()
for (i in 1:N) {
  fdp[[i]] = apply(fdpps[[i]], 2, mean)
  tpp[[i]] = apply(tppps[[i]], 2, mean)
}

save(tpp, fdp, file = paste0(sim_name,".RData"))

########################################################################################

## plot
#load(file = paste0(sim_name,".RData"))
source("q_lower.R")
## Get the data of upper boundary
delta = n / p
eps = s / p
L=1000
tppmax = powermax(delta, eps)
xu = seq(0, tppmax, length.out = L)
yu = rep(0,L)
for (i in 1:L){
  yu[i] = min(1-eps, 1-eps/delta*xu[i])
}
## Get the data of lower boundary
delta = n / p
eps = s / p
tppmax = powermax(delta, eps)
xx = seq(0, tppmax, length.out = L)
yy = rep(0,L)
for (i in 1:L){
  yy[i] = fdrlasso(xx[i], delta, eps)
}

## Plot the simulation data for various levels of noise.
par(mar=c(4.5,5.5,1,2),xpd=FALSE)
plot(xx,yy,type='l',lwd=3,ylim=c(0,1),xlim=c(0,1),cex.lab=2.5,cex.axis=1.5,
     xaxs = "i", yaxs = "i",ylab='FDP',xlab='TPP')
lines(xu,yu,lwd=3)
colors=rev(brewer.pal(n = 8, name = "Dark2"))
for(i in 1:N){
  lines(tpp[[i]],fdp[[i]],lty=5,lwd=4.5,col=colors[i])
}
