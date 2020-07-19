########################################################################################
## This R code is used for plotting Figure 2 and 3 in paper "The Complete Lasso Diagram"
##--------------------------------------------------------------------------------------
## Copyright @ Hua Wang, Yachong Yang, Zhiqi Bu, and Weijie Su, 2020
##--------------------------------------------------------------------------------------
########################################################################################
library(pracma)
#pdf("figure_dt_combine.pdf", width=20, height=5)
#par(mfrow=c(1, 3))
color = 'indianred2'
color_inner = 'skyblue2'

source("q_lower.R")

## Case 1
## Get the upper boundary
delta = 1
eps = 0.3
L=1000
tppmax = powermax(delta, eps)
xu = seq(0, tppmax, length.out = L)
yu = rep(0,L)
for (i in 1:L){
  yu[i] = min(1-eps, 1-eps/delta*xu[i])
}


## Get the data of lower boundary
xx = seq(0, tppmax, length.out = L)
yy = rep(0,L)
for (i in 1:L){
  yy[i] = fdrlasso(xx[i], delta, eps)
}

par(mar=c(5,6,1,2),xpd=FALSE)
plot(0,0,xlim = c(0, 1), ylim=c(0, 1), xaxs='i', yaxs='i',xlab='TPP',ylab='FDP', cex.lab=2.5,cex.axis=1.5, type="n")
polygon(c(0,1,1,0),c(0,0,1,1),col=color)
polygon(c(xx,rev(xu)),c(yy,rev(yu)),col=color_inner)
text(0.5, 0.85, 'Unachievable', cex=2)
text(0.5,0.45,'Achievable',cex=2)
text(0.8, 0.08, 'Unachievable', cex=2)


##################################################
## Case 2
## Get the upper boundary
delta = 0.7
eps = 0.3
L=1000
tppmax = min(powermax(delta, eps), 1)
xu = seq(0, tppmax, length.out = L)
yu = rep(0,L)
for (i in 1:L){
  yu[i] = min(1-eps, 1-eps/delta*xu[i])
}


## Get the data of lower boundary
xx = seq(0, tppmax, length.out = L)
yy = rep(0,L)
for (i in 1:L){
  yy[i] = fdrlasso(xx[i], delta, eps)
}

par(mar=c(5,6,1,2),xpd=FALSE)
plot(0,0,xlim = c(0, 1), ylim=c(0, 1), xaxs='i', yaxs='i',xlab='TPP',ylab='FDP', cex.lab=2.5,cex.axis=1.5, type="n")
polygon(c(0,1,1,0),c(0,0,1,1),col=color)
polygon(c(xx,rev(xu)),c(yy,rev(yu)),col=color_inner)
text(0.5, 0.85, 'Unachievable', cex=2)
text(0.5,0.45,'Achievable',cex=2)
text(0.8,0.15,'Unachievable',cex=2)

##################################################

## Case 3
## Get the upper boundary
delta = 0.5
eps = 0.3
L=1000
tppmax = powermax(delta, eps)
xu = seq(0, tppmax, length.out = L)
yu = rep(0,L)
for (i in 1:L){
  yu[i] = min(1-eps, 1-eps/delta*xu[i])
}


## Get the data of lower boundary
xx = seq(0, tppmax, length.out = L)
yy = rep(0,L)
for (i in 1:L){
  yy[i] = fdrlasso(xx[i], delta, eps)
}

par(mar=c(5,6,1,2),xpd=FALSE)
plot(0,0,xlim = c(0, 1), ylim=c(0, 1), xaxs='i', yaxs='i',xlab='TPP',ylab='FDP', cex.lab=2.5,cex.axis=1.5, type="n")
polygon(c(0,1,1,0),c(0,0,1,1),col=color)
polygon(c(xx,rev(xu)),c(yy,rev(yu)),col=color_inner)
abline(v =tppmax, lwd=2,lty=2)
text(0.5, 0.85, 'Unachievable', cex=2)
text(0.4,0.48,'Achievable',cex=2)
text(0.7,0.2,'Unachievable',cex=2)
#dev.off()
