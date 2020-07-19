##################################################################################
## This R code is used for plotting Figure 1 in paper "The Complete Lasso Diagram"
##--------------------------------------------------------------------------------
## Copyright @ Hua Wang, Yachong Yang, Zhiqi Bu, and Weijie Su, 2020
##--------------------------------------------------------------------------------
##################################################################################

library(pracma)
#pdf("figure_su2017_combine.pdf", width=20, height=5)
#par(mfrow=c(1, 3))
color = 'indianred2'

source("q_lower.R")

len = 100
n = len
delta = 1
eps = 0.3
tppmax = powermax(delta, eps)
tpp = seq(0, tppmax, length.out = len)
fdp = rep(0, len)
for (i in 1:len){
  fdp[i] = fdrlasso(tpp[i], delta, eps)
}



par(mar=c(5,6,1,2),xpd=FALSE)
plot(0,0,xlim = c(0, 1), ylim=c(0, 1), xaxs='i', yaxs='i',xlab='TPP',ylab='FDP', cex.lab=2.5,cex.axis=1.5, type="n")
polygon(c(0, tpp, 1, 1), c(0, fdp, 2*fdp[n]-fdp[n-1], 0), col=color)
#axis(1,cex.axis=1.8)
#axis(2, las=1,cex.axis=1.8)
#box()
#title(xlab='TPP', ylab='FDP', cex.lab=2.2)
text(0.8, 0.08, 'Unachievable', cex=2)
text(0.4, 0.5, 'Possibly Achievable', cex=2)

###########################
len = 100
n = len
delta = 0.7
eps = 0.3
tppmax = powermax(delta, eps)
tpp = seq(0, tppmax, length.out = len)
fdp = rep(0, len)
for (i in 1:len){
  fdp[i] = fdrlasso(tpp[i], delta, eps)
}



par(mar=c(5,6,1,2),xpd=FALSE)
plot(0,0,xlim = c(0, 1), ylim=c(0, 1), xaxs='i', yaxs='i',xlab='TPP',ylab='FDP', cex.lab=2.5,cex.axis=1.5, type="n")
polygon(c(0, tpp, 1, 1), c(0, fdp, 2*fdp[n]-fdp[n-1], 0), col=color)
#axis(1,cex.axis=1.8)
#axis(2, las=1,cex.axis=1.8)
#box()
#title(xlab='TPP', ylab='FDP', cex.lab=2.2)
text(0.8, 0.15, 'Unachievable', cex=2)
text(0.4, 0.5, 'Possibly Achievable', cex=2)


##################################################
len = 100
n = len
delta = 0.5
eps = 0.3
tppmax = powermax(delta, eps)
tpp = seq(0, tppmax, length.out = len)
fdp = rep(0, len)
for (i in 1:len){
  fdp[i] = fdrlasso(tpp[i], delta, eps)
}

par(mar=c(5,6,1,2),xpd=FALSE)
plot(0,0,xlim = c(0, 1), ylim=c(0, 1), xaxs='i', yaxs='i',xlab='TPP',ylab='FDP', cex.lab=2.5,cex.axis=1.5, type="n")
polygon(c(0, tpp, tpp[n], 1, 1), c(0, fdp, 1,  1, 0), col=color)
abline(v =tppmax, lwd=2,lty=2)
#axis(1,cex.axis=1.8)
#axis(2, las=1,cex.axis=1.8)
#box()
#title(xlab='TPP', ylab='FDP', cex.lab=2.2)
text(0.7, 0.2, 'Unachievable', cex=2)
text(0.4, 0.5, 'Possibly Achievable', cex=2)
dev.off()

