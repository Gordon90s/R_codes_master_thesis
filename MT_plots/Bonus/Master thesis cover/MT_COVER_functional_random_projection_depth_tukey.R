# plot functional random projection depth with Tukey depth marginals

print <- T






n.draws <- 20 # Number of stochastic process drawn
n.values <- 100 # Number of values for each stochastic process
# WARNING: Don't forget te resimulate data after chaning the simulation parameters

  
# =============================================================================
# SMOOTH GAUSSIAN PROCESS -----------------------------------------------------

# found here
# http://stats.stackexchange.com/questions/30652/how-to-simulate-functional-data

set.seed(2)
calcSigma <- function(X1,X2,l=1){
  Sigma <- matrix(rep(0, length(X1)*length(X2)), nrow = length(X1))
  for(i in 1:nrow(Sigma)){
  for (j in 1:ncol(Sigma)) Sigma[i,j]<-exp(-1/2*(abs(X1[i]-X2[j])/l)^2)
  }
  return(Sigma)
  }
# The standard deviation of the noise
x.star <- seq(0,10, len=n.values)
nval <- 4
f <- data.frame(x=seq(-0,10,l=nval), y=rnorm(nval,0,10))
sigma.n <- 2
# Recalculate the mean and covariance functions
k.xx <- calcSigma(f$x,f$x)
k.xxs <- calcSigma(f$x, x.star)
k.xsx <- calcSigma(x.star,f$x)
k.xsxs <- calcSigma(x.star,x.star)
f.bar.star <- k.xsx%*%solve(k.xx+sigma.n^2*diag(1,ncol(k.xx)))%*%f$y
cov.f.star <- k.xsxs-k.xsx%*%solve(k.xx+sigma.n^2*diag(1,ncol(k.xx)))%*%k.xxs
values <- matrix (rep(0,length(x.star)*n.draws), ncol=n.draws)
for (i in 1:n.draws)  values[,i] <- mvrnorm(1,f.bar.star,cov.f.star)
smooth.gp.data <- t(values) # save data, need to transpose for having good matrix

n.draws.BM <- 2
function.to.choose <- 8


set.seed(4)
brown.motion <- rproc2fdata(n.draws.BM,t=x.star,sigma="brownian")

if(print == T){pdf("MT_COVER_functional_random_tukey.pdf",width = 7.5, height = 3)}

par(mfrow = c(1,1), mar = c(0,0,0,0) + 0.1, cex = 1.25)
rainbow.plot(smooth.gp.data, x.star, col = "grey", ylab="", xlab="", main="", xaxt="n", yaxt="n")
lines(brown.motion, col = "red", lty = 1, lwd = 2)
lines(x.star, smooth.gp.data[function.to.choose,], col = "black", lwd = 2)

# arrows(7.2,-2,6.7,0.2, length = 0.15, col = "red")
# text(7.37,-2.38, label = expression("v"[1]), col = "red")
# 
# arrows(8.3,4,8,2.1, length = 0.15, col = "red")
# text(8.35,4.3, label = expression("v"[2]), col = "red")
# 
# arrows(1.42,2.8,1.9,1.299, length = 0.15)
# text(1.30,3.1, label = expression("x"))

if(print == T){dev.off()}



