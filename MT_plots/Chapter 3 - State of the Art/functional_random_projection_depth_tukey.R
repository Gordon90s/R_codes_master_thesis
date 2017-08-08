# plot functional random projection depth with Tukey depth marginals

print <- F






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

if(print == T){pdf("functional_random_tukey.pdf",width = 10, height = 5.5)}

par(mfrow = c(1,1), mar = c(4.2, 4, 1, 2) + 0.1, cex = 1.25)
rainbow.plot(smooth.gp.data, x.star, col = "grey", ylab = "x(t)")
lines(brown.motion, col = "red", lty = 1)
lines(x.star, smooth.gp.data[function.to.choose,], col = "black", lwd = 2)

arrows(7.2,-2,6.7,0.2, length = 0.15, col = "red")
text(7.37,-2.38, label = expression("v"[1]), col = "red")

arrows(8.3,4,8,2.1, length = 0.15, col = "red")
text(8.35,4.3, label = expression("v"[2]), col = "red")

arrows(1.42,2.8,1.9,1.299, length = 0.15)
text(1.30,3.1, label = expression("x"))

if(print == T){dev.off()}

# arrows(3.0,3.2,2.05,2.28, length = 0.15)
# text(3,3.5, label = expression("D"["T1"]*"(x("*"t"[1]*"),"*"P"["t"[1]*",n"]*")"))

int.x.v1 <- int.x.v2 <- rep(0,n.draws)
for(i in 1:n.draws){
int.x.v1[i] <- sum(brown.motion[1,]*smooth.gp.data[i,])/10
int.x.v2[i] <- sum(brown.motion[2,]*smooth.gp.data[i,])/10
}

des.pol <- 30
cex <- 1.5

rep.zero <- rep(0,n.draws)

density.x.v1 <- density(int.x.v1)
length.density.x.v1 <- length(density.x.v1$x)
rep.zero.density <- rep(0, length.density.x.v1)
x.range.v1 <- c(1:230)
polygon.x.v1 <- polygon.fct(density.x.v1$x[x.range.v1], rbind(density.x.v1$y[x.range.v1],rep.zero.density[x.range.v1]))

density.x.v2 <- density(int.x.v2)
length.density.x.v2 <- length(density.x.v2$x)
rep.zero.density <- rep(0, length.density.x.v2)
x.range.v2 <- c(1:139)
polygon.x.v2 <- polygon.fct(density.x.v2$x[x.range.v2], rbind(density.x.v2$y[x.range.v2],rep.zero.density[x.range.v2]))

graph.dim <- 1.1

if(print == T){pdf("functional_random_tukey_marginals.pdf",width = 10, height = 4.7)}

par(mfrow = c(1,2), mar = c(4.2, 4, 1, 2) + 0.1, cex = 0.95)

plot(int.x.v1, rep.zero, pch = "*", ylab = "density(y)", xlab = "y", ylim = c(0,0.1), xlim = range(density.x.v1$x))
polygon(x = polygon.x.v1$x, y = polygon.x.v1$y, density = des.pol, border = NULL)
abline(v = int.x.v1[function.to.choose], lwd = 2, lty = 1, col = "white")

lines(density(int.x.v1), col = "red", lwd = 2)
abline(v = int.x.v1[function.to.choose], lwd = 1, lty = 2, col = "red")

points(int.x.v1, rep.zero, pch = "*", cex = cex)
points(int.x.v1[function.to.choose], 0, pch = "*", cex = 2.5, col = "red")

arrows(5,0.009,1.4,0.0025, length = 0.15, col = "red")
text(6.4,0.0105, label = expression("y"["v"[1]]), col = "red")
arrows(8.5,0.092,-1.5,0.052, length = 0.15)
text(8.50,0.096, label = expression("min{F"["n,v"[1]]*"(y"["v"[1]]*"),1-F"["n,v"[1]]*"(y"["v"[1]]*")}"))

arrows(7.5+ shift2,0.040,4.5 + shift2,0.034, length = 0.15, col = "red")
text(8.72+ shift2,0.0418, label = expression("f"["v"[1]*",n"]), col = "red")

shift <- 1.4
shift2 <- 3.8

plot(int.x.v2, rep.zero, pch = "*", ylab = "density(y)", xlab = "y", ylim = c(0,0.1), xlim = range(density.x.v2$x))
polygon(x = polygon.x.v2$x, y = polygon.x.v2$y, density = des.pol, border = NULL)
abline(v = int.x.v2[function.to.choose], lwd = 1.5, lty = 1, col = "white")
lines(density(int.x.v2), col = "red", lwd = 2)
abline(v = int.x.v2[function.to.choose], lwd = 1, lty = 2, col = "red")
points(int.x.v2, rep.zero, pch = "*", cex = cex)
points(int.x.v2[function.to.choose], 0, pch = "*", cex = 2.5, col = "red")
arrows(5+shift,0.009,1.22+shift,0.002, length = 0.15, col = "red")
text(6.7+shift,0.0109, label = expression("y"["v"[2]]), col = "red")
arrows(16.5+ shift2,0.040,10.5 + shift2,0.034, length = 0.15, col = "red")
text(18.72+ shift2,0.042, label = expression("f"["v"[2]*",n"]), col = "red")

arrows(16.5,0.072,0.5,0.032, length = 0.15)
text(17.590,0.077, label = expression("min{F"["n,v"[2]]*"(y"["v"[2]]*"),1-F"["n,v"[2]]*"(y"["v"[2]]*")}"))



if(print == T){dev.off()}



