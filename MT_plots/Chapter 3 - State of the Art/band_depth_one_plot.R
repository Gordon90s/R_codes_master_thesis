# plot band depth


print <- F

n.draws <- 15 # Number of stochastic process drawn
n.values <- 400 # Number of values for each stochastic process
# WARNING: Don't forget te resimulate data after chaning the simulation parameters

  
# =============================================================================
# SMOOTH GAUSSIAN PROCESS -----------------------------------------------------

# found here
# http://stats.stackexchange.com/questions/30652/how-to-simulate-functional-data

set.seed(6)
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

smooth.gp.data[1,] <- smooth.gp.data[1,]+1.8
smooth.gp.data[4,] <- smooth.gp.data[4,]-2
smooth.gp.data[12,] <- smooth.gp.data[12,] + 2


new.smooth.gp.data <- rbind(smooth.gp.data[1,],smooth.gp.data[4,],smooth.gp.data[6,],smooth.gp.data[12,])

smooth.gp.data <- new.smooth.gp.data

cex.title <- 1.35
width.fig <- 10
height.fig <- 5



if(print == T){pdf("band_depth.pdf",width = 9, height = 9)}

par(mfrow = c(1,1), cex = 1.6)
rainbow.plot(smooth.gp.data, x.star, col = "grey", ylab = "x(t)", ylim = c(-4,6.5))
polygon(polygon.fct(x.star, smooth.gp.data)$x,polygon.fct(x.star, smooth.gp.data)$y , col = "red", border = NA)
# lines(x.star, smooth.gp.data[3,], col = "black", lwd = 3)
lines(x.star, smooth.gp.data[1,], col = "blue", lwd = 4, lty = 1)
lines(x.star, smooth.gp.data[2,], col = "blue", lwd = 4, lty = 1)
lines(x.star, smooth.gp.data[4,], col = "blue", lwd = 4)


text(8.4,2, label = expression("B"[3]*"[x"[1]*",x"[2]*",x"[3]*"]"), cex = cex.title)
text(0.55,-1.7, label = expression("x"[3]), cex = 1)
text(0.55,4.1, label = expression("x"[1]), cex = 1)
text(0.55,1.3, label = expression("x"[2]), cex = 1)
# arrows(6.1,-2,4.3,-0.8, length = 0.15)
# text(6.4,-2.2, label = expression("x"))



# arrows(7.2,-2,6.7,0.2, length = 0.15)
# text(7.35,-2.3, label = expression("v"[1]))
# 
# arrows(1.42,2.8,1.9,1.299, length = 0.15)
# text(1.34,3.05, label = expression("x"))


if(print == T){dev.off()}