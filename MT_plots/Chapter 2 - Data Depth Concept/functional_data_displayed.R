# plot to explain functional data and its marginal distribution

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

density.value <- c(20,50,80)

density.1 <- density(smooth.gp.data[,density.value[1]])
density.2 <- density(smooth.gp.data[,density.value[2]])
density.3 <- density(smooth.gp.data[,density.value[3]])
density.upscale <- 3

function.of.choice <- 12

if(print == T){pdf("functional_data_displayed.pdf", width = 10, height = 5.5)}

par(mfrow = c(1,1), mar = c(4.2, 4, 1, 2) + 0.1, cex = 1.25)
rainbow.plot(smooth.gp.data, x.star, col = "grey", ylab = "x(t)")
lines(density.1$y*density.upscale+x.star[density.value[1]], density.1$x, col = "red")
abline(v = x.star[density.value[1]], col = "red", lty = 2)
lines(density.2$y*density.upscale+x.star[density.value[2]], density.2$x, col = "red")
abline(v = x.star[density.value[2]], col = "red", lty = 2)
lines(density.3$y*density.upscale+x.star[density.value[3]], density.3$x, col = "red")
abline(v = x.star[density.value[3]], col = "red", lty = 2)
points(rep(x.star[density.value[1]],n.draws),smooth.gp.data[,density.value[1]], col = "black", pch = "*")
points(rep(x.star[density.value[2]],n.draws),smooth.gp.data[,density.value[2]], col = "black", pch = "*")
points(rep(x.star[density.value[3]],n.draws),smooth.gp.data[,density.value[3]], col = "black", pch = "*")

text(density.value[1]/10+0.2, -3.5, label = expression("t"[1]), col = "red")
text(density.value[2]/10+0.2, -3.5, label = expression("t"[2]), col = "red")
text(density.value[3]/10+0.2, -3.5, label = expression("t"[3]), col = "red")


arrows(3.15,3.3,2.25,2.45, length = 0.15, col = "red")
text(3.45,3.5, label = expression("f"["t"[1]*",n"]), col = "red")

arrows(7.0,-2,7.7,0.3, length = 0.15)
text(6.7,-2.5, label = expression("x(t"[3]*")"))



if(print == T){dev.off()}



