# plot to explain difference in P4 between depths for multivariate and functional data




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

function.of.choice <- 8
expansion.factor <- c(0.4,0.6,1.5,2.3,3.5)
area.curve.expansion <- 45:80
length.area <- length(area.curve.expansion)
expansion.function <- rep(0, length.area)
seq.expansion.function <- seq(-2.8,2.8, length = length.area)
expansion.function <- dnorm(seq.expansion.function)
normalized.expansion.function <- expansion.function/max(expansion.function)

x_n_k_1 <- smooth.gp.data[function.of.choice,]
x_n_k_1[area.curve.expansion] <- x_n_k_1[area.curve.expansion]*(normalized.expansion.function*expansion.factor[1]+1)
x_n_k_2 <- x_n_k_1
x_n_k_2[area.curve.expansion] <- x_n_k_1[area.curve.expansion]*(normalized.expansion.function*expansion.factor[2]+1)
x_n_k_3 <- x_n_k_2
x_n_k_3[area.curve.expansion] <- x_n_k_1[area.curve.expansion]*(normalized.expansion.function*expansion.factor[3]+1)
x_n_k_4 <- x_n_k_3
x_n_k_4[area.curve.expansion] <- x_n_k_1[area.curve.expansion]*(normalized.expansion.function*expansion.factor[4]+1)
x_n_k_5 <- x_n_k_4
x_n_k_5[area.curve.expansion] <- x_n_k_1[area.curve.expansion]*(normalized.expansion.function*expansion.factor[5]+1)

expansion.addition <- c(1,2.5,4,7,12)
x2_n_k_1 <- smooth.gp.data[function.of.choice,]
x2_n_k_1 <- x2_n_k_1 +expansion.addition[1]
x2_n_k_2 <- smooth.gp.data[function.of.choice,]
x2_n_k_2 <- x2_n_k_2 +expansion.addition[2]
x2_n_k_3 <- smooth.gp.data[function.of.choice,]
x2_n_k_3 <- x2_n_k_3 +expansion.addition[3]
x2_n_k_4 <- smooth.gp.data[function.of.choice,]
x2_n_k_4 <- x2_n_k_4 +expansion.addition[4]
x2_n_k_5 <- smooth.gp.data[function.of.choice,]
x2_n_k_5 <- x2_n_k_5 +expansion.addition[5]

if(print == T){pdf("vanishing_at_infinity_property.pdf", width = 10, height = 4.7)}

par(mfrow = c(1,2), mar = c(4.2, 4, 1, 2) + 0.1, cex = 1)

rainbow.plot(smooth.gp.data, x.star, col = "grey", ylab = "x(t)", ylim = c(min(smooth.gp.data),max(x_n_k_5)*1.1))
lines(x.star, x2_n_k_1, col = "black", lwd = 1, lty = 2)
lines(x.star, x2_n_k_2, col = "black", lwd = 1, lty = 2)
lines(x.star, x2_n_k_3, col = "black", lwd = 1, lty = 2)
lines(x.star, x2_n_k_4, col = "black", lwd = 1, lty = 2)
lines(x.star, x2_n_k_5, col = "black", lwd = 1, lty = 2)

arrows(2.8,9.5,4.6,8.3, length = 0.15)
text(1.7,9.8, label = expression("x"["m,1"]))

rainbow.plot(smooth.gp.data, x.star, col = "grey", ylab = "x(t)", ylim = c(min(smooth.gp.data),max(x_n_k_5)*1.1))
lines(x.star, smooth.gp.data[function.of.choice,], col = "black", lwd = 1, lty = 2)
lines(x.star, x_n_k_1, col = "black", lwd = 1, lty = 2)
lines(x.star, x_n_k_2, col = "black", lwd = 1, lty = 2)
lines(x.star, x_n_k_3, col = "black", lwd = 1, lty = 2)
lines(x.star, x_n_k_4, col = "black", lwd = 1, lty = 2)
lines(x.star, x_n_k_5, col = "black", lwd = 1, lty = 2)

arrows(3.15,8,5.2,7, length = 0.15)
text(2.1,8.5, label = expression("x"["m,2"]))

#lines(x.star, smooth.gp.data[function.of.choice.2,], col = "blue", lwd = 2)






# arrows(3.15,3.3,2.25,2.45, length = 0.15)
# text(3.52,3.86, label = expression("f"["x(t"[1]*")"]))
# 
# arrows(7.0,-2,7.7,0.3, length = 0.15)
# text(6.7,-2.5, label = expression("x(t"[3]*")"))



if(print == T){dev.off()}



