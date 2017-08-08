# plot integrated depth with random tukey marginals

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

density.value <- c(20,50,80)

density.1 <- density(smooth.gp.data[,density.value[1]])
density.2 <- density(smooth.gp.data[,density.value[2]])
density.3 <- density(smooth.gp.data[,density.value[3]])
density.upscale <- 3

function.of.choice <- 12

if(print == T){pdf("MT_COVER_integrated_depth_plot.pdf",width = 7.5, height = 3)}

par(mfrow = c(1,1), mar = c(0,0,0,0) + 0.1, cex = 1.25)
rainbow.plot(smooth.gp.data, x.star, col = "grey", ylab="", xlab="", main="", xaxt="n", yaxt="n")


des.pol <- 30

cut.off.1 <- smooth.gp.data[function.of.choice,density.value[1]]
x.polygon.1 <- density.1$x[density.1$x > cut.off.1]
y.polygon.1 <- (density.1$y*density.upscale+x.star[density.value[1]])[1:length(x.polygon.1)]
polygon.1.1 <- cbind((rep(x.star[density.value[1]],length(x.polygon.1))), rev(x.polygon.1))
polygon.1.2 <- cbind(rev(y.polygon.1),(x.polygon.1))
polygon.1 <- rbind(polygon.1.1,polygon.1.2)

polygon(polygon.1[,1],polygon.1[,2], density = des.pol, border = NULL)

cut.off.2 <- smooth.gp.data[function.of.choice,density.value[2]]
x.polygon.2 <- density.2$x[density.2$x > cut.off.2]
y.polygon.2 <- ((density.2$y*density.upscale+x.star[density.value[2]])[max(density.2$x)-length(x.polygon.2):((max(density.2$x)))])[5:260]
polygon.2.1 <- cbind((rep(x.star[density.value[2]],length(x.polygon.2))),rev(x.polygon.2))
polygon.2.2 <- cbind((y.polygon.2),(x.polygon.2))
polygon.2 <- rbind(polygon.2.1,polygon.2.2)

polygon(polygon.2[,1],polygon.2[,2], density = des.pol, border = NULL)

cut.off.3 <- smooth.gp.data[function.of.choice,density.value[3]]
x.polygon.3 <- density.3$x[density.3$x < cut.off.3]
y.polygon.3 <- (density.3$y*density.upscale+x.star[density.value[3]])[1:length(x.polygon.3)]
polygon.3.1 <- cbind((rep(x.star[density.value[3]],length(x.polygon.3))),rev(x.polygon.3))
polygon.3.2 <- cbind((y.polygon.3),(x.polygon.3))
polygon.3 <- rbind(polygon.3.1,polygon.3.2)

polygon(polygon.3[,1],polygon.3[,2], density = des.pol, border = NULL)

lines(density.1$y*density.upscale+x.star[density.value[1]], density.1$x, col = "red", lwd = 2)
abline(v = x.star[density.value[1]], col = "red", lty = 2)
lines(density.2$y*density.upscale+x.star[density.value[2]], density.2$x, col = "red", lwd = 2)
abline(v = x.star[density.value[2]], col = "red", lty = 2)
lines(density.3$y*density.upscale+x.star[density.value[3]], density.3$x, col = "red", lwd = 2)
abline(v = x.star[density.value[3]], col = "red", lty = 2)
lines(x.star, smooth.gp.data[function.of.choice,], col = "black", lwd = 2)


# text(density.value[1]/10+0.2, -3.5, label = expression("t"[1]), col = "red")
# text(density.value[2]/10+0.2, -3.5, label = expression("t"[2]), col = "red")
# text(density.value[3]/10+0.2, -3.5, label = expression("t"[3]), col = "red")
# 
# arrows(0.45,2,0.95,0.35, length = 0.15)
# text(0.35,2.3, label = expression("x"))
# 
# arrows(3.0,3.2,2.05,2.28, length = 0.15)
# text(3,3.5, label = expression("D"["T1"]*"(x("*"t"[1]*"),"*"P"["t"[1]*",n"]*")"))


depth.x <- rep(NA, n.values)

for(i in 1:n.values){
  depth.x[i] <- depth(smooth.gp.data[function.of.choice,i],smooth.gp.data[,i])
}

# plot(x.star, depth.x, type = "l", xlab = "t", ylab = expression("D"["T1"]*"(x(t),P"["t,n"]*")"), ylim = c(0,0.55))
# polygon.depth.x <- polygon.fct(x.star, rbind(rep(0, n.values),depth.x))
# polygon(polygon.depth.x, border = NULL, density = des.pol)
# arrows(3.1,0.48,4.2,0.2, length = 0.15)
# text(3.1,0.51, label = expression("D"["ID"]*"(x,P"["n"]*",D"["T1"]*")"))
# 
# abline(v = x.star[density.value[1]], col = "red", lty = 2)
# abline(v = x.star[density.value[2]], col = "red", lty = 2)
# abline(v = x.star[density.value[3]], col = "red", lty = 2)

# text(density.value[1]/10+0.2, -3.5, label = expression("t"[1]), col = "red")
# text(density.value[2]/10+0.2, -3.5, label = expression("t"[2]), col = "red")
# text(density.value[3]/10+0.2, -3.5, label = expression("t"[3]), col = "red")

if(print == T){dev.off()}
