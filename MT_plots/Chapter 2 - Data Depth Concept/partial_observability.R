# plot partial_observability

print <- T






n.draws <- 40 # Number of stochastic process drawn
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

smooth.gp.data <- smooth.gp.data[c(24,22),]



t.max <- c(10,14)

set.seed(19)
T1 <- sample(1:100,t.max[1])
set.seed(3)
T2 <- sample(1:100,t.max[2])

T1 <- sort(T1)
T2 <- sort(T2)



if(print == T){pdf("partial_observability.pdf", width = 10, height = 4.7)}

points.to.plot <- 1:60
cex <- 1.8

T1 <- T1[T1 < max(points.to.plot)]
T2 <- T2[T2 < max(points.to.plot)]


T2 <- c(7,T2, 48)
T2 <- T2[-8]

T1 <- T1[-3]
T1 <- c(T1,20)


par(mfrow = c(1,1), mar = c(4.2, 4, 1, 2) + 0.1, cex = 1.00)
rainbow.plot(smooth.gp.data[,points.to.plot], x.star[points.to.plot], col = "red", ylab = "x(t)")
lines(x.star[points.to.plot],smooth.gp.data[1,points.to.plot], col = "blue")
points(T1/10-0.06, smooth.gp.data[1,T1], pch = "*", cex = cex, col = "blue")
points(T2/10-0.06, smooth.gp.data[2,T2], pch = "*", cex = cex, col = "red")
unique.T1.T2 <- unique(c(T1,T2))

x.white <- rep(unique.T1.T2[7]/10-0.06,2)
y.white <- c(-1.5,-1.9)

#for(i in 1:length(unique.T1.T2)){
  abline(v = unique.T1.T2[i]/10-0.06, lty = 2)
  lines(x.white,y.white, col = "white", lwd = 4)
}





T1.T2 <- c(7,12,28,48,54)


for(i in 1:length(T1.T2)){
  abline(v = T1.T2[i]/10-0.06, lty = 2, col = "black", lwd = 1.0)
}


arrows(0.33,0.43,0.88,-0.72, length = 0.15)
text(0.23,0.56, label = expression(tilde("x")[1]), col = "red") 

arrows(2.15,-1.4,2.02,-0.93, length = 0.15)
text(2.2,-1.65, label = expression("(t"[2]*","[3]*", x"[2]*"(t"[2]*","[3]*"))"), col = "blue") 


# 
# arrows(0.45,2,0.95,0.35, length = 0.15)
# text(0.35,2.3, label = expression("x"))
# 
# arrows(3.0,3.2,2.05,2.28, length = 0.15)
# text(3,3.5, label = expression("D"["T1"]*"(x("*"t"[1]*"),"*"P"["t"[1]*",n"]*")"))


if(print == T){dev.off()}



