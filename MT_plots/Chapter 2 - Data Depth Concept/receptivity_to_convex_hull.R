# plot to explain difference in P4 between depths for multivariate and functional data

print <- F

n.draws <- 1 # Number of stochastic process drawn
n.values <- 200 # Number of values for each stochastic process
# WARNING: Don't forget te resimulate data after chaning the simulation parameters

  
# =============================================================================
# SMOOTH GAUSSIAN PROCESS -----------------------------------------------------

# found here
# http://stats.stackexchange.com/questions/30652/how-to-simulate-functional-data

set.seed(8)
calcSigma <- function(X1,X2,l=20){
  Sigma <- matrix(rep(0, length(X1)*length(X2)), nrow = length(X1))
  for(i in 1:nrow(Sigma)){
  for (j in 1:ncol(Sigma)) Sigma[i,j]<-exp(-1/2*(abs(X1[i]-X2[j])/l)^2)
  }
  return(Sigma)
  }
# The standard deviation of the noise
x.star <- 1:n.values
nval <- 5
f <- data.frame(x=seq(-10,n.values,l=nval), y=rnorm(nval,0,10))
sigma.n <- 1
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

# rainbow.plot(smooth.gp.data, x.star, col = "black")
rainbow.plot(smooth.gp.data, x.star)


anker <- smooth.gp.data

set.seed(1)
n.draws <- 800 # Number of stochastic process drawn

# WARNING: Don't forget te resimulate data after chaning the simulation parameters

  
# =============================================================================
# SMOOTH GAUSSIAN PROCESS -----------------------------------------------------

# found here
# http://stats.stackexchange.com/questions/30652/how-to-simulate-functional-data

values <- matrix (rep(0,length(x.star)*n.draws), ncol=n.draws)
for(i in 1:n.draws){
calcSigma <- function(X1,X2,l=20){
  Sigma <- matrix(rep(0, length(X1)*length(X2)), nrow = length(X1))
  for(i in 1:nrow(Sigma)){
  for (j in 1:ncol(Sigma)) Sigma[i,j]<-exp(-1/2*(abs(X1[i]-X2[j])/l)^2)
  }
  return(Sigma)
  }
# The standard deviation of the noise
x.star <- 1:n.values
nval <- 50

seq.nval <- seq(1,n.values, len = nval)
seq.nval1 <- seq(-10,n.values+50, len = nval)+rnorm(nval,0,1)
seq.nval2 <- c(rnorm(26,0,0.1),-abs(rnorm(24,-5,25)))


f <- data.frame(x=seq.nval1, y= seq.nval2)
sigma.n <- 1
# Recalculate the mean and covariance functions
k.xx <- calcSigma(f$x,f$x)
k.xxs <- calcSigma(f$x, x.star)
k.xsx <- calcSigma(x.star,f$x)
k.xsxs <- calcSigma(x.star,x.star)
f.bar.star <- k.xsx%*%solve(k.xx+sigma.n^2*diag(1,ncol(k.xx)))%*%f$y
cov.f.star <- k.xsxs-k.xsx%*%solve(k.xx+sigma.n^2*diag(1,ncol(k.xx)))%*%k.xxs

values[,i] <- anker + mvrnorm(1,f.bar.star,cov.f.star)
}

smooth.gp.data <- t(values) # save data, need to transpose for having good matrix

# rainbow.plot(smooth.gp.data, x.star, col = "black")
rainbow.plot(smooth.gp.data, x.star)


data <- smooth.gp.data
smooth.gp.data.new <- smooth.gp.data



# depth.new.ED <- main.depth.function(smooth.gp.data.new, depth.type = 6, graph = F)$depth
depth.new.ED <- MDF(smooth.gp.data.new, x.star = x.star, ED = T)$function.depth.x # apply extremal depth
# depth.new.ED <- main.depth.function(smooth.gp.data.new, depth.type = 1, nproj = 1000, graph = F)$depth
depth.new.alpha.volume <- MDF(smooth.gp.data.new, n.projections = n.proj, x.star = x.star)$function.depth.x

order.depth.new.ED <- order(depth.new.ED)
order.depth.new.alpha.volume <- order(depth.new.alpha.volume)

col <- viridis_pal(alpha = 1, option = "A", end = 1)(n.draws)

smooth.gp.data.new.RT <- smooth.gp.data.new.alpha.volume <- matrix(0, nrow = n.draws, ncol = n.values) 
for(i in 1:n.draws){
  smooth.gp.data.new.RT[i,] <- smooth.gp.data.new[order.depth.new.ED[i],]
  smooth.gp.data.new.alpha.volume[i,] <- smooth.gp.data.new[order.depth.new.alpha.volume[i],]
}




if(print == T){pdf("receptivity_to_convex_hull.pdf", width = 10, height = 4.7)}

par(mfrow = c(1,2), mar = c(4.2, 4, 1, 2) + 0.1, cex = 1.0)

rainbow.plot(smooth.gp.data.new.RT, seq(0,10, length = 200), col = col, ylab = "x(t)", ylim = range(rbind(smooth.gp.data.new,smooth.gp.data)))
rainbow.plot(smooth.gp.data.new.alpha.volume, seq(0,10, length = 200), col = col, ylab = "x(t)", ylim = range(rbind(smooth.gp.data.new,smooth.gp.data)))



#lines(x.star, smooth.gp.data[function.of.choice.2,], col = "blue", lwd = 2)

# arrows(3.15,3.3,2.25,2.45, length = 0.15)
# text(3.52,3.86, label = expression("f"["x(t"[1]*")"]))
# 
# arrows(7.0,-2,7.7,0.3, length = 0.15)
# text(6.7,-2.5, label = expression("x(t"[3]*")"))

if(print == T){dev.off()}



