# rainbowplot for functional data to illustrate depth

print <- T
resimulate.data <- T



if(resimulate.data == T){

n.draws <- 50 # Number of stochastic process drawn
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
  for (j in 1:ncol(Sigma)) Sigma[i,j]<-1*exp(-1/2*(abs(X1[i]-X2[j])/l)^2)
  }
  return(Sigma)
  }
# The standard deviation of the noise
x.star <- seq(0,10, len=n.values)
nval <- 4
#f <- data.frame(x=seq(-0,10,l=nval), y=rnorm(nval,0,10))
f <- data.frame(x=c(0,4,7,10), y=c(-1,0,3,0))
sigma.n <- 0.4
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



if(print == T){pdf("gaussian_process_explained_sigma_T_not_zero.pdf",width = 10, height = 5.5)}

par(mfrow = c(1,1), mar = c(4.2, 4, 1, 2) + 0.1, cex = 1.25)


sorted.smooth.gp.data <- main.depth.function(smooth.gp.data, nproj = 10000, graph = F)$sorted.data
}

# col <- (palette(gray(seq(0.2,1,len = nrow(sorted.smooth.gp.data))))) # 50 shades of grey

# col <- c(viridis_pal(alpha = 1, option = "D", begin = 0.1, end = 0.5, direction = 1)(nrow(sorted.smooth.gp.data)*(3/4)),viridis_pal(alpha = 1, option = "C", begin = 0.0, end = 0.38, direction = 1)(nrow(sorted.smooth.gp.data)*(1/4)))

# col <- hsv(h = normalize(1:nrow(sorted.smooth.gp.data)),
           # s = normalize(1:nrow(sorted.smooth.gp.data)),
           # v = normalize(1:nrow(sorted.smooth.gp.data)), alpha = 1)

# col <- c(viridis_pal(alpha = 1, option = "B", begin = 0, end = 0.5)(nrow(sorted.smooth.gp.data)/2),viridis_pal(alpha = 0.6, option = "B", begin = 0.5, end = 1)(nrow(sorted.smooth.gp.data)/2))

# col <- c(viridis_pal(alpha = 1, option = "D", begin = 0, end = 0.5)(nrow(sorted.smooth.gp.data)/2),viridis_pal(alpha = 1, option = "D", begin = 0.5, end = 1)(nrow(sorted.smooth.gp.data)/2))

col <- viridis_pal(alpha = 1, option = "D", end = 1)(nrow(sorted.smooth.gp.data))

col <- col[sample(n.row, replace = T)]

      i <- 1
      plot(x.star, sorted.smooth.gp.data[i,], type = "l", ylim = c(min(sorted.smooth.gp.data), max(sorted.smooth.gp.data)), xlab = "t", ylab = "x(t)", col = "white")
      # plot with different greys, deeper -> white, less deep -> black
      for(i in 1:nrow(sorted.smooth.gp.data)){
      lines(x.star, sorted.smooth.gp.data[i,], col = "darkgrey")
      }
      lines(x.star, f.bar.star, col = "red", lwd = 2)
      points(f, pch = "*", cex = 2)
            arrows(3.8,2.5,3.97,0.4, length = 0.15)
      text(3.8,2.75,label = expression("y"[2]))
      arrows(7,-0.3,5.75,0.75, length = 0.15, col = "red")
      text(7.1,-0.50,label = expression("mean"), col = "red")

if(print == T){dev.off()}



