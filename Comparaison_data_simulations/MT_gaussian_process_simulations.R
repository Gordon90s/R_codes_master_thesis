# set appropriate directory

main.directory <- "/Users/gordonschucker/Dropbox/MASTER/Masterarbeit & Paper/R & LaTeX"
setwd(main.directory)

  calcSigma <- function(X1,X2,l=1){
  Sigma <- matrix(rep(0, length(X1)*length(X2)), nrow = length(X1))
  for(i in 1:nrow(Sigma)){
  for (j in 1:ncol(Sigma)) Sigma[i,j]<-exp(-1/2*(abs(X1[i]-X2[j])/l)^2)
  }
  return(Sigma)
  }

#--------------------------------------------------------

ptm.gauss <- proc.time()
# Monte Carlo simulations of Gaussian processes


set.seed(1)

n.draws <- 1000 # Number of stochastic process drawn
n.values <- 250 # Number of values for each stochastic process
mc <- 1 # Number of Monte Carlo simulations

x.star <- seq(-5,5, len=n.values)

gaussian.process.array <- array(NA, dim = c(n.draws, n.values, mc))
gaussian.process.data <- matrix(rep(0,n.values*n.draws),nrow=n.draws)                  
                  
for(k in 1:mc){

# The standard deviation of the noise

nval <- 3
f <- data.frame(x=seq(-5,5,l=nval), y=rnorm(nval,0,10))
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

# rainbow.plot(smooth.gp.data, x.star)

gaussian.process.data <- smooth.gp.data  
gaussian.process.array[,,k] <- gaussian.process.data

}                  

saveRDS(gaussian.process.array, "./MT_functions/Comparaison_data_simulations/gaussian_process_array.rds_1000_250_1000.rds")

time.gauss <- proc.time() - ptm.gauss
time.gauss

# ========================================================================
# ========================================================================
# ========================================================================


set.seed(2)

n.draws <- 1000 # Number of stochastic process drawn
n.values <- 100 # Number of values for each stochastic process
mc <- 1 # Number of Monte Carlo simulations

x.star <- seq(-5,5, len=n.values)

gaussian.process.array <- array(NA, dim = c(n.draws, n.values, mc))
gaussian.process.data <- matrix(rep(0,n.values*n.draws),nrow=n.draws)                  
                  
for(k in 1:mc){

# The standard deviation of the noise

nval <- 3
f <- data.frame(x=seq(-5,5,l=nval), y=rnorm(nval,0,10))
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

# rainbow.plot(smooth.gp.data, x.star)

gaussian.process.data <- smooth.gp.data  
gaussian.process.array[,,k] <- gaussian.process.data

}

saveRDS(gaussian.process.array, "./MT_functions/Comparaison_data_simulations/gaussian_process_array.rds_1000_100_1000.rds")

time.gauss <- proc.time() - ptm.gauss
time.gauss


# ========================================================================
# ========================================================================
# ========================================================================


set.seed(3)

n.draws <- 1000 # Number of stochastic process drawn
n.values <- 50 # Number of values for each stochastic process
mc <- 1 # Number of Monte Carlo simulations

x.star <- seq(-5,5, len=n.values)

gaussian.process.array <- array(NA, dim = c(n.draws, n.values, mc))
gaussian.process.data <- matrix(rep(0,n.values*n.draws),nrow=n.draws)                  
                  
for(k in 1:mc){

# The standard deviation of the noise

nval <- 3
f <- data.frame(x=seq(-5,5,l=nval), y=rnorm(nval,0,10))
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

# rainbow.plot(smooth.gp.data, x.star)

gaussian.process.data <- smooth.gp.data  
gaussian.process.array[,,k] <- gaussian.process.data

}                 

saveRDS(gaussian.process.array, "./MT_functions/Comparaison_data_simulations/gaussian_process_array.rds_1000_50_1000.rds")

time.gauss <- proc.time() - ptm.gauss
time.gauss


# for testing purposes
# n.draws <- 2# Number of stochastic process drawn
# n.values <- 5 # Number of values for each stochastic process
# mc <- 2 # Number of Monte Carlo simulatio