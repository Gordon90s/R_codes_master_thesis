# ===========================================================
# -------------  DATA SIMULATION ----------------------------
# ===========================================================

# WARNING: Please load packages from R Masterthesis MAIN

# Initialize simulation parameters

n.draws <- 1000 # Number of stochastic process drawn
n.values <- 100 # Number of values for each stochastic process
# WARNING: Don't forget te resimulate data after chaning the simulation parameters

# save all simulated data, arra[,,i] corresponds to data from
# Gaussian
# Gaussian with varying amplitudes
# Brown-Resnick
# extremal t
data.array <- array(NA, dim = c(n.draws, n.values, 4))


# =============================================================================
# SMOOTH GAUSSIAN PROCESS -----------------------------------------------------

# found here
# http://stats.stackexchange.com/questions/30652/how-to-simulate-functional-data

set.seed(1)
calcSigma <- function(X1,X2,l=1){
  Sigma <- matrix(rep(0, length(X1)*length(X2)), nrow = length(X1))
  for(i in 1:nrow(Sigma)){
  for (j in 1:ncol(Sigma)) Sigma[i,j]<-exp(-1/2*(abs(X1[i]-X2[j])/l)^2)
  }
  return(Sigma)
  }
# The standard deviation of the noise
x.star <- seq(-5,5, len=n.values)
nval <- 5
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

# rainbow.plot(smooth.gp.data, x.star, col = "black")
# rainbow.plot(smooth.gp.data, x.star)
# lines(x.star,f.bar.star,col="red",lwd=2)

data.array[,,1] <- smooth.gp.data


# =============================================================================
# SMOOTH GAUSSIAN PROCESS WITH AMPLITUDE---------------------------------------

n.draws <- 1 # Number of stochastic process drawn
n.values <- 200 # Number of values for each stochastic process
# WARNING: Don't forget te resimulate data after chaning the simulation parameters


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
# rainbow.plot(smooth.gp.data, x.star)

anker <- smooth.gp.data

set.seed(1)
n.draws <- 1000 # Number of stochastic process drawn

# WARNING: Don't forget te resimulate data after chaning the simulation parameters



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
# rainbow.plot(smooth.gp.data, x.star)

n.values <- 100
n.values.200.seq <- 1:100*2 
data.array[,,2] <- smooth.gp.data[,n.values.200.seq] # save only every second point to only have n.values = 100


BR.data <- readRDS("./MT_functions/Comparaison_data_simulations/BR_array_1000_100_100.rds")
data.array[,,3] <- BR.data[,,1]

ET.data <- readRDS("./MT_functions/Comparaison_data_simulations/extremal_t_array.rds_1000_100_1000.rds")
data.array[,,4] <- ET.data[,,1]

saveRDS(data.array, "./MT_functions/Comparaison_data_simulations/main_data_comparaison.rds")



