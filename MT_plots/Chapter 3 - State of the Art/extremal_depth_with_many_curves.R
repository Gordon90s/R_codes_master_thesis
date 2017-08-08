# Extremal depth simulations



print <- F


if(print == T){pdf("extremal_depth_many_curves.pdf", height = 5, width = 9)}



n.draws <- 400 # Number of stochastic process drawn
n.values <- 400 # Number of values for each stochastic process
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



x <- y <- smooth.gp.data



n1.x <- dim(x)[1]
n2.x <- dim(x)[2]
n1.y <- dim(y)[1]
n2.y <- dim(y)[2]


x.depth <- MDF(x, ED = T)

order.depth <- order(x.depth$function.depth.x)

col <- viridis_pal(alpha = 1, option = "A", end = 0.95)(dim(x)[1])

par(mfrow = c(1,2))



ED = 1
if(ED == T){

    r <- 0:n1.y/n1.y
    r.n <- length(r)
    
    D_x <- matrix(0, nrow = n1.x, ncol = n2.x)
  
   for(i in 1:n1.x){
     for(j in 1:n2.x){
      D_x[i,j] <- 1 - abs(sum(y[,j] < x[i,j]) - sum(y[,j] > x[i,j]))/n1.y
     }
    }
  
    Phi.x.r <- matrix(0, nrow = n1.x, ncol = r.n)
    for(i in 1:n1.x){
      for(j in 1:r.n){
        Phi.x.r[i,j] <- sum(D_x[i,] <= r[j])/n2.x
      }
    }
  
 # sort phi.x.r
    
    new.Phi.x.r <- cbind(Phi.x.r,(1:nrow(Phi.x.r))) # add index as last column
    
    # sort rows according to algorithm (sort first rows by values in first column, then values in second column, etc)
    new.Phi.x.r.sorted <- new.Phi.x.r[do.call(order, as.data.frame(new.Phi.x.r)),] 
    # retrieve index vector in last column and reverse it as sorting as not descending
    order.depth.x <- rev(new.Phi.x.r.sorted[,ncol(new.Phi.x.r.sorted)]) #function numbers are sorted from most to extrem to least extreme
    # for example order.depth.x = c(4,2,1,3) means x_4 is most extreme, x_3 is deepest.
    
    function.depth.x <- sort(order.depth.x, index.return = T)$ix # assigning each x_i = depth rank

} # end extremal depth  
  


# 
# for(i in 1:dim(Phi.x.r)[1]){
#   if(i == 1){plot(Phi.x.r[i,], type = "l")}else{lines(Phi.x.r[i,])}
# }






# fancy test Phi.x.r



col <- viridis_pal(alpha = 1, option = "A", end = 0.95)(dim(x)[1])

rainbow.plot(Phi.x.r[order.depth,], normalize(c(0,x.star)), col = col, ylab = expression(Phi["x"["i"]*",n"]*"(q)"), xlab = "q")





rainbow.plot(x[order.depth,],x.star, col = col, ylab = expression("x"["i"]*"(t)"))
#lines(1:6, c(2,2.1,1.9,1.3,0.8,-0.4), lwd = 4)


if(print == T){dev.off()}





