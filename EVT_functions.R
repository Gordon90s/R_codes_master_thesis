# ===========================================================
# -----------------------------------------------------------
# --------------------- EVT FUNCTIONS -----------------------
# -----------------------------------------------------------
# ===========================================================

# Here are gathered all EVT related R functions

# OVERVIEW

# 1) gamma.function.estimator
#    gamma.estimator
#        - estimation of the Hill and moment estimator
#        - corresponding Hill plot and final gamma estimation

# 2) extreme1
#        - estimation of the right tail probability

# 3) tilde.Fn1.general
#        - estimation of the EVT modified cumulative empirical distribution function (one and two sided)

# 4) univariate.tukey.depth
#    univariate.simplicial.depth
#        - univariate Tukey depth and univariate simplicial depth

# 5) RTD.EVT_and_RPD
#        - estimation of the multivariate random Tukey depth
#        - estimation of the multivariate random Tukey depth with EVT modification
#        - estimation of the multivariate random projection depth
#        - estimation of the multivariate random projection depth with EVT modification

# 6) rcauchy
#    rclover
#    relliptical
#        - simulation of select bivariate extreme distribution functions


# ++++++++++++++++++++++++
# PART 1
# ++++++++++++++++++++++++
# ++++++++++++++++++++++++
# ++++++++++++++++++++++++
# ++++++++++++++++++++++++
# ++++++++++++++++++++++++


# Hill & Gamma estimators for range k, Hill-plot

# install.packages("fExtremes")
# library(fExtremes)
# library(graphics)



gamma.estimator <- function(x,k.range = (floor(length(x)/10)):(ceiling(length(x)/2)), hill.plot = T, estimate.hill = T, estimate.moment = T, legend = T, legend.position = "topleft", ylim = NULL, xlab = NULL, ylab = NULL, main = NULL){
# ===========================================================
# ------------------FUNCTION DESCRIPTION --------------------
# ===========================================================
# Extreme value index estimators applying the hill estimator and the moment estimator
# gamma.estimator(x, k.range, hill.plot, estimate.hill, estimate.moment, legend = T, legend.position = "topleft") 
# Input:
#  1) x
#     :: sample in R^1 for which extreme value index gamma is to be computed 
#  2) k.range = (floor(length(x)/10)):(ceiling(length(x)/2))
#     :: range of k's for which the gamma's are estimated
#     :: default is for 10% to 50% of data's lenght
#  3) hill.plot = T
#     :: plot hill plot, i.e. plot gamma w.r.t. to k; default is true
#  4) estimate.hill = T
#     :: estimate hill? default is true
#  5) estimate.moment = T
#     :: estimate moment? default is true
#  6) legend = T
#     :: add legend to graph? default is true  
#  7) position.legend = "topleft"
#     :: where the legend should be positioned  
# Output:
#  1) gamma.hill
#     :: hill estimator for the chosen range k.range (if estimate.hill == T)
#  2) gamma.moment
#     :: moment estimator for the chosen range k.range (if estimate.moment == T)
#  3) hill plot
#     :: hill plot for the chosen range k.range (if hill.plot == T)  
#     
#     WARNING: For some low k's or ranges of k, the estimators might not be defined

# ===========================================================
# ------------------BEGIN FUNCTION --------------------------
n1 <- length(x)
n.k.range <- length(k.range)

x.sorted <- sort(x)
gamma.hill <- gamma.moment <- rep(NA,n1)

for(k in k.range){

  M1 <- sum(log(x.sorted[(n1-k+1):n1]/x.sorted[n1-k]))/k
  M2 <- sum((log(x.sorted[(n1-k+1):n1]/x.sorted[n1-k]))^2)/k
  gamma.hill[k] <- M1
  gamma.moment[k] <- M1 + 1-0.5/(1-(M1)^2/M2)
}

gamma.hill <- gamma.hill[!is.na(gamma.hill)]
names(gamma.hill) <- k.range
gamma.moment <- gamma.moment[!is.na(gamma.moment)]
names(gamma.moment) <- k.range


if(is.null(ylim)){ylim <- range(c(gamma.hill,gamma.moment))}

if(hill.plot == T){
  if(estimate.hill == T & estimate.moment == T){
    plot(k.range, gamma.hill, type = "l", ylim = ylim, xlab = xlab, ylab = ylab, main = main)
    lines(k.range, gamma.moment, type = "l", col = "red")
    if(legend == T){
    legend(legend.position, c("Hill estimator", "Moment estimator"), col = c("black", "red"), lty = 1)
    }
  }

  if(estimate.hill == T & estimate.moment == F){
    plot(k.range, gamma.hill, type = "l", ylim = ylim, xlab = xlab, ylab = ylab, main = main)
  }
 
  if(estimate.hill == F & estimate.moment == T){
    plot(k.range, gamma.moment, type = "l", ylim = ylim, xlab = xlab, ylab = ylab, main = main)
  }   
}  
  
  if(estimate.hill == T & estimate.moment == T){
    return(list(gamma.hill = gamma.hill, gamma.moment = gamma.moment))
  }
   
  if(estimate.hill == T & estimate.moment == F){
    return(list(gamma.hill = gamma.hill))
  }
  
  if(estimate.hill == F & estimate.moment == T){
    return(list(gamma.moment = gamma.moment))
  }     
  
}



# ===========================================================
# ------------------ GAMMA FUNCTION ESTIMATOR ---------------
# ===========================================================

# -----------------------------------------------------------
# Estimation of gamma for each dimension of a stochastic process
# basically applying gamma.estimator to each dimension
# gamma.function.estimator(x, x.star = 1:ncol(x), chosen.k = rep(round(nrow(x))/4, ncol(x)), true.gamma, plot = T, estimate.moment = T, estimate.hill = T, legend = T, legend.position = "center")
# Input
#  1) x
#     :: data matrix for which the extreme value index function gamma is to be calculated.
#        dim = n.draws*n.values
#  2) x.star = 1:ncol(x)
#     :: values assigned to each dimension, default is 1:ncol(x)     
#  3) chosen.k = rep(ceiling(nrow(x)/4)
#     :: k value for which gamma is to be calculated for each dimension. default is 1/4 of n.draws for each dimension
#        WARNING: theoretically each dimension needs to be anlayzed with a hill plot, to dermine its corresponding k value
#  4) true.gamma
#     :: vector of length ncol(x) representing the true gamma function (for example gamma = rep(1, ncol(x))) which can be added if known     
#  5) plot = T
#     :: plot gamma function(s), i.e. plot gamma w.r.t. to each dimension; default is true
#  6) estimate.hill = T
#     :: estimate function with hill estimator? default is true
#  7) estimate.gamma = T
#     :: estimate function with moment estimator? default is true
#  8) legend = T
#     :: add legend to graph? default is true  
#  9) position.legend = "center"
#     :: where the legend should be positioned  
# Output:
#  1) gamma.function.hill
#     :: gamma function estimated with hill estimator (if estimate.hill == T)
#  2) gamma.function.moment
#     :: gamma function estimated with moment estimator (if estimate.moment == T)
#  3) plot
#     :: plot of gamma function(s) (if plot == T)        



gamma.function.estimator <- function(x, x.star = 1:ncol(x), chosen.k = ceiling(nrow(x)/4, ncol(x)), true.gamma, plot = T, estimate.moment = T, estimate.hill = T, legend = T, legend.position = "bottomleft"){
  gamma.function.hill <- gamma.function.moment <- rep(0, ncol(x))
  for (i in 1:ncol(x)){
    gamma.function.hill[i] <- mean(gamma.estimator(x[,i], chosen.k, estimate.hill = T, estimate.moment = F, hill.plot = F))
    gamma.function.moment[i] <- mean(gamma.estimator(x[,i], chosen.k, estimate.hill = F, estimate.moment = T, hill.plot = F))
  }

  if(plot == T){ #begin plot
    if(estimate.hill == T & estimate.moment == T){
      plot(x.star, gamma.function.hill, type = "l", ylim = range(gamma.function.hill, gamma.function.moment))
      lines(x.star, gamma.function.moment, col = "red")
      abline(h = mean(gamma.function.hill), lty = 2)
      abline(h = mean(gamma.function.moment), lty = 2, col = "red")
      if(!missing(true.gamma)){
        lines(x.star, true.gamma, col = "blue")
      }
      if(legend == T){
        if(!missing(true.gamma)){
          legend(legend.position, c("hill fct","moment fct","mean(hill)","mean(moment)","true gamma"), 
               col = c("black","red","black","red","blue"), lty = c(1,1,2,2,1))
        }else{
                  legend(legend.position, c("hill fct","moment fct","mean(hill)","mean(moment)"), 
               col = c("black","red","black","red"), lty = c(1,1,2,2))
        }
      }
    }
    
    if(estimate.hill == T & estimate.moment == F){
      plot(x.star, gamma.function.hill, type = "l")
      abline(h = mean(gamma.function.hill), lty = 2)
      if(!missing(true.gamma)){
        lines(x.star, true.gamma, col = "blue")
      }
      if(legend == T){
        if(!missing(true.gamma)){
          legend(legend.position, c("hill fct","mean(hill)","true gamma"), 
               col = c("black","black","blue"), lty = c(1,2,1))
        }else{
                  legend(legend.position, c("hill fct","mean(hill)"), 
               col = c("black","black"), lty = c(1,2))
        }
      }
    }    
    
    if(estimate.hill == F & estimate.moment == T){
      plot(x.star, gamma.function.moment, type = "l", col = "red")
      abline(h = mean(gamma.function.moment), lty = 2, col = "red")
      if(!missing(true.gamma)){
        lines(x.star, true.gamma, col = "blue")
      }
      if(legend == T){
        if(!missing(true.gamma)){
          legend(legend.position, c("moment fct","mean(moment)","true gamma"), 
               col = c("red","red","blue"), lty = c(1,2,1))
        }else{
                  legend(legend.position, c("moment fct","mean(moment)"), 
               col = c("red","red"), lty = c(1,2))
        }
      }
    }    
   } #end plot
  
  
  if(estimate.hill == T & estimate.moment == T){
    return(rbind(gamma.function.hill, gamma.function.moment))
  }

  if(estimate.hill == T & estimate.moment == F){
    return(gamma.function.hill)
  } 
 
  if(estimate.hill == F & estimate.moment == T){
    return(gamma.function.moment)
  }    
  
} # end function



# ++++++++++++++++++++++++
# PART 2
# ++++++++++++++++++++++++
# ++++++++++++++++++++++++
# ++++++++++++++++++++++++
# ++++++++++++++++++++++++
# ++++++++++++++++++++++++



# 1 - EVT cdf (tail function) w.r.t. to x0 at y given k and gamma (moment method)
extreme1 <- function(k,x0,y,gamma){
   x <- sort(x0)  
   n <- length(x)
   F.y <- 1-k/n*(y/x[n-k])^(-1/max(0.01,gamma)) # make sure gamma is positive
   return(F.y)
}



# ++++++++++++++++++++++++
# PART 3
# ++++++++++++++++++++++++
# ++++++++++++++++++++++++
# ++++++++++++++++++++++++
# ++++++++++++++++++++++++
# ++++++++++++++++++++++++


# tilde.Fn1.general is the EVT adjusted function of min(F(w[i]), 1-F(w[i])), using the distribution of x, for a chosen k and EVI gamma

tilde.Fn1.general <- function(w,x,k,gamma_right, gamma_left = gamma_right, one.sided = T){ # gamma_right is the right tail gamma estimator and gamma left is the left tail gamma estimator
      hd <- rep(0,n)
      x <- sort(x)
      m <- length(x)
      for (i in 1:n){
         if (w[i]>=x[m-k]) {  # right tail
            F.w <- extreme1(k,x,w[i],gamma_right)
            hd[i] <- min(F.w,1-F.w)
         } else if(w[i]<=x[k+1] & one.sided == F) {  # left tail (modified only if one.sided ==T)
		 x1 <- -x
		 w1 <- -w[i]
		 F.w <- extreme1(k,x1,w1,gamma_left)
  hd[i] <- min(F.w,1-F.w)
             }	else hd[i] <- min(sum(x<=w[i])/length(x),sum(x>=w[i])/length(x)) # standard empirical cdf in the middle of the distribution

      }
            return(hd)
}


# EVT modified ecdf using absolute values

EVT.two.sided.absolute <- function(w,x,k,gamma){ # in practice w represents projected vector `x' and x (data to calculate depth with) projected vector `y' (original data)
       n <- length(w)
      hd <- rep(0,n)
      x <- sort(x)
      x.sort.abs <- sort(abs(x))
      m <- length(x)
      for (i in 1:n){
         if (abs(w[i])> x.sort.abs[m-k]) {  # right tail
            F.w <- extreme1(k,x.sort.abs,abs(w[i]),gamma)
            hd[i] <- min(F.w, 1-F.w)
            
         }else hd[i] <- min(sum(x<=w[i])/length(x),sum(x>=w[i])/length(x)) # standard empirical cdf in the middle of the distribution

      }
            return(hd)
}







# ++++++++++++++++++++++++
# PART 4
# ++++++++++++++++++++++++
# ++++++++++++++++++++++++
# ++++++++++++++++++++++++
# ++++++++++++++++++++++++
# ++++++++++++++++++++++++


univariate.tukey.depth <- function(w,x){ # calculate univariate Tukey depth w.r.t to x
       n <- length(w)
      hd <- rep(0,n)
      x <- sort(x)
      m <- length(x)
      for (i in 1:n){
hd[i] <- min(sum(x<=w[i])/length(x),sum(x>=w[i])/length(x)) # standard empirical cdf in the middle of the distribution
      }
            return(hd)
}

univariate.simplicial.depth <- function(w,x){ # calculate univariate simplicial depth w w.r.t to x
       n <- length(w)
      hd <- rep(0,n)
      x <- sort(x)
      m <- length(x)
      for (i in 1:n){
hd[i] <- (sum(x<=w[i])/length(x))*(sum(x>=w[i])/length(x)) # standard empirical cdf in the middle of the distribution
      }
            return(hd)
}


# ++++++++++++++++++++++++
# PART 5
# ++++++++++++++++++++++++
# ++++++++++++++++++++++++
# ++++++++++++++++++++++++
# ++++++++++++++++++++++++
# ++++++++++++++++++++++++


# CALCULATE RPD, EVT modified RPD (though virtually identical to RPD), RTD and EVT modified RTD
# for multivariate data

RTD.EVT_and_RPD <- function(x,xx, n.proj = 100,k,gamma, one.sided = F){
        n <- dim(xx)[1]
        d <- dim(xx)[2]
        u <- matrix(runif(d*n.proj,-1,1),n.proj,d)
        norm <- sqrt(rowSums(u*u))
        arg <- u/norm
        z <- arg %*% t(xx)
        m <- length(x[, 1])
        z1 <- arg %*% t(x)
        out1 <- matrix(0, n.proj,m)
        out2 <- matrix(0, n.proj,m)
        for(i in 1:n.proj) { # for each projection do...
            out1[i,] <- tilde.Fn1.general(z1[i,], z[i,],k,gamma, one.sided = one.sided) 
                }
        RTD.EVT = as.vector(apply(out1,2,min,na.rm=TRUE)) # take the minimum over all projections
        RPD.EVT = as.vector(apply(out1,2,mean,na.rm=TRUE)) # same with mean
        for(i in 1:n.proj) { # for each projection do...
            out2[i,] <- univariate.tukey.depth(z1[i,], z[i,])
        }
        RPD = as.vector(apply(out2,2,mean,na.rm=TRUE))  # random projection depth
        RTD = as.vector(apply(out2,2,min,na.rm=TRUE))   # traditional random Tukey depth
        results = list(RPD = RPD, RPD.EVT = RPD.EVT, RTD = RTD, RTD.EVT = RTD.EVT)
 return(results)
}




# ++++++++++++++++++++++++
# PART 6
# ++++++++++++++++++++++++
# ++++++++++++++++++++++++
# ++++++++++++++++++++++++
# ++++++++++++++++++++++++
# ++++++++++++++++++++++++

rcauchy = function(n) {
  r = sqrt(runif(n)^(-2)-1)
  theta = (2*pi) * runif(n)
  #Z = array(c(r * cos(theta), r * sin(theta)), dim = c(n, 2))
  #Z
  array(c(r ,theta), dim = c(n, 2))
}


rclover = function(n){
  U = runif(n)
  V = runif(n)
  R = (U^(-2)-1)^(1/6)
  Theta = rep(0,n)
  
for(i in 1:n){
    # f = function(t){(5*t+sin(4*t))/(10*pi)-runif(1)}
    f = function(t){((3*R[i])/(10*pi*((1+R[i]^6)^(3/2))))*(5*t+sin(4*t))/(3*R[i]/(1+R[i]^6)^(3/2))-V[i]}
    Theta[i]=uniroot(f,lower=0,upper=2*pi)$root
}
   # cbind(R*cos(Theta), R*sin(Theta))
    cbind(R,Theta)
}



relliptical = function(n) {
  r = (runif(n)^(-2)-1)^(1/6)
  theta = (2*pi)*runif(n)
  Z = array(c(2*r*cos(theta), r*sin(theta)), dim = c(n, 2))
  Z
}










# ++++++++++++++++++++++++
# BONUS
# ++++++++++++++++++++++++
# ++++++++++++++++++++++++
# ++++++++++++++++++++++++
# ++++++++++++++++++++++++
# ++++++++++++++++++++++++


# BONUS FUNCTION  | PROBABLY NOT OF VERY IMPORTANT USE



#-----------------------------------------
# frechet function with evi gamma
frechet <- function(x,gamma){
  
  f.value <- rep(0,length(x))
  for(i in 1:length(x))
    
    if(1+gamma*x[i] > 0){f.value[i] <- exp(-(1+gamma*x[i])^(-1/gamma))}else{f.value[i] <- 0}
  return(f.value)
}
