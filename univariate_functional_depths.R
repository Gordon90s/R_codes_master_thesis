# ===========================================================
# ------ UNIVARIATE FUNCTIONAL DEPTH FUNCTIONS --------------
# ===========================================================


# This R code groups all written FUNCTIONAL depth functions during this thesis
# Because all functions have not been implemented all at once and before finishing the thesis,
# the notation of the original papers is often used
# making it, unfortunately, highly inconsistant over all functions

# Here is the overview of all RELEVANT functional depth functions available in this R document:

# OVERVIEW

# 1) main.depth.function
#   - Band depth [computational time increasing fast for larger n]
#   - Modified band depth [computational time increasing fast for larger n, use integrated depth with simplicial marginals]
#   - Half-region depth
#   - Modified half-region depth

# WARNING: REQUIRES INSTALLATION OF (YET IMCOMPLETE) R PACKAGE ``ldfun'' by Agostinelli (2015) ``Local Half-Region Depth for Functional Data''


# 2) MFD    (originally standing for multivariate functional depth)
#   - Integrated depth (univariate only) with and without weight function with Tukey or simplicial marginals
#   - Integrated depth (univariate only) with EVT modification for the Tukey marginals (though virtually not changing anything to the final depth values)
#   - Extremal depth

# 3) MFRPD   (originally standing for multivariate functional random projection depth, in the end only univariate implementation)
#   - functional random projection depth with mean
#   - functional random projection depth with minimum
#   - functional random projection depth with mean & EVT (though virtually not changing anything to the final depth values)
#   - functional random projection depth with minimum & EVT

# More details can be gained from the detailed comments for each R function bellow

# ++++++++++++++++++++++++
# PART 1
# ++++++++++++++++++++++++
# ++++++++++++++++++++++++
# ++++++++++++++++++++++++
# ++++++++++++++++++++++++
# ++++++++++++++++++++++++


# ===========================================================
# -----------------------------------------------------------
# --------------FIRST IMPLEMENTATIONS -----------------------
# -----------------------------------------------------------
# ===========================================================


# -----------------------------------------------------------
# main.depth.function(data, depth.type = 1, tau, save.ordered.data = T, graph = T)
# -----------------------------------------------------------
# Input: 
#  1) x 
#     :: functional data set with which the depths are calculated
#  2) y = NULL
#     :: functional data the depths are calculated for. By default x is taken and y = NULL
#  3) depth.type = type.number  
#     :: choose function depth type
#        type.number in {1,...,9}
#        1 = Random Tukey Depth / now called functional random projection depth with marginal Depth = Tukey depth (depending on the paper also called functional Random Tukey depth or double random projection depth)
#        2 = Band Depth [computational time increasing fast for larger n]
#        3 = Local Band Depth [ignoring for now]
#        4 = Modified Band Depth [computational time increasing fast for larger n, use integrated depth with simplicial marginals]
#        5 = Modified Local Band Depth [ignoring for now]
#        6 = Half Region Depth
#        7 = Local Half Region Depth [ignoring for now]
#        8 = Modified Half Region Depth
#        9 = Modified Local Half Region Depth [ignoring for now]


#  4) tau
#       :: for local band depths, default is the 30%-quantile of the supnorm 
#          between the curves (see more in Agostinelli (2015) p29)
#          the local depths by Agostinelli always calculate depths and local depths.
#  5) save.ordered.data = T 
#       :: save ordered data? True or False
#  6) graph = T 
#       :: graph wanted? True or False
#  7) nproj
#       :: number of random projections for the RTD      

# Output:
#  1) $depth
#       :: depth using one of 10 available types
#  2) $normalized.depth
#       :: normalized depth
#  3) $order
#       :: data ordering by depth
#  4) $depth.type
#       :: depth type used
#  5) $sorted.data (if save.ordered.data == T)
#       :: ordered data via depth
#  6) colored graph by depth (if graph == T)
#       :: Black = extreme, white = center, red = deepest curve

# Assumption:
#  1) no missing values, functions saved at equidistent points





# ===========================================================
# ----------------FUNCTION IMPLEMENTATION -------------------
# ===========================================================

# normalize function
normalize <- function(x) {
  x <- (x - min(x))/(max(x) - min(x))
  return(x)
}

# main depth function
main.depth.function <-function(data, y = NULL, depth.type = 1, tau = NULL, save.ordered.data = T, graph = T, nproj = 10){
  
was.y.null <- F # initializing checking variable
  
# was y set to NULL?  
  if(is.null(y)){
    y <- data # if yes, set y = data
    was.y.null <- T
  }
  
  if(depth.type == 1){
    
    if(was.y.null == F){print("ERROR, depth.RT function can't handle y != x (yet?)")}
    
    fdatobj <- fdata(data)
    fdataobj.RT <- depth.RT(fdatobj, nproj = nproj) # calculate Tukey depth
    fdataobj.RT.normalized <- normalize(fdataobj.RT$dep)
    order <- order(fdataobj.RT$dep) # order depth
    sorted.data <- data[order,] # sort data according to depth
    
    result <- list(depth = fdataobj.RT$dep, normalized.depth = fdataobj.RT.normalized, order = order, depth.type = depth.type, sorted.data = sorted.data)
    
    if(graph == T){
      col <- (palette(gray(seq(0.2,1,len = nrow(y))))) # 50 shades of grey
     # col <- magent(nrow(y))
      # initialize plot
      i <- 1
      plot(sorted.data[i,], type = "l", ylim = c(min(sorted.data), max(sorted.data)), xlab = "t", ylab = "x(t)")
      # plot with different greys, deeper -> white, less deep -> black
      for(i in 1:nrow(sorted.data)){
      lines(sorted.data[i,], col = col[i])
      }
     #lines(fdataobj.RT$median, col = "red")
    }  
  }
  
  if(depth.type == 2){
    
    if (is.null(tau)) {
      tau <- quantile.localdepth.functional(data, probs = 0.3)
    }
    
    # calculate modified half region depth
    data.localdepth.band <- localdepth.band(data, y, tau = tau)
    # the higher tau is chosen, the more the local depth becomes global
    normalized.depth <- normalize(data.localdepth.band$depth) # normalize depth
    order <- order(data.localdepth.band$depth) # order depth
    sorted.data <- y[order,] # sort data according to depth
    
    result <- list(depth = data.localdepth.band$depth, normalized.depth = normalized.depth, order = order, depth.type = depth.type, sorted.data = sorted.data)
    
    if(graph == T){
      col <- (palette(gray(seq(0.1,1,len = nrow(y))))) # 50 shades of gray
      
      # initialize plot
      i <- 1
      plot(sorted.data[i,], type = "l", ylim = c(min(sorted.data), max(sorted.data)), xlab = "x", ylab = "y", main = "BD")
      # plot with different grays, deeper -> white, less deep -> black
      for(i in 1:nrow(sorted.data)){
        lines(sorted.data[i,], col = col[i])
      }
      
      #plot deepest curve
      lines(sorted.data[nrow(sorted.data),], col = "red")
    }
  } 
  
  if(depth.type == 3){
    
    if (is.null(tau)) {
      tau <- quantile.localdepth.functional(data, probs = 0.3)
    }
    
    # calculate modified half region depth
    data.localdepth.band <- localdepth.band(data, y, tau = tau)
    # the higher tau is chosen, the more the local depth becomes global
    normalized.depth <- normalize(data.localdepth.band$localdepth) # normalize depth
    order <- order(data.localdepth.band$localdepth) # order depth
    sorted.data <- y[order,] # sort data according to depth
    
    result <- list(depth = data.localdepth.band$localdepth, normalized.depth = normalized.depth, order = order, depth.type = depth.type, sorted.data = sorted.data)
    
    if(graph == T){
      col <- (palette(gray(seq(0.1,1,len = nrow(y))))) # 50 shades of gray
      
      # initialize plot
      i <- 1
      plot(sorted.data[i,], type = "l", ylim = c(min(sorted.data), max(sorted.data)), xlab = "x", ylab = "y", main = "LBD")
      # plot with different grays, deeper -> white, less deep -> black
      for(i in 1:nrow(sorted.data)){
        lines(sorted.data[i,], col = col[i])
      }
      
      #plot deepest curve
      lines(sorted.data[nrow(sorted.data),], col = "red")
    }
  }   
  
  if(depth.type == 4){
    
    if (is.null(tau)) {
      tau <- quantile.localdepth.functional(data, probs = 0.3)
    }
    
    # calculate modified half region depth
    data.localdepth.modband <- localdepth.modband(data, y, tau = tau)
    # the higher tau is chosen, the more the local depth becomes global
    normalized.depth <- normalize(data.localdepth.modband$depth) # normalize depth
    order <- order(data.localdepth.modband$depth) # order depth
    sorted.data <- y[order,] # sort data according to depth
    
    result <- list(depth = data.localdepth.modband$depth, normalized.depth = normalized.depth, order = order, depth.type = depth.type, sorted.data = sorted.data)
    
    if(graph == T){
      col <- (palette(gray(seq(0.1,1,len = nrow(y))))) # 50 shades of gray
      
      # initialize plot
      i <- 1
      plot(sorted.data[i,], type = "l", ylim = c(min(sorted.data), max(sorted.data)), xlab = "x", ylab = "y", main = "MBD")
      # plot with different grays, deeper -> white, less deep -> black
      for(i in 1:nrow(sorted.data)){
        lines(sorted.data[i,], col = col[i])
      }
      
      #plot deepest curve
      lines(sorted.data[nrow(sorted.data),], col = "red")
    }
  }   
  
  if(depth.type == 5){
    
    if (is.null(tau)) {
      tau <- quantile.localdepth.functional(data, probs = 0.3)
    }
    
    # calculate modified half region depth
    data.localdepth.modband <- localdepth.modband(data, y, tau = tau)
    # the higher tau is chosen, the more the local depth becomes global
    normalized.depth <- normalize(data.localdepth.modband$localdepth) # normalize depth
    order <- order(data.localdepth.modband$localdepth) # order depth
    sorted.data <- y[order,] # sort data according to depth
    
    result <- list(depth = data.localdepth.modband$localdepth, normalized.depth = normalized.depth, order = order, depth.type = depth.type, sorted.data = sorted.data)
    
    if(graph == T){
      col <- (palette(gray(seq(0.1,1,len = nrow(y))))) # 50 shades of gray
      
      # initialize plot
      i <- 1
      plot(sorted.data[i,], type = "l", ylim = c(min(sorted.data), max(sorted.data)), xlab = "x", ylab = "y", main = "LMBD")
      # plot with different grays, deeper -> white, less deep -> black
      for(i in 1:nrow(sorted.data)){
        lines(sorted.data[i,], col = col[i])
      }
      
      #plot deepest curve
      lines(sorted.data[nrow(sorted.data),], col = "red")
    }
  }   
  
  if(depth.type == 6){
    
    if (is.null(tau)) {
      tau <- quantile.localdepth.functional(data, probs = 0.3)
    }
    
    # calculate modified half region depth
    data.localdepth.halfregion <- localdepth.halfregion(data, y, tau = tau)
    # the higher tau is chosen, the more the local depth becomes global
    normalized.depth <- normalize(data.localdepth.halfregion$depth) # normalize depth
    order <- order(data.localdepth.halfregion$depth) # order depth
    sorted.data <- y[order,] # sort data according to depth
    
    result <- list(depth = data.localdepth.halfregion$depth, normalized.depth = normalized.depth, order = order, depth.type = depth.type, sorted.data = sorted.data)
    
    if(graph == T){
      col <- (palette(gray(seq(0.1,1,len = nrow(y))))) # 50 shades of gray
      
      # initialize plot
      i <- 1
      plot(sorted.data[i,], type = "l", ylim = c(min(sorted.data), max(sorted.data)), xlab = "x", ylab = "y", main = "HD")
      # plot with different grays, deeper -> white, less deep -> black
      for(i in 1:nrow(sorted.data)){
        lines(sorted.data[i,], col = col[i])
      }
      
      #plot deepest curve
      lines(sorted.data[nrow(sorted.data),], col = "red")
    }
  } 
  
  if(depth.type == 7){
    
    if (is.null(tau)) {
      tau <- quantile.localdepth.functional(data, probs = 0.3)
    }
    
    # calculate modified half region depth
    data.localdepth.halfregion <- localdepth.halfregion(data, y, tau = tau)
    # the higher tau is chosen, the more the local depth becomes global
    normalized.depth <- normalize(data.localdepth.halfregion$localdepth) # normalize depth
    order <- order(data.localdepth.halfregion$localdepth) # order depth
    sorted.data <- y[order,] # sort data according to depth
    
    result <- list(depth = data.localdepth.halfregion$localdepth, normalized.depth = normalized.depth, order = order, depth.type = depth.type, sorted.data = sorted.data)
    
    if(graph == T){
      col <- (palette(gray(seq(0.1,1,len = nrow(y))))) # 50 shades of gray
      
      # initialize plot
      i <- 1
      plot(sorted.data[i,], type = "l", ylim = c(min(sorted.data), max(sorted.data)), xlab = "x", ylab = "y", main = "LHD")
      # plot with different grays, deeper -> white, less deep -> black
      for(i in 1:nrow(sorted.data)){
        lines(sorted.data[i,], col = col[i])
      }
      
      #plot deepest curve
      lines(sorted.data[nrow(sorted.data),], col = "red")
    }
  } 
   
  if(depth.type == 8){
    
    if (is.null(tau)) {
      tau <- quantile.localdepth.functional(data, probs = 0.3)
    }
    
    # calculate modified half region depth
    data.localdepth.modhalfregion <- localdepth.modhalfregion(data, y, tau = tau)
    # the higher tau is chosen, the more the local depth becomes global
    normalized.depth <- normalize(data.localdepth.modhalfregion$depth) # normalize depth
    order <- order(data.localdepth.modhalfregion$depth) # order depth
    sorted.data <- y[order,] # sort data according to depth
    
    result <- list(depth = data.localdepth.modhalfregion$depth, normalized.depth = normalized.depth, order = order, depth.type = depth.type, sorted.data = sorted.data)
    
    if(graph == T){
     col <- (palette(gray(seq(0.1,1,len = nrow(y))))) # 50 shades of gray
     
     # initialize plot
     i <- 1
      plot(sorted.data[i,], type = "l", ylim = c(min(sorted.data), max(sorted.data)), xlab = "x", ylab = "y", main = "MHRD")
     # plot with different grays, deeper -> white, less deep -> black
        for(i in 1:nrow(sorted.data)){
       lines(sorted.data[i,], col = col[i])
        }
    
     #plot deepest curve
      lines(sorted.data[nrow(sorted.data),], col = "red")
     }
  } 
  
  if(depth.type == 9){
    
    if (is.null(tau)) {
      tau <- quantile.localdepth.functional(data, probs = 0.3)
    }
    
    # calculate modified half region depth
    data.localdepth.modhalfregion <- localdepth.modhalfregion(data, y, tau = tau)
    # the higher tau is chosen, the more the local depth becomes global
    normalized.depth <- normalize(data.localdepth.modhalfregion$localdepth) # normalize depth
    order <- order(data.localdepth.modhalfregion$localdepth) # order depth
    sorted.data <- y[order,] # sort data according to depth
    
    result <- list(depth = data.localdepth.modhalfregion$localdepth, normalized.depth = normalized.depth, order = order, depth.type = depth.type, sorted.data = sorted.data)
    
    if(graph == T){
      col <- (palette(gray(seq(0.1,1,len = nrow(y))))) # 50 shades of gray
      
      # initialize plot
      i <- 1
      plot(sorted.data[i,], type = "l", ylim = c(min(sorted.data), max(sorted.data)), xlab = "x", ylab = "y", main = "LMHRD")
      # plot with different grays, deeper -> white, less deep -> black
      for(i in 1:nrow(sorted.data)){
        lines(sorted.data[i,], col = col[i])
      }
      
      #plot deepest curve
      lines(sorted.data[nrow(sorted.data),], col = "red")
    }
  }   
  
  


return(result)    
}


# ++++++++++++++++++++++++
# PART 2
# ++++++++++++++++++++++++
# ++++++++++++++++++++++++
# ++++++++++++++++++++++++
# ++++++++++++++++++++++++
# ++++++++++++++++++++++++


# ===========================================================
# -----------------------------------------------------------
# INTEGRATED DEPTH WITH WEIGHT FUNCTION, EVT MODIFICATION and EXTREMAL DEPTH
# -----------------------------------------------------------
# ===========================================================

# -----------------------------------------------------------
# MFD(y, x = y, x.star = 1:(dim(y)[2]), K = 1, marginal.depth = "random Tukey", n.projections = 500, weight = "alpha volume", 
#                      range = 1:n.values, alpha.coverage = 0.25, gamma.fct, chosen.k = rep(round(nrow(x))/4, ncol(x)), estimate.hill = T, ED = F)
# -----------------------------------------------------------
# Input: 
#  1) y 
#     :: multivariate functional data with which the depths are to be calculated, 3-dimensional array dim= c(n.draws, n.values, K) 
#  2) x = y
#     :: multivariate functional data for which the depths are to be calculated, 3-dimensional array dim= c(n.draws, n.values, K)  
#        default is x = y
#  3) T = 1:(dim(y)[2])
#     :: time points t_1 < t_2 < ... < t_{dim(y)[2])} at which the n.values are evaluted, default is just 1:(dim(y)[2])   
#  4) marginal.depth = "random Tukey"
#     :: depth chosen for the marginal depth, default is the "random Tukey"
#        other possibility is "EVT random Tukey", the EVT refined random Tukey depth, explained in detail in the master thesis
#        or "simplicial" for the simplicial depth
#  5) one.sided.EVT = T
#     :: if(marginal.depth == "EVT random Tukey") is selected, then either have both tails EVT refined with one.sided.EVT = F
#        or only the right tail with one.sided.EVT = T (the default value)  
#  6) constant.gamma = F
#     :: if(marginal.depth == "EVT random Tukey") is selected, then one can decide whether one wants to assume the same
#        extreme value index gamma for all marginals or not. Default is false.
#  7) n.projections = 500
#     :: if chosen depth is "random Tukey" and if K>1 then n.projections gives the number or random projections that are carried out 
#        in the algorithm. Default is 500, which is a 'good' value for K=2
#  8) weight = "alpha volume"
#     :: the MFD allows for different weight functions for the integral (c.f. depth formula). Available are 
#        i) "constant"
#          :: constant weighting over a specified range cf. 6)
#        ii) "alpha volume"
#          :: the integral is weighted with the volume obtained by marginal.depth >= alpha.coverage cf. depth formula and 7) for more infos 
#        default is "alpha volume"    
#  9) range = 1:n.values
#     :: range over which the functional depth is to be calculed [NOT IMPLEMENTED YET]
#        default is over all n.values          
# 10) alpha.coverage = 0.25
#     :: parameter with which the weight of the integral is determined. A value of 0.25 represents an approximate 50% coverage for univariate normal data
#        and a value of 1/8 represents an approximate 50% coverage for bivariate normal data. These coverages are only for the marginals. [QUESTION CONCERNING FUNCTIONAL COVERAGE]
# 11) gamma.fct
#     :: value of the gamma function for each dimension 1:ncol(x). If not provided, it is automatically calculated using gamma.function.estimator function and chosen.k
# 12) chosen.k = rep(round(nrow(y))/4, ncol(y))
#     :: k value for which gamma is to be calculated for each dimension. default is 1/4 of n.draws for each dimension
#        WARNING: theoretically each dimension needs to be anlayzed with a hill plot, to dermine its corresponding k value   
# 13) estimate.hill = T
#     :: estimate gamma function using hill estimator if estimate.hill == T, moment estimator otherwise    
# 14) ED = F
#     :: if ED = T, then the extremal depth is calculated. Weight functions are ignored if ED = T.        
# Output:
#  1) function.depth.x
#     :: depth values for all functions of data set x
#  2) alpha.range
#     :: alpha functional band obtained which which the depth function was weighted if weight = "alpha volume"  
#  3) depth.x
#     :: marginal depth for each point of the data matrix x
#  4) marginal.depth
#     :: marginal depth decided upon in Input 4)  


# --------------------------------
# BEGIN DEBUG
# ----------------------------------

# BR.data <- smooth.gp.data
# BR.data <- parallel.data
# 
# # test x & y's
# y <- BR.data
# x <- BR.data[1:(nrow(BR.data)/2),]
# 
# 
# 
# 
# 
# # test x = y
# y <- x <- BR.data
# 
# 
# # test x & y's
# y <- array(c(1:20,11:30), dim = c(20,2))
# x <- array(1:20, dim = c(10,2))
# 
# # test x & y's
# y <- array(c(1:20,41:50), dim = c(20,2))
# x <- array(1:20, dim = c(10,2))
# 
# 
# x <- rbind(1:200,1:200)
# y <- x
# 
# alpha.coverage = 0.25
# K = 1
# x.star = 1:(dim(y)[2])
# weight = "alpha volume"
# marginal.depth = "EVT random Tukey"


# x <- y <- data.array[,,1]

# --------------------------------
# END DEBUG
# ----------------------------------

# WARNING: Works only for K = 1
MDF <- function(y, x = y, x.star = 1:(dim(y)[2]), K = 1, marginal.depth = "random Tukey", one.sided.EVT = T, constant.gamma = F, 
                n.projections = 500, weight = "alpha volume", range = 1:(dim(y)[2]), alpha.coverage = 0.25,
                gamma.fct, chosen.k = rep(ceiling(nrow(y)/4), ncol(y)), estimate.hill = T, ED = F){
  
if(K == 1){ 
    
n1.x <- dim(x)[1]
n2.x <- dim(x)[2]
n1.y <- dim(y)[1]
n2.y <- dim(y)[2]

# add more warnings??
  if(n2.x != n2.y){
    warning("Error: matrix dimensions not appropriate!!")
  }

# calculate marginal depth for each point in data matrix applying the random Tukey depth
if(marginal.depth == "random Tukey"){
  # corresponding marginal depth matrix for x
  depth.x <- matrix(0, nrow = n1.x, ncol = n2.x)
    for(j in 1:n2.x){
      for(i in 1:n1.x){  
       depth.x[i,j] <- depth(x[i,j],y[,j])
      }
    }
  # corresponding marginal depth matrix for y (to calculate alpha-volume weight function)
  depth.y <- matrix(0, nrow = n1.y, ncol = n2.y)
  for(j in 1:n2.y){
    for(i in 1:n1.y){  
      depth.y[i,j] <- depth(y[i,j],y[,j])
    }
  }
}


# calculate marginal depth for each point in data matrix applying the random Tukey depth
if(marginal.depth == "simplicial"){
  # corresponding marginal depth matrix for x
  depth.x <- matrix(0, nrow = n1.x, ncol = n2.x)
    for(j in 1:n2.x){
      for(i in 1:n1.x){
       depth.x[i,j] <- univariate.simplicial.depth(x[i,j],y[,j])
      }
    }
  # corresponding marginal depth matrix for y (to calculate alpha-volume weight function)
  depth.y <- matrix(0, nrow = n1.y, ncol = n2.y)
  for(j in 1:n2.y){
    for(i in 1:n1.y){  
      depth.y[i,j] <- univariate.simplicial.depth(y[i,j],y[,j])
    }
  }
}



if(marginal.depth == "EVT random Tukey"){
  # calculate marginal depth for each point in data matrix applying the EVT refined random Tukey depth
  if(missing(gamma.fct)){
  gamma.fct <- gamma.function.estimator(y, x.star, chosen.k, estimate.hill = estimate.hill, estimate.moment = !estimate.hill, plot = F)
  }
  
  if(constant.gamma == F){
    depth.x <- matrix(0, nrow = n1.x, ncol = n2.x)
    for(j in 1:n2.x){
       depth.x[,j] <- tilde.Fn1.general(x[,j],y[,j],chosen.k[j],gamma.fct[j], one.sided = one.sided.EVT)
    }
    # corresponding marginal depth matrix for y (to calculate alpha-volume weight function)
   depth.y <- matrix(0, nrow = n1.y, ncol = n2.y)
   for(j in 1:n2.y){
       depth.y[,j] <- tilde.Fn1.general(y[,j],y[,j],chosen.k[j],gamma.fct[j], one.sided = one.sided.EVT)
    }  
  } # end constant.gamma
  
  if(constant.gamma == T){
    gamma.fct.constant <- mean(gamma.fct)
    depth.x <- matrix(0, nrow = n1.x, ncol = n2.x)
    for(j in 1:n2.x){
       depth.x[,j] <- tilde.Fn1.general(x[,j],y[,j],chosen.k[j],gamma.fct.constant, one.sided = one.sided.EVT)
    }
    # corresponding marginal depth matrix for y (to calculate alpha-volume weight function)
   depth.y <- matrix(0, nrow = n1.y, ncol = n2.y)
   for(j in 1:n2.y){
       depth.y[,j] <- tilde.Fn1.general(y[,j],y[,j],chosen.k[j],gamma.fct.constant, one.sided = one.sided.EVT)
    }  
  } # end constant.gamma  
}

if(marginal.depth == "Frechet(1)"){
  depth.x <- matrix(0, nrow = n1.x, ncol = n2.x)
  for(j in 1:n2.x){
    depth.x[,j] <- frechet1.tilde(x[,j])
  }
}

if(weight == "alpha volume"){
# weight function alpha-volume    
depth.larger.alpha <- matrix(0, nrow = n1.y, ncol = n2.y)
depth.larger.alpha <- depth.y >= alpha.coverage
alpha.range <- matrix(0, nrow = 2, ncol = n2.y)
w_j.m <- rep(0, n2.y) # initilize marginal weight

# for each marginal j, get the corresponding range
for(j in 1:n2.y){
  alpha.range[,j] <- range(y[depth.larger.alpha[,j],j])
  if(j == 1){
    w_j.m[j] <- diff(alpha.range[,j])*(x.star[j+1]-x.star[j])
  }else if(j == n2.y){
    w_j.m[j] <- diff(alpha.range[,j])*(x.star[j]-x.star[j-1])
  }else{
    w_j.m[j] <- diff(alpha.range[,j])*(x.star[j+1]-x.star[j-1])
  }
}
w_j.alpha <- w_j.m/sum(w_j.m)
} # end if weight

if(weight == "constant"){
# constant weight function
w_j.contant.m <- rep(0, n2.y) # initilize marginal weight
for(j in 1:n2.y){
  if(j == 1){
    w_j.contant.m[j] <- (x.star[j+1]-x.star[j])
  }else if(j == n2.y){
    w_j.contant.m[j] <- (x.star[j]-x.star[j-1])
  }else{
    w_j.contant.m[j] <- (x.star[j+1]-x.star[j-1])
  }
}
w_j.constant <- w_j.contant.m/sum(w_j.contant.m)    
}

if(weight == "alpha volume"){
  w_j <- w_j.alpha
} else if(weight == "constant"){
  w_j <- w_j.constant
}

function.depth.x <- rep(0, n1.x)
for(i in 1:n1.x){
  function.depth.x[i] <- sum(depth.x[i,]*w_j)
  }
} # end if(K == 1)
  
  
# EXTREMAL DEPTH  
  
if(ED == T){
  if(marginal.depth == "random Tukey"){
    
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

    
    
    } # end random Tukey
  
  
  if(marginal.depth == "EVT random Tukey"){
    r.n <- 1000 # the bigger r.n, the more refined r will be
    r <- 0:r.n/r.n
  } # end EVT random Tukey
 
} # end extremal depth
  
  
  



if(ED == F){    
  if(weight == "alpha volume"){
  results <- list(function.depth.x = function.depth.x, alpha.range = alpha.range, depth.x = depth.x, marginal.depth = marginal.depth)  
  }else if(weight == "constant"){
  results <- list(function.depth.x = function.depth.x, depth.x = depth.x, marginal.depth.used = marginal.depth)  
  }
} # End ED == F  

if(ED == T){
  results <- list(function.depth.x = function.depth.x, Phi.x.r = Phi.x.r)
  
}
  
  
return(results)
  
} # end MFD function
  

  
# Live tests FOR DEBUG K>2
# MDF(data,y, K=2)




# ++++++++++++++++++++++++
# PART 3
# ++++++++++++++++++++++++
# ++++++++++++++++++++++++
# ++++++++++++++++++++++++
# ++++++++++++++++++++++++
# ++++++++++++++++++++++++




# =============================================================
# ------- MULTIVARIATE FUNCTIONAL PROJECTION DEPTH ------------
# =============================================================


# WARNING: IMPLEMENTED ONLY FOR K = 1

# -----------------------------------------------------------
# Input: 
#  1) data
#     :: only univariate functional data atm. no ``y'' available at the moment, i.e. y = data for now)
#  2) x.star
#     :: discretization points  
#  3) n.proj
#     :: number of Gaussian random projections, for more information see the function ``rproc2fdata'' of the package ``fda.usc''. Default is n.proj = 1000
#  4) proj.type = "normal"
#     :: for EVT analysis, it is sometimes important to choose projections that only have positive values
#     :: possibilities are 
#     :: "normal" for traditional Gaussan random projections or
#     :: "abs" for absolute values taken on the projections
#  5) method
#     :: choice between "mean" or "min". Default is "mean"
#  6) EVT = T
#     :: calculate EVT modifications
#  7) range.kappa.estimation
#     :: range over which gamma is estimated. Not trivial, see EVT in general for more informations.
#  8) set.seed.value
#     :: if seed is seeked. Default is set.seed.value = NULL 
#  9) more.than.one.depth
#     :: bonus feature, default is more.than.one.depth = NULL
#     :: you can try the following 3 combinations
#     ::   more.than.one.depth == 2 & EVT == F
#     ::   more.than.one.depth == 2 & EVT == T  
#     ::   more.than.one.depth == 4 & EVT == T  
#  10) constant.gamma = T
#     :: assume that gamma is the same for every projection v_i
#  11) one.sided = T
#     :: modify the ecdf with EVT only on the right side

  
# Output  
#  1) depth
#     :: depth values for data, saved as list.
#     :: if more.than.one.depth != NULL, then the function possibly returns more than one depth. See ``more.than.one.depth'' above.  


MFRPD <- function(data, x.star = 1:ncol(data), n.proj = 1000, proj.type = "normal", method = "mean", 
                  EVT = F, range.kappa.estimation = (round(nrow(data)/10)):(3*round(nrow(data)/10)), set.seed.value = NULL, 
                  more.than.one.depth = 0){

  set.seed(set.seed.value) # sets seed if seeked
  
  n.draws <- nrow(data)
  n.values <- ncol(data)
  
  k <- n.proj # number of `functional' random projections
  proj.functions <- rproc2fdata(k, t = x.star) # regular Gaussian processes
  
  if(proj.type == "normal"){
  proj <- proj.functions$data
  }else if(proj.type == "abs"){
  proj <- abs(proj.functions$data)
  }

  # calculate inner products
  scalars <- array(0, dim = c(n.draws, k)) 
  for(i in 1:n.draws){
    for(j in 1:k){
      scalars[i,j] <- mean(data[i,]*proj[j,])
    }                      
  }
  
  if(EVT == F | more.than.one.depth == 4){
  # apply univariate Tukey depth to inner products
  depth.k <- matrix(0, nrow = n.draws, ncol = k)
  for(j in 1:k){
    depth.k[,j] <- univariate.tukey.depth(scalars[,j],scalars[,j])
  }
  
  FRPD.depth.mean <- apply(depth.k, 1, mean) # calculate final functional random projection depth with mean 
  FRPD.depth.min <- apply(depth.k, 1, min) # calculate final functional random projection depth with min
  } 
  
  if(EVT == T){
    # calculate for all k inner products with EVT modification
    gamma.estimator.fct <- matrix(NA, nrow = k, ncol = max(range.kappa.estimation))
    for(j in 1:k){
      
      # taking absolute values as data is symmetric around zero
      # for other see for example using
      # apply(scalars,2,median)
      gamma.estimator.fct[j,range.kappa.estimation] <- gamma.estimator(abs(scalars[,j]),k.range = range.kappa.estimation, legend = F, 
                                                            legend.position = "bottomright", hill.plot = F, 
                                                            main = paste0("k = ",j), xlab = expression(kappa), estimate.hill = T, 
                                                            estimate.moment = F, ylim = c(0.45,1.45))$gamma.hill
    }
    
    
    # estimate gamma for each k
    gamma.estimator.final <- apply(gamma.estimator.fct[,range.kappa.estimation],1,mean)
    # estimate gamma by averaging over all gamma_k
    gamma.estimator.final.mean <- mean(gamma.estimator.final)
    
    # calculate EVT with variable gamma over k
    if(more.than.one.depth == 4){
      TD.EVT.marginal <- matrix(0,nrow = n.draws, ncol = k)
      for(j in 1:k){
        TD.EVT.marginal[,j] <- EVT.two.sided.absolute(scalars[,j],scalars[,j],
                                                 round(mean(range.kappa.estimation)), # what final kappa to use
                                                 gamma.estimator.final[j]) # which gamma to use
      }
      FRPD.EVT.depth.mean.variable.gamma <- apply(TD.EVT.marginal, 1, mean) # calculate final functional random projection depth with EVT modification with mean
      FRPD.EVT.depth.min.variable.gamma <- apply(TD.EVT.marginal, 1, min) # calculate final functional random projection depth with EVT modification with min
    }
    
    # calculate EVT with constant gamma over k
    # if(constant.gamma == T | more.than.one.depth == 6){
    # TD.EVT.marginal.constant.gamma <- matrix(0,nrow = n.draws, ncol = k)
    # for(j in 1:k){
    #   TD.EVT.marginal.constant.gamma[,j] <- EVT.two.sided.absolute(scalars[,j],scalars[,j],
    #                                            round(mean(range.kappa.estimation)), # what final kappa to use, factor 1/2 comes because two sided (k/2 on the right, k/2 on the left)
    #                                            gamma.estimator.final.mean) # which gamma to use / apply one sided EVT
    # }
    # FRPD.EVT.depth.mean <- apply(TD.EVT.marginal.constant.gamma, 1, mean) # calculate final functional random projection depth with EVT modification with mean
    # FRPD.EVT.depth.min <- apply(TD.EVT.marginal.constant.gamma, 1, min) # calculate final functional random projection depth with EVT modification with min
    # }
    
    
        
  } # end EVT
  
  # function returns under different senarios
  if(more.than.one.depth == 0){
    if(EVT == F & method == "mean"){return(list(depth = FRPD.depth.mean))}
    if(EVT == F & method == "min"){return(list(depth = FRPD.depth.min))}
    if(EVT == T & method == "mean"){return(list(depth = FRPD.EVT.depth.mean.variable.gamma))}
    if(EVT == T & method == "min"){return(list(depth = FRPD.EVT.depth.min.variable.gamma))}    
    
  }else{
    if(more.than.one.depth == 2 & EVT == F){
      return(list(depth.mean = FRPD.depth.mean, depth.min = FRPD.depth.min))
    }else if(more.than.one.depth == 2 & EVT == T){
      return(list(depth.mean = FRPD.depth.mean, 
                  depth.min = FRPD.depth.min, 
                  depth.EVT.mean = FRPD.EVT.depth.mean.variable.gamma, 
                  depth.EVT.min = FRPD.EVT.depth.min.variable.gamma, 
                  final.gamma = gamma.estimator.final))
    }else if(more.than.one.depth == 4 & EVT == T){
      return(list(depth.mean = FRPD.depth.mean, 
                  depth.min = FRPD.depth.min, 
                  depth.EVT.mean = FRPD.EVT.depth.mean.variable.gamma, 
                  depth.EVT.min = FRPD.EVT.depth.min.variable.gamma, 
                  final.gamma = gamma.estimator.final))
    }
  }
  
} # end MFRPD function




