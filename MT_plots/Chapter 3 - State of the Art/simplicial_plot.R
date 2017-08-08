# simplicial plot using jitter technique

library(mvtnorm)
library(scales)
library(fMultivar)
library(base)


print = F

n.x <- 60

#set.seed(18)
set.seed(22)

x.norm <- rmvnorm(n.x, mean = c(0,0)) 
#plot(x.norm)

#+ rcauchy2d(n.x, rho = 0)/4
jitter.factor <- 0.05





                  
if(print == T){pdf("simplicial_depth.pdf", width = 9, height = 9)}
par(mar = c(6,5,5,2)+0.1)
par(mfrow = c(1,1), cex = 1.3)

  plot(x.norm, pch = ".", xlab = expression("x"[1]), ylab = expression("x"[2])) 
  for(i in 1:n.x){
    for(j in 1:n.x){
      for(k in 1:n.x){
        if(i < j & j < k){
        triangle <- cbind(jitter(x.norm[i,], factor = jitter.factor),jitter(x.norm[j,], factor = jitter.factor),jitter(x.norm[k,], factor = jitter.factor))
        polygon(triangle[1,],triangle[2,], border = alpha("black", 0.005))
        print(cbind(i,j,k))
        }
      }
    }
  }

if(print == T){dev.off()}
