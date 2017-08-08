# Extremal depth simulations



print <- T



x <- y <- rbind(c(2,2.1,1.9,1.3,0.8,-0.4),c(1.5,1.5,1.5,1.5,1.5,1.5),c(1,1,1.1,1,1,1.1),c(0.9,0.8,0.9,1.6,2.1,3.5),c(0.5,0.4,-0.1,-0.8,-1.6,-2.5),c(0,0.3,0,0.2,0,-0.1),c(-0.8,-0.7,-0.3,-0.4,-0.8,-0.8),c(-1.5,-1.2,-1.2,-1.4,-1.4,-1.5))


n1.x <- dim(x)[1]
n2.x <- dim(x)[2]
n1.y <- dim(y)[1]
n2.y <- dim(y)[2]


x.depth <- MDF(x, ED = F)

order.depth <- order(x.depth$function.depth.x)



if(print == T){pdf("extremal_depth_8_curves.pdf", height = 5, width = 9)}

par(mfrow = c(1,2))

position.text.left <- rep(1.5,8)
position.text.left2 <- rep(5.2,8)

for(i in 1:8){
  if(i == 1){
#text(position.text.left[i],x[i,position.text.left[i]]+0.2, label = expression("x"), cex = 1)
        text(position.text.left2[i],x[i,position.text.left2[i]]-0.1, label = expression("x"), cex = 1)
#text(position.text.left[i]+0.1,x[i,position.text.left[i]]-0.1+0.2, label = i, cex = 0.7) 
text(position.text.left2[i]+0.1,x[i,position.text.left2[i]]-0.1-0.1, label = i, cex = 0.7)  
    }else{
#    text(position.text.left[i],x[i,position.text.left[i]], label = expression("x"), cex = 1)
#text(position.text.left[i]+0.1,x[i,position.text.left[i]]-0.1, label = i, cex = 0.7)  

      if(i == 5){
        text(position.text.left2[i]+0.1,x[i,position.text.left2[i]]-3.3, label = i, cex = 0.7)
        text(position.text.left2[i],x[i,position.text.left2[i]], label = expression("x"), cex = 1)
      }else{
        text(position.text.left2[i]+0.1,x[i,position.text.left2[i]]-0.1, label = i, cex = 0.7)
        text(position.text.left2[i],x[i,position.text.left2[i]], label = expression("x"), cex = 1)
        }    

  }
}

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
  



# fancy test Phi.x.r





#col <- viridis_pal(alpha = 1, option = "A", end = 0.9)(dim(x)[1])

col1 <- viridis_pal(alpha = 1, option = "A", end = 0.9, begin = 0.4)(4)
col2 <- viridis_pal(alpha = 1, option = "D", end = 0.9, begin = 0.4)(4)

col <- c(col1[1],col2[1], col1[2],col2[2], col1[3],col2[3], col1[4],col2[4])

set.seed(10)

#col <- iwanthue(8, random = T)

phi.x.r.order.for.plot <- jitter(Phi.x.r[order.depth,], factor = 0.1)


rainbow.plot(phi.x.r.order.for.plot, normalize(0:8),  col = col, ylab = expression(Phi["x"["i"]*",n"]*"(q)"), xlab = "q", lwd = 2)

text(0.135,0.822, label = expression(Phi["x"["8"]*",n"]), cex = 1.00, col = "black")
arrows(0.17,0.78,0.25,0.72, length = 0.08, col = "black")

text(0.95,0.467, label = expression(Phi["x"["6"]*",n"]), cex = 1.00, col = "black")
arrows(0.88,0.45,0.8,0.4, length = 0.08, col = "black")

# position.text.left <- rep(1.5,8)
# position.text.left2 <- rep(5.2,8)
# 
# for(i in 1:8){
#   if(i == 0){
# text(position.text.left[i],Phi.x.r[i,position.text.left[i]]+0.2, label = expression("x"), cex = 1)
#         text(position.text.left2[i],Phi.x.r[i,position.text.left2[i]]-0.1, label = expression("x"), cex = 1)
# text(position.text.left[i]+0.1,Phi.x.r[i,position.text.left[i]]-0.1+0.2, label = i, cex = 0.7)
# text(position.text.left2[i]+0.1,Phi.x.r[i,position.text.left2[i]]-0.1-0.1, label = i, cex = 0.7)  
#     }else{
#    text(position.text.left[i],Phi.x.r[i,position.text.left[i]], label = expression("x"), cex = 1)
# text(position.text.left[i]+0.0,Phi.x.r[i,position.text.left[i]]-0.0, label = i, cex = 0.7)
# 
# text(position.text.left2[i],Phi.x.r[i,position.text.left2[i]], label = expression("x"), cex = 1)
# text(position.text.left2[i]+0.0,Phi.x.r[i,position.text.left2[i]]-0.0, label = i, cex = 0.7)  
#   }
# }







rainbow.plot(x[order.depth,], col = col, ylab = expression("x"["i"]*"(t)"), lwd = 2)
for(i in 1:8){
  if(i == 1){
#text(position.text.left[i],x[i,position.text.left[i]]+0.2, label = expression("x"), cex = 1)
        text(position.text.left2[i],x[i,position.text.left2[i]]-0.1, label = expression("x"), cex = 1)
#text(position.text.left[i]+0.1,x[i,position.text.left[i]]-0.1+0.2, label = i, cex = 0.7) 
text(position.text.left2[i]+0.1,x[i,position.text.left2[i]]-0.1-0.1, label = i, cex = 0.7)  
    }else{
#    text(position.text.left[i],x[i,position.text.left[i]], label = expression("x"), cex = 1)
#text(position.text.left[i]+0.1,x[i,position.text.left[i]]-0.1, label = i, cex = 0.7)  

      if(i == 5){
        text(position.text.left2[i]+0.1,x[i,position.text.left2[i]]-0.5, label = i, cex = 0.7)
        text(position.text.left2[i],x[i,position.text.left2[i]]-0.4, label = expression("x"), cex = 1)
      }else{
        text(position.text.left2[i]+0.1,x[i,position.text.left2[i]]-0.1, label = i, cex = 0.7)
        text(position.text.left2[i],x[i,position.text.left2[i]], label = expression("x"), cex = 1)
        }    
 
  }
}


if(print == T){dev.off()}

Phi.x.r[,c(1,2,4,6,8)]





