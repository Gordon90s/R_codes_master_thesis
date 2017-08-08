#=============================================================
# -------- correlation plots for depth analysis ---------------------


#=============================================================
# ------------- load MT_MAIN_depth_function_for_analysis.R ---------------------
source("./MT_functions/MT_MAIN_depth_function_for_analysis.R")



print = F # if one wants to print pdfs
directory.for.plots <- "//Users//gordonschucker//Dropbox//MASTER//Masterarbeit & Paper//R & LaTeX//n_proj_FPRD_analysis//"


# read finished data
data.array <- readRDS("./MT_functions/Comparaison_data_simulations/main_data_comparaison.rds")
data.names <- c("Gaussian", "Gaussian with amplitude", "Brown-Resnick","Extremal-t")
# data.array <- data.array[1:200,,]


n.draws <- dim(data.array)[1]
n.values <- dim(data.array)[2]

col <- viridis_pal(alpha = 1, option = "A", end = 1)(n.draws) # simulate beautiful color palettes

#=============================================================
# -------- INTEGRATED DEPTH ANALYSIS ---------------------

depth.names <- c("IDT1","IDS1","IDT1 alpha","IDS1 alpha","MHRD")


time.depth <- proc.time()
all.depths <- vector("list",4)
for(data.i in 1:4){
all.depths[[data.i]] <- general.int.depth(data.array[,,data.i])  # calculates all the depths given in data.names
}
(proc.time() - time.depth)[3]/60
number.of.depths <- length(depth.names)
#################################

directory.for.plots <- "//Users//gordonschucker//Dropbox//MASTER//Masterarbeit & Paper//R & LaTeX//"


for(data.i in 1:4){
  if(print == T){pdf(paste0(directory.for.plots,"comparison_int_new_data",data.i,".pdf"), width = 10, height = 10)}
  
  n <- number.of.depths
  par(mfrow = c(n,n), mar = c(0,0,0,0) + 0.0)
  for(depth.j in 1:n){
    for(depth.k in 1:n){
      if(depth.j>depth.k){
         plot(normalize(all.depths[[data.i]][[depth.j]]$depth),normalize(all.depths[[data.i]][[depth.k]]$depth), pch = ".",
           xaxt='n',yaxt='n', ann=FALSE, xlim = c(-0.1,1.1), ylim = c(-0.1,1.1))
      }else if(depth.j == depth.k){
        
        if(data.i < 3){
          rainbow.plot(all.depths[[data.i]][[depth.j]]$data.sorted, col = col,
                     xaxt = 'n', yaxt = 'n', ann=F,
               ylim = c(min(all.depths[[data.i]][[depth.j]]$data.sorted),
                          max(all.depths[[data.i]][[depth.j]]$data.sorted)+diff(range(all.depths[[data.i]][[depth.j]]$data.sorted))*1/5))
        } else{
          rainbow.plot(log(all.depths[[data.i]][[depth.j]]$data.sorted), col = col,
                     xaxt = 'n', yaxt = 'n', ann=F,
               ylim = c(min(log(all.depths[[data.i]][[depth.j]]$data.sorted)),
                          max(log(all.depths[[data.i]][[depth.j]]$data.sorted))+diff(range(log(all.depths[[data.i]][[depth.j]]$data.sorted)))*1/5))
          
        }
        
        if(data.i < 3){
        text(x = n.values/2, y = max(all.depths[[data.i]][[depth.j]]$data.sorted)+diff(range(all.depths[[data.i]][[depth.j]]$data.sorted))*1/10, 
             labels = bquote(.(depth.names[depth.j])), 
             cex = 1.3, col = "black")
        }else{
        text(x = n.values/2, y = max(log(all.depths[[data.i]][[depth.j]]$data.sorted))+diff(log(range(all.depths[[data.i]][[depth.j]]$data.sorted)))*1/10, 
             labels = bquote(.(depth.names[depth.j])), 
             cex = 1.3, col = "black")  
        }
      }else{
              plot(c(0, 1), c(0, 1), ann = T, type = 'n', xaxt = 'n', yaxt = 'n', xlim = c(-0.1,1.1), ylim = c(-0.1,1.1))
        text(x = 0.5, y = 0.5, labels = bquote(rho == .(round(cor(normalize(all.depths[[data.i]][[depth.j]]$depth),normalize(all.depths[[data.i]][[depth.k]]$depth)),3))), 
       cex = 1.2, col = "black")
      }
    }
  }
  if(print == T){dev.off()}
}#end for(data.i in 1:4)

#=============================================================
# -------- GENERAL DEPTH ANALYSIS ---------------------

depth.names <- c("IDT1","IDT1 alpha","FRPDT1 mean","FRPDT1 min","ED")

n.proj <- 1000
time.depth <- proc.time()
all.depths <- vector("list",4)
for(data.i in 1:4){
all.depths[[data.i]] <- general.depth.analysis(data.array[,,data.i], n.proj = n.proj, proj.type = "normal")  # calculates all the depths given in data.names
}
(proc.time() - time.depth)[3]/60
number.of.depths <- length(depth.names)
#################################

directory.for.plots <- "//Users//gordonschucker//Dropbox//MASTER//Masterarbeit & Paper//R & LaTeX//"


for(data.i in 1:4){
  if(print == T){pdf(paste0(directory.for.plots,"comparison_new_data",data.i,".pdf"), width = 10, height = 10)}
  
  n <- number.of.depths
  par(mfrow = c(n,n), mar = c(0,0,0,0) + 0.0)
  for(depth.j in 1:n){
    for(depth.k in 1:n){
      if(depth.j>depth.k){
         plot(normalize(all.depths[[data.i]][[depth.j]]$depth),normalize(all.depths[[data.i]][[depth.k]]$depth), pch = ".",
           xaxt='n',yaxt='n', ann=FALSE, xlim = c(-0.1,1.1), ylim = c(-0.1,1.1))
      }else if(depth.j == depth.k){
        
        if(data.i < 3){
          rainbow.plot(all.depths[[data.i]][[depth.j]]$data.sorted, col = col,
                     xaxt = 'n', yaxt = 'n', ann=F,
               ylim = c(min(all.depths[[data.i]][[depth.j]]$data.sorted),
                          max(all.depths[[data.i]][[depth.j]]$data.sorted)+diff(range(all.depths[[data.i]][[depth.j]]$data.sorted))*1/5))
        } else{
          rainbow.plot(log(all.depths[[data.i]][[depth.j]]$data.sorted), col = col,
                     xaxt = 'n', yaxt = 'n', ann=F,
               ylim = c(min(log(all.depths[[data.i]][[depth.j]]$data.sorted)),
                          max(log(all.depths[[data.i]][[depth.j]]$data.sorted))+diff(range(log(all.depths[[data.i]][[depth.j]]$data.sorted)))*1/5))
          
        }
        
        if(data.i < 3){
        text(x = n.values/2, y = max(all.depths[[data.i]][[depth.j]]$data.sorted)+diff(range(all.depths[[data.i]][[depth.j]]$data.sorted))*1/10, 
             labels = bquote(.(depth.names[depth.j])), 
             cex = 1.3, col = "black")
        }else{
        text(x = n.values/2, y = max(log(all.depths[[data.i]][[depth.j]]$data.sorted))+diff(log(range(all.depths[[data.i]][[depth.j]]$data.sorted)))*1/10, 
             labels = bquote(.(depth.names[depth.j])), 
             cex = 1.3, col = "black")  
        }
      }else{
              plot(c(0, 1), c(0, 1), ann = T, type = 'n', xaxt = 'n', yaxt = 'n', xlim = c(-0.1,1.1), ylim = c(-0.1,1.1))
        text(x = 0.5, y = 0.5, labels = bquote(rho == .(round(cor(normalize(all.depths[[data.i]][[depth.j]]$depth),normalize(all.depths[[data.i]][[depth.k]]$depth)),3))), 
       cex = 1.2, col = "black")
      }
    }
  }
  if(print == T){dev.off()}
}#end for(data.i in 1:4)



#=============================================================
# -------- FUNCTIONAL EVT DEPTH ANALYSIS ---------------------


depth.names <- c("FRPDT1 mean","FRPDT1 min","FRPDT1 min EVT")

n.proj <- 1000 # number of projections for the FRPD
time.depth <- proc.time()
all.depths <- vector("list",4)
for(data.i in 3:4){
all.depths[[data.i]] <- general.proj.EVT.FRV.depth(data.array[,,data.i],n.proj = n.proj)  # calculates all the depths given in data.names
}
(proc.time() - time.depth)[3]/60
number.of.depths <- length(depth.names)




for(data.i in 3:4){
  if(print == T){pdf(paste0(directory.for.plots,"EVTFRPD_new_data",data.i,".pdf"), width = 10, height = 10)}
  
  n <- number.of.depths
  par(mfrow = c(n,n), mar = c(0,0,0,0) + 0.0)
  for(depth.j in 1:n){
    for(depth.k in 1:n){
      if(depth.j>depth.k){
         plot(normalize(all.depths[[data.i]][[depth.j]]$depth),normalize(all.depths[[data.i]][[depth.k]]$depth), pch = ".",
           xaxt='n',yaxt='n', ann=FALSE, xlim = c(-0.1,1.1), ylim = c(-0.1,1.1))
      }else if(depth.j == depth.k){
        
        if(data.i < 3){
          rainbow.plot(all.depths[[data.i]][[depth.j]]$data.sorted, col = col,
                     xaxt = 'n', yaxt = 'n', ann=F,
               ylim = c(min(all.depths[[data.i]][[depth.j]]$data.sorted),
                          max(all.depths[[data.i]][[depth.j]]$data.sorted)+diff(range(all.depths[[data.i]][[depth.j]]$data.sorted))*1/5))
        } else{
          rainbow.plot(log(all.depths[[data.i]][[depth.j]]$data.sorted), col = col,
                     xaxt = 'n', yaxt = 'n', ann=F,
               ylim = c(min(log(all.depths[[data.i]][[depth.j]]$data.sorted)),
                          max(log(all.depths[[data.i]][[depth.j]]$data.sorted))+diff(range(log(all.depths[[data.i]][[depth.j]]$data.sorted)))*1/5))
          
        }
        
        if(data.i < 3){
        text(x = n.values/2, y = max(all.depths[[data.i]][[depth.j]]$data.sorted)+diff(range(all.depths[[data.i]][[depth.j]]$data.sorted))*1/10, 
             labels = bquote(.(depth.names[depth.j])), 
             cex = 1.3, col = "black")
        }else{
        text(x = n.values/2, y = max(log(all.depths[[data.i]][[depth.j]]$data.sorted))+diff(log(range(all.depths[[data.i]][[depth.j]]$data.sorted)))*1/10, 
             labels = bquote(.(depth.names[depth.j])), 
             cex = 1.3, col = "black")  
        }
      }else{
              plot(c(0, 1), c(0, 1), ann = T, type = 'n', xaxt = 'n', yaxt = 'n', xlim = c(-0.1,1.1), ylim = c(-0.1,1.1))
        text(x = 0.5, y = 0.5, labels = bquote(rho == .(round(cor(normalize(all.depths[[data.i]][[depth.j]]$depth),normalize(all.depths[[data.i]][[depth.k]]$depth)),3))), 
       cex = 1.2, col = "black")
      }
    }
  }
  if(print == T){dev.off()}
}#end for(data.i in 1:4)



data.i <- 3
depth.j <- 1 # depth with which to order and used as "base" order
most.extremes <- 1:20
order.depth.j <- all.depths[[data.i]][[depth.j]]$order
extrem.depths1 <- all.depths[[data.i]][[depth.j]]$depth[order.depth.j][most.extremes]
depth.j <- 3
extrem.depths2 <- all.depths[[data.i]][[depth.j]]$depth[order.depth.j][most.extremes] 

par(mfrow = c(1,1), mar = c(4.2, 4, 1, 2) + 0.1, cex = 1)
plot(extrem.depths1,extrem.depths2, pch = ".", ylab = "EVTFRPD,min", xlab = "FRPD,mean")
cor(extrem.depths1,extrem.depths2)
cbind(extrem.depths1,extrem.depths2)

rainbow.plot(all.depths[[data.i]][[depth.j]]$data.sorted, col = "black")
lines(all.depths[[data.i]][[depth.j]]$data.sorted[1,], col = "red")
lines(all.depths[[data.i]][[depth.j]]$data.sorted[2,], col = "red")
lines(all.depths[[data.i]][[depth.j]]$data.sorted[3,], col = "red")
lines(all.depths[[data.i]][[depth.j]]$data.sorted[4,], col = "red")
lines(all.depths[[data.i]][[depth.j]]$data.sorted[5,], col = "red")

if(print == T){pdf("cor_most_extremes_new.pdf", width = 10, height = 6)}
par(mfrow= c(1,2), mar = c(4.2, 4, 1, 2) + 0.1)
data.i <- 3
most.extremes.max.seq <- 1:n.draws
cor.most.extremes <- rep(NA, length(most.extremes.max.seq))
depth.j <- 1
order.depth.j <- all.depths[[data.i]][[depth.j]]$order
depth.k <- 3
for(most.extremes.max in most.extremes.max.seq){
cor.most.extremes[most.extremes.max] <- cor(all.depths[[data.i]][[depth.j]]$depth[order.depth.j][1:most.extremes.max],all.depths[[data.i]][[depth.k]]$depth[order.depth.j][1:most.extremes.max])
}


plot(most.extremes.max.seq[10:n.draws],cor.most.extremes[10:n.draws], type = "l",
     xlab = "k most extreme depths", ylab = "correlation") 




# copy paste for data.i <- 4



data.i <- 4
most.extremes.max.seq <- 1:n.draws
cor.most.extremes <- rep(NA, length(most.extremes.max.seq))
depth.j <- 1
order.depth.j <- all.depths[[data.i]][[depth.j]]$order
depth.k <- 3
for(most.extremes.max in most.extremes.max.seq){
cor.most.extremes[most.extremes.max] <- cor(all.depths[[data.i]][[depth.j]]$depth[order.depth.j][1:most.extremes.max],all.depths[[data.i]][[depth.k]]$depth[order.depth.j][1:most.extremes.max])
}


plot(most.extremes.max.seq[10:n.draws],cor.most.extremes[10:n.draws], type = "l",
     xlab = "k most extreme depths", ylab = "correlation") 













if(print == T){dev.off()}


