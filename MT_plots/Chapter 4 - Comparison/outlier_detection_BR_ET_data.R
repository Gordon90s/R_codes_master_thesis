# MT_MAIN_outlier_detection

#=============================================================
# ------------- load MT_MAIN_most_extreme_depths_comparaison.R ---------------------
source("./MT_functions/MT_MAIN_depth_function_for_analysis.R")

print = T

BR.data <- readRDS("./MT_functions/Comparaison_data_simulations/BR_array_1000_100_100.rds")
ET.data <- readRDS("./MT_functions/Comparaison_data_simulations/extremal_t_array.rds_1000_100_1000.rds")


n <- 5
simulation.names <- c("simulation 1", "simulation 2", "simulation 3", "simulation 4", "simulation 5")
most.extremes <- 5

data.to.use.BR <- BR.data[1:1000,,1:n]
data.to.use.ET <- ET.data[1:1000,,1:n]

set.seed(5)



start.timing <- proc.time()

most.extremes.depth.BR <- vector("list", n) 
most.extremes.depth.ET <- vector("list", n)

for(i in 1:n){
most.extremes.depth.BR[[i]] <- most.extremes.depth(data.to.use.BR[,,i],most.extremes)
most.extremes.depth.ET[[i]] <- most.extremes.depth(data.to.use.ET[,,i],most.extremes)
}

proc.time() - start.timing


number.of.depths <- 4
depth.names <- c("IDT1 alpha","FRPDT1 mean","FRPDT1 min","ED")



cex = 1

if(print == T){pdf("most_extremes_BR.pdf", height = 10, width = 10)}

#par(mfrow = c(5,n), mar = c(5, 4, 4, 2) + 0.1)
par(mfrow = c(5,n), oma = c(0,2.2,2.2,0))

for(j in 1:n){
  par(mar = c(0,0,0,0))

rainbow.plot(data.to.use.BR[,,j], col = most.extremes.depth.BR[[j]]$col, ylim = range(data.to.use.BR[,,j]), ylab="", xlab="", main="", xaxt="n", yaxt="n")
  mtext(text=simulation.names[j],side=3,line=1, cex = cex)
  if(j == 1) mtext(text="all data",side=2,line=1, cex = cex)
}
for(i in 1:number.of.depths){
  for(j in 1:n){
  rainbow.plot(most.extremes.depth.BR[[j]]$sorted.data.outlier[,,i], col = most.extremes.depth.BR[[j]]$col.outlier[i,], ylim = range(data.to.use.BR[,,j]), ylab="", xlab="", main="", xaxt="n", yaxt="n")  
  if(j == 1) mtext(text=depth.names[i],side=2,line=1, cex = cex)
    }
}

if(print == T){dev.off()}


if(print == T){pdf("most_extremes_ET.pdf", height = 10, width = 10)}

#par(mfrow = c(5,n), mar = c(5, 4, 4, 2) + 0.1)
par(mfrow = c(5,n), oma = c(0,2.2,2.2,0))

for(j in 1:n){
  par(mar = c(0,0,0,0))

rainbow.plot(data.to.use.ET[,,j], col = most.extremes.depth.ET[[j]]$col, ylim = range(data.to.use.ET[,,j]), ylab="", xlab="", main="", xaxt="n", yaxt="n")
  mtext(text=simulation.names[j],side=3,line=1, cex = cex)
  if(j == 1) mtext(text="all data",side=2,line=1, cex = cex)
}
for(i in 1:number.of.depths){
  for(j in 1:n){
  rainbow.plot(most.extremes.depth.ET[[j]]$sorted.data.outlier[,,i], col = most.extremes.depth.ET[[j]]$col.outlier[i,], ylim = range(data.to.use.ET[,,j]), ylab="", xlab="", main="", xaxt="n", yaxt="n")  
  if(j == 1) mtext(text=depth.names[i],side=2,line=1, cex = cex)
    }
}

if(print == T){dev.off()}




