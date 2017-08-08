# 5 most extremes MT COVER

#=============================================================
# ------------- load MT_MAIN_most_extreme_depths_comparaison.R ---------------------
source("./MT_functions/MT_MAIN_most_extreme_depths_comparaison.R")

read.data.again <- F

if(read.data.again == T){
BR.data <- readRDS("./MT_functions/Comparaison_data_simulations/BR_array_1000_100_100.rds")
ET.data <- readRDS("./MT_functions/Comparaison_data_simulations/extremal_t_array.rds_1000_100_1000.rds")
}

print = T

data.to.use.BR <- BR.data[1:1000,,6]


set.seed(16)
col = iwanthue(nrow(data.to.use.BR))[sample(nrow(data.to.use.BR), replace = T)]

if(print == T){pdf("MT_COVER_5_most_extremes.pdf", width = 3.5/2*3, height = 3)}
  par(mar = c(0,0,0,0)+0.1)
  rainbow.plot(data.to.use.BR, col = col, ylab="", xlab="", main="", xaxt="n", yaxt="n", lwd = 2)
if(print == T){dev.off()}


  
  
# set.seed(5)
# depth <- MFRPD(data.to.use.BR)
# depth.order <- order(depth$depth)
# 
# set.seed(16)
# col = iwanthue(nrow(data.to.use.BR))[sample(nrow(data.to.use.BR), replace = T)]
# 
# if(print == T){pdf("MT_COVER_5_most_extremes.pdf", width = 3.5/2*3, height = 3)}
#   par(mar = c(0,0,0,0)+0.1)
#   rainbow.plot(data.to.use.BR[depth.order[1:100],], col = col, ylab="", xlab="", main="", xaxt="n", yaxt="n", lwd = 2)
# if(print == T){dev.off()}
  

