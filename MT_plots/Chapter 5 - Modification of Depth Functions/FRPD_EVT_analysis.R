#=============================================================
# -------- n.proj projection depth analysis ---------------------


#=============================================================
# ------------- load MT_MAIN_depth_function_for_analysis.R ---------------------
source("./MT_functions/MT_MAIN_depth_function_for_analysis.R")



print = F # if one wants to print pdfs
directory.for.plots <- "//Users//gordonschucker//Dropbox//MASTER//Masterarbeit & Paper//R & LaTeX//n_proj_FPRD_analysis//"


# read finished data
data.array <- readRDS("./MT_functions/Comparaison_data_simulations/main_data_comparaison.rds")
# data.array <- data.array[1:250,,]


n.draws <- dim(data.array)[1]
n.values <- dim(data.array)[2]


col <- viridis_pal(alpha = 1, option = "A", end = 1)(n.draws) # simulate beautiful color palettes


# general correlation plots, last analysis


depth.names <- c("FRPDT1 mean","FRPDT1 min", "FRPDT1 min EVT")

n.proj <- 1000 # number of projections for the FRPD
time.depth <- proc.time()
all.depths <- vector("list",length(depth.names))
for(data.i in 1:2){
all.depths[[data.i]] <- general.proj.EVT.depth(data.array[,,data.i],n.proj = n.proj)  # calculates all the depths given in data.names
}
(proc.time() - time.depth)[3]/60
number.of.depths <- length(depth.names)



# Plot DD plots "FRPDT1 mean" vs "FRPDT1 min" vs "FRPDT1 min EVT"
# BR Data

if(print == T){pdf("DD_plots_FRPD_EVT_BR.pdf", width = 9, height = 5.7)}
data.i = 3

cex = 0.9
par(mfrow = c(2,3), mar = c(4.2, 4, 1, 2) + 0.1, cex = cex)
plot(all.depths[[data.i]][[1]]$depth, all.depths[[data.i]][[3]]$depth, pch = ".", xlab = expression("D"["FRPDT1,mean,1000"]), ylab = expression("D"["FRPDT1,min,EVT,1000"]))
text(0.088,0.09, labels = bquote(rho*" = "*.(round(cor(all.depths[[data.i]][[1]]$depth, all.depths[[data.i]][[3]]$depth),3))))
plot(all.depths[[data.i]][[1]]$depth, all.depths[[data.i]][[2]]$depth, pch = ".", xlab = expression("D"["FRPDT1,mean,1000"]), ylab = expression("D"["FRPDT1,min,1000"]))
text(0.088,0.09, labels = bquote(rho*" = "*.(round(cor(all.depths[[data.i]][[1]]$depth, all.depths[[data.i]][[2]]$depth),3))))
plot(all.depths[[data.i]][[2]]$depth, all.depths[[data.i]][[3]]$depth, pch = ".", xlab = expression("D"["FRPDT1,min,1000"]), ylab = expression("D"["FRPDT1,min,EVT,1000"]))
text(0.027,0.09, labels = bquote(rho*" = "*.(round(cor(all.depths[[data.i]][[2]]$depth, all.depths[[data.i]][[3]]$depth),3))))
plot.lim.x <- c(0,0.2)
plot.lim.y <- c(0,0.027)
plot(all.depths[[data.i]][[1]]$depth, all.depths[[data.i]][[3]]$depth, pch = ".", xlab = expression("D"["FRPDT1,mean,1000"]), ylab = expression("D"["FRPDT1,min,EVT,1000"]), xlim = plot.lim.x, ylim = plot.lim.y)
plot(all.depths[[data.i]][[1]]$depth, all.depths[[data.i]][[2]]$depth, pch = ".", xlab = expression("D"["FRPDT1,mean,1000"]), ylab = expression("D"["FRPDT1,min,1000"]), xlim = plot.lim.x, ylim = plot.lim.y)
plot.lim.x <- c(0,0.04)
plot.lim.y <- c(0,0.04)
plot(all.depths[[data.i]][[2]]$depth, all.depths[[data.i]][[3]]$depth, pch = ".", xlab = expression("D"["FRPDT1,min,1000"]), ylab = expression("D"["FRPDT1,min,EVT,1000"]), xlim = plot.lim.x, ylim = plot.lim.y)

if(print == T){dev.off()}

# same as above with ET data

if(print == T){pdf("DD_plots_FRPD_EVT_ET.pdf", width = 9, height = 5.7)}
data.i = 4

cex = 0.9
par(mfrow = c(2,3), mar = c(4.2, 4, 1, 2) + 0.1, cex = cex)
plot(all.depths[[data.i]][[1]]$depth, all.depths[[data.i]][[3]]$depth, pch = ".", xlab = expression("D"["FRPDT1,mean,1000"]), ylab = expression("D"["FRPDT1,min,EVT,1000"]))
text(0.086,0.115, labels = bquote(rho*" = "*.(round(cor(all.depths[[data.i]][[1]]$depth, all.depths[[data.i]][[3]]$depth),3))))
plot(all.depths[[data.i]][[1]]$depth, all.depths[[data.i]][[2]]$depth, pch = ".", xlab = expression("D"["FRPDT1,mean,1000"]), ylab = expression("D"["FRPDT1,min,1000"]))
text(0.086,0.115, labels = bquote(rho*" = "*.(round(cor(all.depths[[data.i]][[1]]$depth, all.depths[[data.i]][[2]]$depth),3))))
plot(all.depths[[data.i]][[2]]$depth, all.depths[[data.i]][[3]]$depth, pch = ".", xlab = expression("D"["FRPDT1,min,1000"]), ylab = expression("D"["FRPDT1,min,EVT,1000"]))
text(0.032,0.115, labels = bquote(rho*" = "*.(round(cor(all.depths[[data.i]][[2]]$depth, all.depths[[data.i]][[3]]$depth),3))))
plot.lim.x <- c(0,0.2)
plot.lim.y <- c(0,0.027)
plot(all.depths[[data.i]][[1]]$depth, all.depths[[data.i]][[3]]$depth, pch = ".", xlab = expression("D"["FRPDT1,mean,1000"]), ylab = expression("D"["FRPDT1,min,EVT,1000"]), xlim = plot.lim.x, ylim = plot.lim.y)
plot(all.depths[[data.i]][[1]]$depth, all.depths[[data.i]][[2]]$depth, pch = ".", xlab = expression("D"["FRPDT1,mean,1000"]), ylab = expression("D"["FRPDT1,min,1000"]), xlim = plot.lim.x, ylim = plot.lim.y)
plot.lim.x <- c(0,0.04)
plot.lim.y <- c(0,0.04)
plot(all.depths[[data.i]][[2]]$depth, all.depths[[data.i]][[3]]$depth, pch = ".", xlab = expression("D"["FRPDT1,min,1000"]), ylab = expression("D"["FRPDT1,min,EVT,1000"]), xlim = plot.lim.x, ylim = plot.lim.y)

if(print == T){dev.off()}

















if(print == T){pdf("cor_most_extremes.pdf", width = 9, height = 9)}
par(mfrow= c(2,2), mar = c(4.2, 4, 1, 2) + 0.1, oma = c(0,0,0,0))
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
     xlab = "k most extreme depth values", ylab = "correlation", main = "ordered by FRPDT1 mean, Brown-Resnick data") 



cor.most.extremes <- rep(NA, length(most.extremes.max.seq))
depth.j <- 3
order.depth.j <- all.depths[[data.i]][[depth.j]]$order
depth.k <- 1
for(most.extremes.max in most.extremes.max.seq){
cor.most.extremes[most.extremes.max] <- cor(all.depths[[data.i]][[depth.j]]$depth[order.depth.j][1:most.extremes.max],all.depths[[data.i]][[depth.k]]$depth[order.depth.j][1:most.extremes.max])
}

plot(most.extremes.max.seq[10:n.draws],cor.most.extremes[10:n.draws], type = "l",
     xlab = "k most extreme depth values", ylab = "correlation", main = "ordered by FRPDT1 min EVT, Brown-Resnick data") 


data.i <- 4
cor.most.extremes <- rep(NA, length(most.extremes.max.seq))
depth.j <- 1
order.depth.j <- all.depths[[data.i]][[depth.j]]$order
depth.k <- 3
for(most.extremes.max in most.extremes.max.seq){
cor.most.extremes[most.extremes.max] <- cor(all.depths[[data.i]][[depth.j]]$depth[order.depth.j][1:most.extremes.max],all.depths[[data.i]][[depth.k]]$depth[order.depth.j][1:most.extremes.max])
}

plot(most.extremes.max.seq[10:n.draws],cor.most.extremes[10:n.draws], type = "l",
     xlab = "k most extreme depth values", ylab = "correlation", main = "ordered by FRPDT1 mean, extremal-t data") 



cor.most.extremes <- rep(NA, length(most.extremes.max.seq))
depth.j <- 3
order.depth.j <- all.depths[[data.i]][[depth.j]]$order
depth.k <- 1
for(most.extremes.max in most.extremes.max.seq){
cor.most.extremes[most.extremes.max] <- cor(all.depths[[data.i]][[depth.j]]$depth[order.depth.j][1:most.extremes.max],all.depths[[data.i]][[depth.k]]$depth[order.depth.j][1:most.extremes.max])
}

plot(most.extremes.max.seq[10:n.draws],cor.most.extremes[10:n.draws], type = "l",
     xlab = "k most extreme depth values", ylab = "correlation", main = "ordered by FRPDT1 min EVT, extremal-t data") 

if(print == T){dev.off()}


# BONUS GRAPHIC


if(print == T){pdf("gamma_FRPD.pdf", width = 10, height = 4.7)}

par(mfrow = c(1,2))
plot(density(all.depths[[3]][[1]]$gamma.final), xlab = "gamma", ylab = "density", main = "Brown-Resnick data")
abline(v = mean(all.depths[[3]][[1]]$gamma.final), col = "red", lty = 2)
plot(density(all.depths[[4]][[1]]$gamma.final), xlab = "gamma", ylab = "density", main = "extremal-t data")
abline(v = mean(all.depths[[4]][[1]]$gamma.final), col = "red", lty = 2)
if(print == T){dev.off()}









