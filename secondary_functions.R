# SECONDARY FUNCTIONS

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# to normalize data, in particular depth (should also exist as R function?)

normalize <- function(x) { 
  x <- (x - min(x))/(max(x) - min(x))
  return(x)
}

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# transforms a matrix into polygon for the function polygon (used for plotting)

polygon.fct <- function(x.star, matrix.for.polygon){
  polygon.l <- apply(matrix.for.polygon, 2, min)
  polygon.u <- apply(matrix.for.polygon, 2, max)
 
  result.x <- c(x.star, rev(x.star))
  result.y <- c(polygon.l, rev(polygon.u))
  
 # results <- rbind(result.x,result.y)
  results <- list(x = result.x, y = result.y)
 return(results)  
}

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# generate nice colors

swatch <- function(x) {
  # x: a vector of colours (hex, numeric, or string)
  par(mai=c(0.2, max(strwidth(x, "inch") + 0.4, na.rm = TRUE), 0.2, 0.4))
  barplot(rep(1, length(x)), col=rev(x), space = 0.1, axes=FALSE, 
          names.arg=rev(x), cex.names=0.8, horiz=T, las=1)  
}

# Example:
# swatch(colours()[1:40])
# swatch(iwanthue(5))
# swatch(1:4)

iwanthue <- function(n, hmin=0, hmax=360, cmin=0, cmax=180, lmin=0, lmax=100, 
                     plot=FALSE, random=FALSE) {
  # Presently doesn't allow hmax > hmin (H is circular)
  # n: number of colours
  # hmin: lower bound of hue (0-360)
  # hmax: upper bound of hue (0-360)
  # cmin: lower bound of chroma (0-180)
  # cmax: upper bound of chroma (0-180)
  # lmin: lower bound of luminance (0-100)
  # lmax: upper bound of luminance (0-100)
  # plot: plot a colour swatch?
  # random: should clustering be random? (if FALSE, seed will be set to 1,
  #         and the RNG state will be restored on exit.) 
  require(colorspace)
  stopifnot(hmin >= 0, cmin >= 0, lmin >= 0, 
            hmax <= 360, cmax <= 180, lmax <= 100, 
            hmin <= hmax, cmin <= cmax, lmin <= lmax,
            n > 0)
  if(!random) {
    if (exists(".Random.seed", .GlobalEnv)) {
      old_seed <- .GlobalEnv$.Random.seed
      on.exit(.GlobalEnv$.Random.seed <- old_seed)
    } else {
      on.exit(rm(".Random.seed", envir = .GlobalEnv))
    }
    set.seed(1)
  }
  lab <- LAB(as.matrix(expand.grid(seq(0, 100, 1), 
                                   seq(-100, 100, 5), 
                                   seq(-110, 100, 5))))
  if (any((hmin != 0 || cmin != 0 || lmin != 0 ||
           hmax != 360 || cmax != 180 || lmax != 100))) {
    hcl <- as(lab, 'polarLUV')
    hcl_coords <- coords(hcl)
    hcl <- hcl[which(hcl_coords[, 'H'] <= hmax & hcl_coords[, 'H'] >= hmin &
                       hcl_coords[, 'C'] <= cmax & hcl_coords[, 'C'] >= cmin & 
                       hcl_coords[, 'L'] <= lmax & hcl_coords[, 'L'] >= lmin), ]
    #hcl <- hcl[-which(is.na(coords(hcl)[, 2]))]
    lab <- as(hcl, 'LAB')    
  }
  lab <- lab[which(!is.na(hex(lab))), ]
  clus <- kmeans(coords(lab), n, iter.max=50)
  if (isTRUE(plot)) {
    swatch(hex(LAB(clus$centers)))
  }
  hex(LAB(clus$centers))
}






#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

# Rainbow plot for univaraite functional data
# 
# data = matrix with functional observations displayed by row 
# col = wanted color (needs to be nrow(data) long) # recommended use with $col from depth function
# T are the value for each evaluation point

rainbow.plot <- function(data, x.star = 1:ncol(data), col, ylab = "X(t)", xlab = "t", main = "", ylim = NULL, ann = T, xaxt='s',yaxt='s', lwd = 1){
  
  if(missing(col)){
    col <- iwanthue(nrow(data)) # random colors 
  }
  
  
  values <- t(data) # data needs to be transposed for matplot
  values <- cbind(x = x.star, as.data.frame(values))
  matplot(x=values[,1], y=values[,-1], lty= 1, type = "l", col = col, ylab = ylab, xlab = xlab, main = main, ylim = ylim, ann = ann, xaxt = xaxt, yaxt = yaxt, lwd = lwd)
}

