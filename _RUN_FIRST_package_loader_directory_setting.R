# PACKAGE LOADER 


# set appropriate directory

main.directory <- "/Users/gordonschucker/Dropbox/MASTER/Masterarbeit & Paper/R & LaTeX"
setwd(main.directory)

# Remove scientific notation (i.e. 4.1994e-05)
options(scipen=999)

# =============================================================
# ------------------INSTALL PACKAGES----------------------
# =============================================================

# install.packages("devtools")
# install.packages("colorspace")
# install.packages("SAPP")
# install.packages("fda.usc")
# install.packages("SpatialExtremes")
# install.packages("plot3Drgl")
# install.packages("roxgen2")
# install.packages("scatterplot3d")
# install.packages("fextRemes")
# install.packages("rgl")
# 
# Packages for RefinedHD.R
# 
# install.packages("Rlab")
# install.packages("mvtnorm")
# install.packages("depth")
# install.packages("abind")
# install.packages("geometry")
# 
# install.packages("RColorBrewer")
# 
# install.packages("fMultivar")
# 
# install.packages("depthTools")
# install.packages("colorRamps")
# 
# 
# install.packages("viridis")
# install.packages("ddalpha")
# install.packages("xtable")


# =============================================================
# ------------------ LOAD LIBRARIES ----------------------
# =============================================================

# not all are necessary for the main functions though

library(graphics)
library(MASS)
library(colorspace)
library(SAPP)
library(cluster)
library(ldfun)
library(plot3D)
library(fda.usc)
library(scatterplot3d)
library(colorRamps)

library(devtools)
# devtools::install_github("klutometis/roxygen")
# library(roxygen2) # not sure what this was for, but not working in R V3.3.2 anymore...
library(SpatialExtremes) # to simulate max stable processes
#library(rgl)
library(plot3Drgl)

library(Rlab)
library(mvtnorm)
library(MASS)
library(depth)
#library(fExtremes)
library(stats)
library(geometry)

library(RColorBrewer)

library(fMultivar)

library(depthTools)

library(grDevices)
library(viridis)
library(ddalpha)
library(xtable) # to create LaTeX tables
options(xtable.floating = FALSE)
options(xtable.timestamp = "")


# WARNING

# install package from Claudio Agostinelli, half-region depth & co
# install("/Users/gordonschucker/Dropbox/MASTER/Masterarbeit & Paper/R & LaTeX/R codes/ldfun")

# WARNING
# make sure the folder MT_functions is in the main directory

# =============================================================
# ------------------ EVT - EVT functions ----------------------
source("./MT_functions/EVT_functions.R")

# =============================================================
# --------------- Univariate functional depths ----------------------
source("./MT_functions/univariate_functional_depths.R")

# ========================================================================
# --------- Secondary functions ----------------------
source("./MT_functions/secondary_functions.R")







