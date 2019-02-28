
# how does


library(tidyverse)
library(raster)
library(scales)
library(colormap)
library(grid)
library(gridExtra)
library(ecoclim)
library(mgcv)

setwd("E:/edges/edgy")

source("code/agg.r") # modified version of raster::aggregate

select <- dplyr::select



##### data setup #####

# load climate data

clim <- stack(c("F:/chelsa/bio19/CHELSA_bio10_1.tif",
                "F:/chelsa/cmip5_bio19/CHELSA_bio_mon_CCSM4_rcp60_r1i1p1_g025.nc_1_2041-2060.tif",
                "F:/chelsa/bio19/CHELSA_bio10_12.tif",
                "F:/chelsa/cmip5_bio19/CHELSA_bio_mon_CCSM4_rcp60_r1i1p1_g025.nc_12_2041-2060.tif"))
names(clim) <- c("temp0", "temp1", "ppt0", "ppt1")

# exclude polar regions, where square landscapes are hardly square
clim <- clim %>%
      crop(extent(-180, 180, -66, 66)) %>%
      reclassify(c(-Inf, -2000, NA))

clim[[3]] <- log10(clim[[3]])
clim[[4]] <- log10(clim[[4]])

writeRaster(clim, "../climate_change.tif")
clim <- stack("../climate_change.tif")
names(clim) <- c("temp0", "temp1", "ppt0", "ppt1")

means <- aggregate(clim, 50)

# 1 if temp and precip are both increasing or decreasing, 0 if signs differ
delta_alignment <- function(x, ...){
      dt <- x[2] - x[1]
      dp <- x[4] - x[3]
      as.integer(sign(dt) == sign(dp))
}
deltas <- calc(means, delta_alignment)

# 1 if positive correlation, 0 if negative
correlation <- function(x, ...){
      m <- matrix(x, ncol=2)
      r <- cor(m[,1], m[,2], use="pairwise.complete.obs")
      as.integer(sign(r) == 1)
}
space <- aggregate(clim[[c(1, 3)]], fact=c(50, 50, 2), fun=correlation)

# 1 if match, 0 if mismatch
m <- stack(space, deltas) %>%
      calc(function(x, ...){x[1] == x[2]}) %>%
      rasterToPoints() %>%
      as.data.frame()

p <- ggplot(m, aes(x, y, fill=factor(layer))) +
      geom_raster() +
      scale_fill_manual(values=c("darkgreen", "red")) +
      theme_void() +
      theme(legend.position="none")
ggsave("figures/divergence/divergence.png", p, width=8, height=5)





