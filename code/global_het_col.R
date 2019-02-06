
# global search for landscapes with different combinations 
# of climate heterogeneity and collinearity


library(tidyverse)
library(raster)
library(scales)
library(colormap)
library(grid)
library(gridExtra)

setwd("E:/edges/edgy")



##### data setup #####

# load climate data
clim <- stack(c("F:/chelsa/bio19/CHELSA_bio10_1.tif",
                "F:/chelsa/bio19/CHELSA_bio10_12.tif"))
names(clim) <- c("temp", "ppt")

# exclude polar regions, where square landscapes are hardly square
clim <- crop(clim, extent(-180, 180, -66, 66))

# standardize
clim$ppt <- log10(clim$ppt)
clim$temp <- scale(clim$temp)
clim$ppt <- scale(clim$ppt)



##### spatial stats #####

collinearity <- function(x, ...){
      m <- matrix(x, ncol=2)
      cor(m[,1], m[,2], use="pairwise.complete.obs")^2
}

heterogeneity <- function(x, ...){
      if(length(x) != 5000) return(NA) # edge of domain
      if(all(is.na(x))) return(NA)
      
      # reconstitute input rasters
      a <- raster(matrix(x[1:2500], nrow=50))
      b <- raster(matrix(x[2501:5000], nrow=50))
      
      # mean difference betwean each cell and its neighbors
      a <- focal(a, matrix(1, 3, 3), function(x) mean(abs(x - x[5])))
      b <- focal(b, matrix(1, 3, 3), function(x) mean(abs(x - x[5])))
      
      # bivariate euclidean difference
      a <- sqrt(a^2 + b^2)
      
      # mean across landscape
      mean(values(a), ...)
}

coverage <- function(x, ...){
      length(na.omit(x))/length(x)
}

col <- aggregate(clim, fact=c(50, 50, 2), fun=collinearity)
het <- aggregate(clim, fact=c(50, 50, 2), fun=heterogeneity)
cvg <- aggregate(clim, fact=c(50, 50, 2), fun=coverage)

s <- stack(col, het, cvg)
writeRaster(s, "data/collinearity_heterogeneity.tif")



##### plot results #####

names(s) <- c("col", "het", "cvg")
d <- s %>%
      rasterToPoints() %>%
      as.data.frame() %>%
      na.omit()

d$hex <- colors2d(dplyr::select(d, col, het), 
                  c("green", "blue", "black", "yellow"),
                  xtrans="rank", ytrans="rank")

scatter <- ggplot(d, aes(col, het)) + 
      geom_point(color=d$hex, size=3) +
      theme_minimal()

map <- ggplot(d, aes(x, y)) + 
      geom_raster(fill=d$hex) +
      theme_void() +
      coord_fixed()

p <- arrangeGrob(scatter, map, nrow=1)
png("figures/global_het_col.png", width=1500, height=750)
grid.draw(p)
dev.off()



##### some examplar landscapes #####
