

library(tidyverse)
library(raster)
library(rgdal)
library(alphahull)
library(rgeos)
library(colormap)
library(grid)
library(gridExtra)


setwd("e:/edges/edgy")

spf <- list.files("F:/little_trees/raw_data", 
                  recursive=T, full.names=T, pattern="shp")


f <- list.files("F:/chelsa/derived", full.names=T)
ext <- extent(-173, -50, 13, 80)
diversity <- raster(f[1]) %>% crop(ext)
diversity[] <- 0

diversify <- function(spp){
      # spp <- "Quercus douglasii"
      message(spp)
      readOGR(spp]) %>%
            rasterize(diversity, 1) %>%
            "+"(diversity)
}

for(x in spf) diversity <- diversify(x)


