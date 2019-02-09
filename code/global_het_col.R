
# global search for landscapes with different combinations 
# of climate heterogeneity and collinearity


library(tidyverse)
library(raster)
library(scales)
library(colormap)
library(grid)
library(gridExtra)
library(ecoclim)

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

# stash
writeRaster(clim, "../climate_het_col.tif")
clim <- stack("../climate_het_col.tif")
names(clim) <- c("temp", "ppt")


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

variation <- function(x, ...){
      m <- matrix(x, ncol=2)
      mean(apply(m, 2, sd, ...))
}

coverage <- function(x, ...){
      length(na.omit(x))/length(x)
}

col <- aggregate(clim, fact=c(50, 50, 2), fun=collinearity)
het <- aggregate(clim, fact=c(50, 50, 2), fun=heterogeneity)
var <- aggregate(clim, fact=c(50, 50, 2), fun=variation)
cvg <- aggregate(clim, fact=c(50, 50, 2), fun=coverage)

s <- stack(col, het, var, cvg)
writeRaster(s, "data/collinearity_heterogeneity.tif", overwrite=T)
s <- stack("data/collinearity_heterogeneity.tif")



##### plot results #####

id <- s[[1]]
id[] <- 1:ncell(id)
s <- stack(s, id)

names(s) <- c("col", "het", "var", "cvg", "id")
d <- s %>%
      rasterToPoints() %>%
      as.data.frame() %>%
      na.omit()

 
for(vars in list(c("col", "var"), c("col", "het"), c("het", "var"))){
      
      v1 <- vars[1]
      v2 <- vars[2]
      d$v1 <- d[,v1]
      d$v2 <- d[,v2]
      
      d$hex <- colors2d(dplyr::select(d, v1, v2), 
                        c("yellow", "red", "black", "green"))
      
      if(all.equal(vars, c("het", "var")) == T){
            d$hex <- colorwheel2d(dplyr::select(d, v1, v2) %>% mutate_all(rescale),
                                  colors=c("gray20",
                                           "orangered", "yellow", "#8cff00", "green", "#063800",
                                           "black", "darkred", "red"),
                                  origin=c(.25, .2), kernel=sqrt)
     }
     
      scatter <- ggplot(d, aes(v1, v2)) + 
            geom_point(color=d$hex, size=1) +
            theme_minimal() +
            theme(text=element_text(size=25)) +
            labs(x=v1,
                 y=v2)
      
      map <- ggplot(d, aes(x, y)) + 
            geom_raster(fill=d$hex) +
            theme_void()
      
      png(paste0("figures/het_col/global_", v1, "_", v2, ".png"), width=2000, height=1000)
      plot(map)
      print(scatter, 
            vp=viewport(x = 0, y = 0, 
                        width = unit(0.3, "npc"), height = unit(0.6, "npc"),
                        just = c("left", "bottom")))
      dev.off()
      
      
      
      
      ##### some examplar landscapes #####
      
      # function extracts high-res cliamte data for a given landscape
      ls_data <- function(data){
            ext <- extent(data$x-1, data$x+1, data$y-1, data$y+1)
            pts <- crop(clim, ext) %>%
                  rasterToPoints() %>% 
                  as.data.frame() %>% 
                  na.omit()
            coordinates(pts) <- c("x", "y")
            crs(pts) <- crs(clim)
            pts$id <- extract(id, pts)
            pts <- as.data.frame(pts) %>%
                  filter(id==data$id) %>%
                  rename(lsx=x, lsy=y) %>%
                  left_join(data, .)
            return(pts)
      }
      
      # select some landscapes and grab their high-res cliamtes
      e <- d %>%
            filter(cvg > .75) %>%
            mutate(col=rescale(v1),
                   het=rescale(v2)) %>%
            mutate(pos = rescale(col+het),
                   neg = rescale(col-het)) %>%
            filter(abs(pos-.5) > .4 | abs(neg-.5) > .4) %>%
            mutate(heterogeneity = ifelse(het > .5, "high", "low"),
                   collinearity = ifelse(col > .5, "high", "low")) %>%
            group_by(heterogeneity, collinearity) %>%
            sample_n(1) %>%
            ungroup() %>%
            split(1:nrow(.)) %>%
            map(ls_data) %>%
            do.call("rbind", .) %>%
            mutate(heterogeneity=factor(heterogeneity, levels=c("low", "high")),
                   collinearity=factor(collinearity, levels=c("low", "high"))) %>%
            group_by(heterogeneity, collinearity) %>%
            mutate(lsx=rescale(lsx),
                   lsy=rescale(lsy),
                   t=rescale(temp),
                   p=rescale(ppt),
                   group=paste(collinearity, heterogeneity)) %>%
            ungroup() %>%
            arrange(heterogeneity, collinearity)
      
      # bivariate climate color ramps
      # weird ordering needed so manual colors plot on correct facets 
      e <- split(e, e$group)[c(4,2,3,1)] %>% 
            lapply(function(x){
                  xd <- dplyr::select(x, t, p) %>% as.matrix()
                  x$hex <- colorwheel2d(xd, kernel=function(x) x^2)
                  return(x)
            }) %>%
            do.call("rbind", .)
      
      
      lblr1 <- as_labeller(function(x) paste(v1, x, sep=": "))
      lblr2 <- as_labeller(function(x) paste(v2, x, sep=": "))
      
      map <- ggplot(e, aes(lsx, lsy)) +
            geom_raster(fill=e$hex) + 
            facet_wrap(vars(heterogeneity, collinearity), 
                       scales="free", as.table=F, nrow=1,
                       labeller=labeller(heterogeneity=lblr1, 
                                         collinearity=lblr2)) +
            theme_void() +
            theme(legend.position="none",
                  strip.text.y=element_text(angle=-90))
      
      scatter <- ggplot(e, aes(temp, ppt)) +
            geom_point(color=e$hex, size=3) +
            theme_minimal() +
            facet_wrap(vars(heterogeneity, collinearity), 
                       scales="free", as.table=F, nrow=1) +
            theme(legend.position="none",
                  strip.text=element_blank())
      
      p <- arrangeGrob(map, scatter, ncol=1)
      ggsave(paste0("figures/het_col/landscape_", v1, "_", v2, ".png"), p, width=12, height=6)
}


##### 3d map #####

d$hex <- dplyr::select(d, col, het, var) %>%
      colors3d(trans="ecdf")

map <- ggplot(d, aes(x, y)) + 
      geom_raster(fill=d$hex) +
      theme_void()

pd <- pairsData(d, c("col", "het", "var"), "hex")

scatters <- ggplot() +
      geom_point(data=pd, aes(x_value, y_value), color=pd$hex) +
      facet_wrap(vars(x_var, y_var), scales="free", nrow=1,
                 labeller=label_both)

p <- arrangeGrob(map, scatters, ncol=1, heights=c(2,1.3))
ggsave("figures/het_col/global_3d.png", p, width=12, height=10)
