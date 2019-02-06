
# geoclimate manifolds


library(dplyr)
library(tidyr)
library(raster)
library(sampSurf)
library(geosphere)
library(ggplot2)
library(grid)
library(gridExtra)
library(rgeos)


setwd("E:/edges/manifold")

# load PRISM rasters
f <- list.files("data", pattern="annual_asc.asc", recursive=T, full.names=T)
s <- stack(f[!grepl("aux", f)])

# log-transform precip
s[[1]] <- log(s[[1]])

# crop to focal areas
circle <- function(x,y,r){
      c <- data.frame(x=x, y=y)
      coordinates(c) <- c("x", "y")
      crs(c) <- crs(s)
      gBuffer(c, width=r, quadsegs=20)
}
crop_circle <- function(p, r) mask(crop(r, p), p)
regions <- list(
      thchp=circle(-118.5, 36, 2), # tehachapi
      grbsn=circle(-115.3, 39.4, 2), # great basin
      grpla=circle(-98, 41, 2), # great plains
      #sfthl=circle(-120.4, 38.5, .5)) # sierra foothills
      sfthl=circle(-120.8, 39.0, .6)) # sierra foothills
r <- lapply(regions, crop_circle, r=s)

# generate species
makespp <- function(s, occ){
      v <- scale(values(s))
      d <- apply(v, 1, function(x)sqrt(sum(x^2)))
      d <- as.integer(ecdf(d)(d)<occ)
      r <- s[[1]]
      r[] <- d
      s <- stack(s, r)
      names(s) <- c("temp", "precip", "occ")
      s
}
r <- lapply(r, makespp, occ=.1)


# reduce size by sampling points
p <- lapply(r, sampleRandom, size=1000, sp=T, xy=T)

# pairwise distances
d <- list(geog=lapply(p, distm),
          clim=lapply(p, function(x)as.matrix(dist(scale(x@data[,3:4])))))

# merge and restructure distance data
dm2df <- function(x){
      x <- as.data.frame(x)
      px <- paste0("p", 1:ncol(x))
      colnames(x) <- px
      x$px1 <- px
      x <- gather(x, px2, dist, -px1)
      return(x)
}
f <- lapply(d, function(x)lapply(x, dm2df))
for(space in names(f)){
      for(region in names(f[[1]])){
            f[[space]][[region]]$space <- space
            f[[space]][[region]]$region <- region
      }
}
f <- lapply(f, function(x)Reduce("rbind",x))
f <- Reduce("rbind", f)

# join to raw data
pd <- lapply(p, function(x){
      x <- x@data
      x$px <- paste0("p", 1:nrow(x))
      return(x)
})
for(region in names(pd)) pd[[region]]$region <- region
pd <- Reduce("rbind", pd)
pd1 <- pd
pd2 <- pd
names(pd1) <- paste0(names(pd1), 1)
names(pd2) <- paste0(names(pd2), 2)
names(pd1) <- sub("region1", "region", names(pd1))
names(pd2) <- sub("region2", "region", names(pd2))

f <- left_join(f, pd1)
f <- left_join(f, pd2)




stop("wootwoot")


eb <- element_blank()
style <- theme(panel.background=eb, panel.grid=eb, 
               axis.text=eb, axis.ticks=eb, axis.title=eb,
               title=element_text(size=30),
               legend.position="none")


# semivariogram data
pd <- spread(f, space, dist)
pd$occ12 <- paste0(pd$occ1, pd$occ2)

# semivariogram plot
semivar <- function(location){
      pdv <- filter(pd, region==location)
      return(ggplot() +
                   geom_point(data=sample_n(pdv, 10000), 
                              aes(geog, clim, color=occ12),
                              size=2, alpha=.3, shape=16) +
                   geom_point(data=filter(pdv, occ12=="11") %>% sample_n(1000), 
                              aes(geog, clim), 
                              size=2, alpha=.3, shape=16) +
                   geom_smooth(data=sample_n(pdv, 10000) %>% as.data.frame(), 
                               aes(geog, clim), method="loess",
                               color="goldenrod", size=1.5, se=F) +
                   geom_smooth(data=filter(pdv, occ12=="11") %>% sample_n(1000) %>% as.data.frame(), 
                               aes(geog, clim), method="loess",
                               color="red3", size=1.5, se=F) +
                   scale_color_manual(values=c("cadetblue", "blue", "blue", "black")) +
                   style)
}


# raster data
rs2df <- function(x, ...){
      require(colormap)
      x <- as.data.frame(rasterToPoints(x))
      x$col <- colorwheel2d(scale(as.matrix(x[,3:4])), ...)
      x
}
rdf <- lapply(r, rs2df, colors = c("black", "yellow", "green", "cyan", "blue", "magenta", "red"))
for(region in names(rdf)) rdf[[region]]$region <- region
rdf <- do.call("rbind", rdf)


# geoclim plots
makeplot <- function(location, space){
      rdfv <- filter(rdf, region==location)
      if(space=="clim"){
            return(ggplot() + 
                         geom_point(data=rdfv, aes(temp, precip), 
                                    color=rdfv$col, shape=16) + 
                         geom_point(data=filter(rdfv, occ==1), aes(temp, precip), 
                                    alpha=.1, color="white", shape=16) + 
                         style +
                         labs(title=location))
      }
      if(space=="geog"){
            return(ggplot() + 
                         geom_raster(data=rdfv, aes(x,y), fill=rdfv$col) + 
                         geom_raster(data=filter(rdfv, occ==1), aes(x,y), alpha=.5, fill="white") +
                         style)
      }
}


p <- arrangeGrob(makeplot("grpla", "clim"), makeplot("thchp", "clim"), makeplot("sfthl", "clim"), makeplot("grbsn", "clim"), 
                 makeplot("grpla", "geog"), makeplot("thchp", "geog"), makeplot("sfthl", "geog"), makeplot("grbsn", "geog"),
                 semivar("grpla"), semivar("thchp"), semivar("sfthl"), semivar("grbsn"),
                 nrow=3)
png("charts/geoclim_spp.png", width=2000, height=1500)
grid.draw(p)
dev.off()


###### same chart as above but without species niche

# semivariogram plot
semivar <- function(location){
      pdv <- filter(pd, region==location)
      return(ggplot() +
                   geom_point(data=sample_n(pdv, 10000), 
                              aes(geog, clim),
                              size=2, alpha=.3, shape=16, color="black") +
                   geom_smooth(data=sample_n(pdv, 10000) %>% as.data.frame(), 
                               aes(geog, clim), method="loess",
                               color="red3", size=1.5, se=F) +
                   style)
}


# raster data
rs2df <- function(x, ...){
      require(colormap)
      x <- as.data.frame(rasterToPoints(x))
      x$col <- colorwheel2d(scale(as.matrix(x[,3:4])), ...)
      x
}
rdf <- lapply(r, rs2df, colors = c("black", "yellow", "green", "cyan", "blue", "magenta", "red"))
for(region in names(rdf)) rdf[[region]]$region <- region
rdf <- do.call("rbind", rdf)


# geoclim plots
makeplot <- function(location, space){
      rdfv <- filter(rdf, region==location)
      if(space=="clim"){
            return(ggplot() + 
                         geom_point(data=rdfv, aes(temp, precip), 
                                    color=rdfv$col, shape=16) + 
                         style +
                         labs(title=location))
      }
      if(space=="geog"){
            return(ggplot() + 
                         geom_raster(data=rdfv, aes(x,y), fill=rdfv$col) + 
                         style)
      }
}


p <- arrangeGrob(makeplot("grpla", "clim"), makeplot("sfthl", "clim"), makeplot("grbsn", "clim"), 
                 makeplot("grpla", "geog"), makeplot("sfthl", "geog"), makeplot("grbsn", "geog"),
                 semivar("grpla"), semivar("sfthl"), semivar("grbsn"),
                 nrow=3)
png("charts/geoclim.png", width=1500, height=1500)
grid.draw(p)
dev.off()







# TODO:
#     version with directional dists on univariate climate (requires reference pixel)

