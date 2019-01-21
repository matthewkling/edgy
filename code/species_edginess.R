
library(tidyverse)
library(raster)
library(rgdal)
library(alphahull)
library(rgeos)
library(colormap)
library(grid)
library(gridExtra)


source("E:/edges/edginess/ah2sp.R")



# load species range maps
spp <- c("Quercus douglasii", "Pseudotsuga menziesii", 
         "Pinus albicaulis", "Pinus coulteri")
f <- list.files("F:/little_trees/raw_data", 
                recursive=T, full.names=T, pattern="shp")
f <- f[grepl(paste(spp, collapse="|"), f)]

d <- lapply(f, readOGR)
names(d) <- basename(dirname(f))

ext <- lapply(d, extent) %>%
      Reduce("merge", .)



# load climate data
f <- list.files("F:/chelsa/derived", full.names=T)
wb <- stack(f[c(1,2,4)]) %>% crop(ext)
f <- list.files("F:/chelsa/monthly48", full.names=T)
jja <- stack(f[c(34,35,36)]) %>% crop(ext) %>% mean()
djf <- stack(f[c(37,40,41)]) %>% crop(ext) %>% mean()

s <- stack(wb, jja, djf)
names(s) <- c("AET", "CWD", "PPT", "JJA", "DJF")

v <- values(s) 
pca <- prcomp(na.omit(v), center=T, scale=T)
pc <- s[[1:2]]
pc[!is.na(pc)] <- pca$x[,1:2]
names(pc) <- c("pc1", "pc2")

l <- as.data.frame(pca$rotation)
for(i in 1:length(l)) l[[i]] <- pca$sdev[i] * l[[i]]
l$input <- rownames(l)


# realized niche data
niche <- function(x){
      
      cc <- stack(s, pc) %>% crop(x)
      xr <- rasterize(x, cc[[1]])
      xd <- stack(cc, xr) %>%
            rasterToPoints() %>%
            as.data.frame() %>%
            na.omit() %>%
            rename(id=layer)
      
      # species-specific PCA
      spca <- xd %>% dplyr::select(AET:DJF) %>% scale(center=T, scale=T) %>% prcomp()
      spc <- spca$x
      colnames(spc) <- paste0("s", colnames(spc))
      xd <- cbind(xd, spc)
      
      ld <- as.data.frame(spca$rotation)
      for(i in 1:length(ld)) l[[i]] <- spca$sdev[i] * ld[[i]]
      ld$input <- rownames(ld)
      
      x <- list(x=xd,
                loadings=ld)
      
      return(x)
}



# alpha hulls
#box <- r[[2]]
hulls <- function(box){
      require(FNN)
      
      x <- box$x
      
      dupe <- duplicated(x[,c("pc1", "pc2", "sPC1", "sPC2")])
      xu <- x[!dupe,]
      
      dupe <- duplicated(x[,c("x", "y")])
      xg <- x[!dupe,]
      
      if(nrow(xu)>10000){
            rng1 <- diff(range(xu$sPC1))
            rng2 <- diff(range(xu$sPC2))
            xu <- xu %>%
                  mutate(pc1b=plyr::round_any(sPC1, rng1/200),
                         pc2b=plyr::round_any(sPC2, rng2/200)) %>%
                  group_by(pc1b, pc2b) %>%
                  sample_n(1) %>%
                  ungroup()
      }
      if(nrow(xg)>10000){
            rng1 <- diff(range(xg$x))
            rng2 <- diff(range(xg$y))
            xg <- xg %>%
                  mutate(pc1b=plyr::round_any(ecdf(x)(x), .005),
                         pc2b=plyr::round_any(ecdf(y)(y), .005)) %>%
                  group_by(pc1b, pc2b) %>%
                  sample_n(1) %>%
                  ungroup()
      }
      
      # alpha hull
      ah <- ahull(xu[,c("pc1", "pc2")], alpha=.5)
      ah2 <- ahull(xu[,c("sPC1", "sPC2")], alpha=.5)
      ahg <- ahull(xg[,c("x", "y")], alpha=2)
      sp <- ah2sp(ah)
      sp2 <- ah2sp(ah2)
      spg <- ah2sp(ahg, rnd=3)
      
      # sample-hull distances
      p <- xu
      p2 <- xu
      pg <- xg
      coordinates(p) <- c("pc1", "pc2")
      coordinates(p2) <- c("sPC1", "sPC2")
      coordinates(pg) <- c("x", "y")
      dst <- as.vector(gDistance(p, as(sp, "SpatialLines"), byid=T))
      dst2 <- as.vector(gDistance(p2, as(sp2, "SpatialLines"), byid=T))
      dstg <- as.vector(gDistance(pg, as(spg, "SpatialLines"), byid=T))
      
      # point-sample neighbors
      nn <- get.knnx(xu[,c("pc1", "pc2")], x[,c("pc1", "pc2")], k=1)
      nn2 <- get.knnx(xu[,c("sPC1", "sPC2")], x[,c("sPC1", "sPC2")], k=1)
      nng <- get.knnx(xg[,c("x", "y")], x[,c("x", "y")], k=1)
      x$dst <- dst[nn$nn.index[,1]]
      x$dst2 <- dst2[nn2$nn.index[,1]]
      x$dstg <- dstg[nng$nn.index[,1]]
      
      
      
      box$x <- x
      box$sp <- sp
      box$sp2 <- sp2
      box$spg <- spg
      return(box)
}



plots <- function(box){
      #box <- r[[2]]
      x <- box$x
      
      
      ### global PC space ###
      
      # continuous colors
      fade <- scales::rescale(x$dst, c(1,0)) ^ 3
      fadeg <- scales::rescale(x$dstg, c(1,0)) ^ 3
      
      x$hex <- x %>%
            dplyr::select(pc1, pc2) %>%
            colorwheel2d(colors = c("gray80", "yellow", "green", "cyan", "blue",
                                    "magenta", "red"),
                         kernel=function(x) x^.01) %>%
            col2rgb() %>%
            t() %>% 
            "*"(fade) %>%
            rgb(maxColorValue=255)
      x$hexg <- x %>%
            dplyr::select(x, y) %>%
            colorwheel2d(colors = c("gray80", "yellow", "green", "cyan", "blue",
                                    "magenta", "red"),
                         kernel=function(x) x^.01) %>%
            col2rgb() %>%
            t() %>% 
            "*"(fadeg) %>%
            rgb(maxColorValue=255)
      
      x <- sample_n(x, nrow(x))
      
      sd <- fortify(box$sp)
      sdg <- fortify(box$spg)
      
      lengthen <- 3
      
      p1 <- ggplot(x, aes(pc1, pc2)) + 
            geom_point(color=x$hex) +
            geom_polygon(data=sd, aes(long, lat),
                         fill=NA, color="black") +
            geom_segment(data=l, aes(x=0, y=0, xend=PC1*lengthen, yend=PC2*lengthen),
                         color="red") +
            geom_text(data=l, aes(x=PC1*lengthen+.15, y=PC2*lengthen+.15, label=input), 
                      color="red", hjust=.5, vjust=1) +
            theme_minimal() +
            coord_fixed()
      
      p1g <- ggplot(x, aes(pc1, pc2)) + 
            geom_point(color=x$hexg) +
            geom_segment(data=l, aes(x=0, y=0, xend=PC1*lengthen, yend=PC2*lengthen),
                         color="red") +
            geom_text(data=l, aes(x=PC1*lengthen+.15, y=PC2*lengthen+.15, label=input), 
                      color="red", hjust=.5, vjust=1) +
            theme_minimal() +
            coord_fixed()
      
      p2 <- ggplot(x, aes(x, y)) +
            geom_raster(fill=x$hex) +
            theme_void() +
            coord_fixed(ratio=1.2)
      
      p2g <- ggplot(x, aes(x, y)) +
            geom_raster(fill=x$hexg) +
            geom_polygon(data=sdg, aes(long, lat, group=piece, order=order),
                         fill=NA, color="black") +
            theme_void() +
            coord_fixed(ratio=1.2)
      
      p <- arrangeGrob(p1, p1g, ncol=1)
      p <- arrangeGrob(p2, p, p2g, ncol=3, widths=c(2,1,2))
      p <- arrangeGrob(textGrob(label=x$spp[1]), p, ncol=1, heights=c(1, 10))
      
      png(paste0("E:/edges/edginess/species/", x$spp[1], ".png"), 
          width=2000, height=800)
      grid.draw(p)
      dev.off()
      
      return("boom")
      
      ### specific PC space ###
      
      # continuous colors
      fade <- scales::rescale(x$dst2, c(1,0)) ^ 3
      
      x$hex <- x %>%
            dplyr::select(sPC1, sPC2) %>%
            colorwheel2d(colors = c("gray80", "yellow", "green", "cyan", "blue",
                                    "magenta", "red"),
                         kernel=function(x) x^.01) %>%
            col2rgb() %>%
            t() %>% 
            "*"(fade) %>%
            rgb(maxColorValue=255)
      
      
      
      lengthen <- 3
      p1 <- ggplot(x, 
                   aes(sPC1, sPC2)) + 
            geom_point(color=x$hex) +
            geom_segment(data=box$loadings, aes(x=0, y=0, xend=PC1*lengthen, yend=PC2*lengthen),
                         color="red") +
            geom_text(data=box$loadings, aes(x=PC1*lengthen+.15, y=PC2*lengthen+.15, label=input), 
                      color="red", hjust=.5, vjust=1) +
            theme_minimal() +
            coord_fixed()
      
      p2 <- ggplot(x, 
                   aes(x, y)) +
            geom_raster(fill=x$hex) +
            theme_minimal() +
            coord_fixed(ratio=1.2)
      
      p <- arrangeGrob(textGrob(label=x$spp[1]), 
                       arrangeGrob(p2, p1, nrow=1, widths=c(3, 1)), 
                       ncol=1, heights=c(1, 10))
      png(paste0("E:/edges/edginess/species/", x$spp[1], " spc.png"), 
          width=1000, height=1000)
      grid.draw(p)
      dev.off()
}

r <- lapply(d, niche)
r <- lapply(r, hulls)
for(i in 1:length(r)) r[[i]]$x$spp <- names(r)[i]
lapply(r, plots)




