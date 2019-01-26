
library(tidyverse)
library(raster)
library(rgdal)
library(alphahull)
library(rgeos)
library(colormap)
library(grid)
library(gridExtra)


setwd("e:/edges/edgy")

source("code/ah2sp.R")



# load species range maps
spf <- list.files("F:/little_trees/raw_data", 
                  recursive=T, full.names=T, pattern="shp")


# load climate data
f <- list.files("F:/chelsa/derived", full.names=T)
ext <- extent(-173, -50, 13, 72)
wb <- stack(f[c(1,2,4)]) %>% crop(ext)
f <- list.files("F:/chelsa/monthly48", full.names=T)
jja <- stack(f[c(34,35,36)]) %>% crop(ext) %>% mean()
djf <- stack(f[c(37,40,41)]) %>% crop(ext) %>% mean()

s <- stack(wb, jja, djf)
names(s) <- c("AET", "CWD", "PPT", "JJA", "DJF")


# state borders, for plotting
usa <- getData("GADM", country="USA", level=1) %>%
      spTransform(crs(s)) %>%
      fortify() %>%
      mutate(group=paste("usa", group))
can <- getData("GADM", country="CAN", level=1) %>%
      spTransform(crs(s)) %>%
      fortify() %>%
      mutate(group=paste("can", group))
mex <- getData("GADM", country="MEX", level=1) %>%
      spTransform(crs(s)) %>%
      fortify() %>%
      mutate(group=paste("mex", group))
borders <- rbind(usa, can, mex)



niche <- function(spp){
      # spp <- "Quercus douglasii"
      
      rng <- readOGR(spf[grepl(spp, spf)])
      xd <- s %>% 
            crop(rng) %>% 
            mask(rng) %>%
            rasterToPoints() %>%
            as.data.frame() %>%
            na.omit()
      
      # species-specific PCA
      spca <- xd %>% dplyr::select(AET:DJF) %>% scale(center=T, scale=T) %>% prcomp()
      spc <- spca$x
      #colnames(spc) <- paste0("s", colnames(spc))
      xd <- cbind(xd, spc)
      
      ld <- as.data.frame(spca$rotation)
      #for(i in 1:length(ld)) l[[i]] <- spca$sdev[i] * ld[[i]]
      ld$input <- rownames(ld)
      
      x <- list(spp=spp,
                x=xd,
                loadings=ld)
      return(x)
}



# alpha hulls and edginess
edges <- function(box){
      #box <- r[[2]]
      
      require(FNN)
      
      x <- box$x
      
      dupe <- duplicated(x[,c("PC1", "PC2")])
      xu <- x[!dupe,]
      
      dupe <- duplicated(x[,c("x", "y")])
      xg <- x[!dupe,]
      
      if(nrow(xu)>10000){
            rng1 <- diff(range(xu$PC1))
            rng2 <- diff(range(xu$PC2))
            xu <- xu %>%
                  mutate(PC1b=plyr::round_any(PC1, rng1/200),
                         PC2b=plyr::round_any(PC2, rng2/200)) %>%
                  group_by(PC1b, PC2b) %>%
                  sample_n(1) %>%
                  ungroup()
      }
      if(nrow(xg)>10000){
            rng1 <- diff(range(xg$x))
            rng2 <- diff(range(xg$y))
            xg <- xg %>%
                  mutate(PC1b=plyr::round_any(ecdf(x)(x), .005),
                         PC2b=plyr::round_any(ecdf(y)(y), .005)) %>%
                  group_by(PC1b, PC2b) %>%
                  sample_n(1) %>%
                  ungroup()
      }
      
      # alpha hull
      ah <- ahull(xu[,c("PC1", "PC2")], alpha=.5)
      ahg <- ahull(xg[,c("x", "y")], alpha=2)
      sp <- ah2sp(ah)
      spg <- ah2sp(ahg, rnd=3)
      
      # sample-hull distances
      p <- xu
      pg <- xg
      coordinates(p) <- c("PC1", "PC2")
      coordinates(pg) <- c("x", "y")
      dst <- as.vector(gDistance(p, as(sp, "SpatialLines"), byid=T))
      dstg <- as.vector(gDistance(pg, as(spg, "SpatialLines"), byid=T))
      
      # point-sample neighbors
      nn <- get.knnx(xu[,c("PC1", "PC2")], x[,c("PC1", "PC2")], k=1)
      nng <- get.knnx(xg[,c("x", "y")], x[,c("x", "y")], k=1)
      x$dst <- dst[nn$nn.index[,1]]
      x$dstg <- dstg[nng$nn.index[,1]]
      
      box$x <- x
      box$sp <- sp
      box$spg <- spg
      box$cor_pearson <- cor(x$dst, x$dstg, use="pairwise.complete.obs", method="pearson")
      box$cor_spearman <- cor(x$dst, x$dstg, use="pairwise.complete.obs", method="spearman")
      
      write(paste0(box$spp, ", ", box$cor_pearson, ", ", box$cor_spearman), 
            file="e:/edges/edgy/data/stats.txt", append=T)
      saveRDS(box, paste0("f:/edges/species_data/", box$spp, ".rds"))
      
      return(box)
}



plots <- function(box){
      #box <- r[[4]]
      x <- box$x
      
      
      ### global PC space ###
      
      # continuous colors
      fade <- scales::rescale(x$dst, c(1,0)) ^ 3
      fadeg <- scales::rescale(x$dstg, c(1,0)) ^ 3
      
      x$hex <- x %>%
            dplyr::select(PC1, PC2) %>%
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
      
      l <- box$loadings
      
      lengthen <- 5
      
      require(shadowtext)
      
      p1 <- ggplot(x, aes(PC1, PC2)) + 
            geom_point(color=x$hex) +
            geom_polygon(data=sd, aes(long, lat),
                         fill=NA, color="black") +
            geom_segment(data=l, aes(x=0, y=0, xend=PC1*lengthen, yend=PC2*lengthen),
                         color="gray80") +
            geom_shadowtext(data=l, aes(x=PC1*lengthen+.15, y=PC2*lengthen+.15, label=input), 
                            color="black", bg.colour="white", hjust=.5, vjust=1) +
            theme_minimal() +
            coord_fixed()
      
      p1g <- ggplot(x, aes(PC1, PC2)) + 
            geom_point(color=x$hexg) +
            geom_segment(data=l, aes(x=0, y=0, xend=PC1*lengthen, yend=PC2*lengthen),
                         color="gray80") +
            geom_shadowtext(data=l, aes(x=PC1*lengthen+.15, y=PC2*lengthen+.15, label=input), 
                            color="black", bg.colour="white", hjust=.5, vjust=1) +
            theme_minimal() +
            coord_fixed()
      
      p2 <- ggplot(x, aes(x, y)) +
            geom_polygon(data=borders, aes(long, lat, group=group),
                         fill="gray80", color="white") +
            geom_raster(fill=x$hex) +
            theme_void() +
            coord_fixed(ratio=1.2, xlim=range(x$x), ylim=range(x$y))
      
      p2g <- ggplot(x, aes(x, y)) +
            geom_polygon(data=borders, aes(long, lat, group=group),
                         fill="gray80", color="white") +
            geom_raster(fill=x$hexg) +
            geom_polygon(data=sdg, aes(long, lat, group=piece, order=order),
                         fill=NA, color="black") +
            theme_void() +
            coord_fixed(ratio=1.2, xlim=range(x$x), ylim=range(x$y))
      
      p <- arrangeGrob(p1, p1g, ncol=1)
      p <- arrangeGrob(p2, p, p2g, ncol=3, widths=c(2,1,2))
      p <- arrangeGrob(textGrob(label=box$spp, gp=gpar(fontsize=40)), 
                       p, ncol=1, heights=c(1, 10))
      
      png(paste0("figures/species/", box$spp, ".png"), 
          width=2000, height=800)
      grid.draw(p)
      dev.off()
}


#sp <- c("Quercus douglasii", "Pseudotsuga menziesii", "Pinus albicaulis", "Pinus coulteri")
#r <- sp[4] %>% niche() %>% hulls() %>% plots()

sp <- basename(dirname(spf))
#done <- list.files("figures/species") %>% sub("\\.png", "", .)
done <- read.csv("e:/edges/edgy/data/stats.txt", stringsAsFactors=F)$species
sp <- sp[! sp %in% done]

library(doParallel)
cl <- makeCluster(4)
registerDoParallel(cl)

r <- foreach(x = sp) %dopar% {
      
      library(tidyverse)
      library(raster)
      library(rgdal)
      library(alphahull)
      library(rgeos)
      library(colormap)
      library(grid)
      library(gridExtra)
      
      try(x %>% niche() %>% edges())# %>% plots())
}

stopCluster(cl)





## histogram of correlation coefficients ##

r <- read.csv("e:/edges/edgy/data/stats.txt", stringsAsFactors=F) %>%
      gather(stat, value, -species) %>%
      mutate(stat=sub("cor_", "", stat))
      
p <- ggplot(r, aes(value, color=stat, fill=stat)) + 
      geom_vline(xintercept=0) +
      geom_density(alpha=.3) +
      theme_minimal()+
      theme(legend.position="bottom") +
      labs(x=paste0("correlation between climatic and geographic edginess\n(n =",
                    length(unique(r$species)), " tree species)"))

ggsave("figures/edginess_cor_hist.png", p, width=6, height=4, units="in")
