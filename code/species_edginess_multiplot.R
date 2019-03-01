
# manuscript figure 3


library(tidyverse)
library(raster)
library(rgdal)
library(alphahull)
library(rgeos)
library(colormap)
library(grid)
library(gridExtra)


done <- list.files("f:/edges/species_data/", full.names=T)


if(F){
      # state borders, for plotting
      usa <- getData("GADM", country="USA", level=1) %>%
            #spTransform(crs(s)) %>%
            fortify() %>%
            mutate(group=paste("usa", group))
      can <- getData("GADM", country="CAN", level=1) %>%
            #spTransform(crs(s)) %>%
            fortify() %>%
            mutate(group=paste("can", group))
      mex <- getData("GADM", country="MEX", level=1) %>%
            #spTransform(crs(s)) %>%
            fortify() %>%
            mutate(group=paste("mex", group))
      borders <- rbind(usa, can, mex)
}

world <- map_data("world")

species <- c("Pinus albicaulis", "Juniperus osteosperma", 
             "Quercus lobata", "Acer rubrum")#, "Rhus microphylla")

if(! all(sapply(species, function(spp) any(grepl(spp, done))))) stop("problematic species list")


################## species plots ######################

spp_plots <- function(spp){
      #spp <- species[2]
      box <- readRDS(done[grepl(spp, done)])
      
      x <- box$x
      
      # assign 2d color ramps
      x$g <- 0
      x$b <- scales::rescale(x$dst, 1:0) ^ 4
      x$r <- scales::rescale(x$dstg, 1:0) ^ 4
      x$hex <- rgb(dplyr::select(x, r, g, b), maxColorValue=1)
      
      x <- sample_n(x, nrow(x))
      
      
      sd <- fortify(box$sp)
      sdg <- fortify(box$spg)
      l <- box$loadings
      
      lengthen <- 5
      
      require(shadowtext)
      
      xs <- x %>%
            mutate(dst=plyr::round_any(dst, diff(range(dst))/20),
                   dstg=plyr::round_any(dstg, diff(range(dstg))/20)) %>%
            group_by(dst, dstg) %>%
            summarize(n=n(),
                      hex=hex[1])
      
      scatter <- ggplot(xs, aes(dst, dstg, size=n)) + 
            geom_point(color=xs$hex) +
            geom_vline(xintercept=0, color="#0080ff", size=1) +
            geom_hline(yintercept=0, color="#ff8000", size=1) +
            #annotate(geom="text", x=max(x$dst)*.05, y=max(x$dstg)*.95, 
            #         color="black", size=5, hjust=0,
            #         label=paste0("r=", round(cor(x$dst, x$dstg), 2))) +
            theme_minimal() +
            theme(legend.position="none") +
            theme(axis.title=element_text(size=20, vjust=0),
                  axis.text=element_blank(), axis.ticks=element_blank()) +
            labs(x = "dist. to climate edge (stdev.)",
                 y = "dist. to geographic edge (deg.)")
      
      climate <- ggplot(x, aes(PC1, PC2)) + 
            geom_point(color=x$hex) +
            geom_polygon(data=sd, aes(long, lat),
                         fill=NA, color="#0080ff", size=1) +
            #geom_segment(data=l, aes(x=0, y=0, xend=PC1*lengthen, yend=PC2*lengthen),
            #             color="gray80") +
            #geom_shadowtext(data=l, aes(x=PC1*lengthen+.15, y=PC2*lengthen+.15, label=input), 
            #                color="black", bg.colour="white", hjust=.5, vjust=1) +
            theme_minimal() +
            theme(axis.title=element_text(size=20, vjust=0),
                  axis.text=element_blank(), axis.ticks=element_blank()) +
            coord_fixed()
      
      mag <- max(diff(range(x$x)), diff(range(x$y)))
      ar <- 1
      map <- ggplot() +
            #geom_polygon(data=borders, aes(long, lat, group=group),
            #             fill="gray80", color="white") +
            geom_polygon(data=world, aes(long, lat, group=group),
                         fill="gray80") +
            geom_raster(data=x, aes(x, y), 
                        fill=x$hex) +
            geom_polygon(data=sdg, aes(long, lat, group=piece),
                         fill=NA, color="#ff8000", size=1) +
            annotate(geom="text", size=15,
                     x=mean(range(x$x))-.45*mag, 
                     y=mean(range(x$y))+.45*mag,
                     label=letters[match(spp, species)]) +
            annotate(geom="text", size=8,
                     x=mean(range(x$x))+.45*mag, 
                     y=mean(range(x$y))+.45*mag,
                     label=spp, hjust=1) +
            theme_void() +
            scale_x_continuous(expand=c(0,0)) +
            scale_y_continuous(expand=c(0,0)) +
            coord_fixed(ratio=ar, 
                        xlim=mean(range(x$x))+mag*c(-.5,.5)*ar, 
                        ylim=mean(range(x$y))+mag*c(-.5,.5))
      
      p <- arrangeGrob(climate, scatter, ncol=1)
      p <- arrangeGrob(map, p, ncol=2, widths=c(2,1))
      
      #png(paste0("figures/species_v2/", box$spp, ".png"), 
      #    width=1200, height=800)
      #grid.draw(p)
      #dev.off()
      return(p)
}
s <- lapply(species, spp_plots)


##################### key map ########################

bbox <- function(spp){
      box <- readRDS(done[grepl(spp, done)])
      x <- box$x
      expand.grid(x=range(x$x), y=range(x$y)) %>%
            mutate(spp=spp,
                   order=c(1,2,4,3)) %>%
            arrange(order)
}
k <- species %>%
      lapply(bbox) %>%
      do.call("rbind", .)

lett <- k %>%
      mutate(letter=letters[match(spp, species)]) %>%
      group_by(letter) %>%
      summarize(x=min(x), y=max(y))


key <- ggplot() + 
      geom_polygon(data=world, aes(long, lat, group=group),
                   fill="gray80") +
      geom_polygon(data=k, 
                   aes(x, y, group=spp),
                   fill=NA, color="black", size=.5) +
      geom_text(data=lett, aes(x, y, label=letter),
                nudge_x=2, nudge_y=-1, size=8) +
      coord_map(projection="ortho", 
                orientation=c(20, mean(range(k$x)), 0),
                ylim=c(0, 90)) +
      scale_y_continuous(breaks=seq(-90, 90, 15)) +
      scale_x_continuous(breaks=seq(-180, 180, 15)) +
      theme_minimal() +
      theme(axis.title=element_blank(), axis.text=element_blank(),
            panel.grid=element_line(color="gray85", size=1))



################## correlation histogram ###############

r <- read.csv("e:/edges/edgy/data/stats.txt", stringsAsFactors=F) %>%
      group_by(species) %>%
      slice(1)
rs <- r[r$species %in% species,] %>%
      ungroup() %>%
      mutate(letter=letters[1:nrow(.)])

hst <- ggplot() + 
      geom_vline(xintercept=0, color="gray", size=1) +
      geom_histogram(data=r, aes(cor_pearson), alpha=.3, bins=20) +
      geom_segment(data=rs, aes(x=cor_pearson, xend=cor_pearson,
                                y=0, yend=10)) +
      geom_text(data=rs, aes(x=cor_pearson, y=15, label=letter),
                size=8) +
      theme_minimal() +
      theme(axis.title=element_text(size=20),
            axis.text=element_text(size=20)) +
      labs(x=paste0("correlation between distances to climatic and geographic edges\n(n = ",
                    length(unique(r$species)), " tree species)"),
           y="number of species")


################ manual color legend ################

lpd <- data.frame(x=c(.25, .25, 
                      1, 1, 1, 1), 
                  y=c(4, 3,
                      4, 3, 2, 1), 
                  text=c("geographic edge", "climate edge", 
                         "near geographic edge", "near climate edge", "near both edges", "near neither edge"),
                  color=c("#ff8000", "#0080ff", 
                          "red", "blue", "magenta", "black"),
                  symbol=c("line", "line", 
                           "point", "point", "point", "point"))

legend <- ggplot() +
      geom_point(data=lpd, aes(x, y, shape=symbol, size=symbol), 
                 color=lpd$color) +
      geom_text(data=lpd, aes(x, y, label=paste("  ", text)), 
                color=lpd$color, hjust=0, size=10) +
      scale_shape_manual(values=c(95, 15)) +
      scale_size_manual(values=c(20, 5)) +
      #annotate(geom="point", x=2, y=1, color="white") +
      xlim(0, 2) + ylim(-3, 5) +
      theme_void() +
      theme(legend.position="none")


################ assemble panels ##################

p <- arrangeGrob(arrangeGrob(s[[1]], s[[1]], nrow=1),
                 textGrob(label=""),
                 arrangeGrob(s[[3]], s[[4]], nrow=1),
                 textGrob(label=""),
                 arrangeGrob(arrangeGrob(legend, hst, nrow=2, heights=c(1, 2)),
                             key, nrow=1),
                 ncol=1, heights=c(10, 1, 10, 1, 10))

png("figures/species_edges/multispecies.png", 
    width=2000, height=2000)
grid.draw(p)
dev.off()