
# Global patterns of landscape climate heterogeneity and collinearity


# to do:
# use geosphere::distm instead of dist



library(tidyverse)
library(raster)
library(scales)
library(colormap)
library(grid)
library(gridExtra)
library(ecoclim)
library(mgcv)

source("code/agg.r") # modified version of raster::aggregate

select <- dplyr::select



##### data setup #####

# load climate data
if(F){
      clim <- stack(c("F:/chelsa/bio19/CHELSA_bio10_1.tif",
                      "F:/chelsa/bio19/CHELSA_bio10_12.tif"))
      names(clim) <- c("temp", "ppt")
      
      # exclude polar regions, where square landscapes are hardly square
      clim <- crop(clim, extent(-180, 180, -66, 66))
      
      # standardize
      clim$ppt <- log10(clim$ppt)
      clim$temp <- scale(clim$temp)
      clim$ppt <- scale(clim$ppt)
      
      writeRaster(clim, "../climate_het_col.tif")
      clim <- stack("../climate_het_col.tif")
      names(clim) <- c("temp", "ppt")
      
      # version with latitude and longitude as layers
      lat <- lon <- clim[[1]]
      ll <- coordinates(lat)
      lon[] <- ll[,1]
      lat[] <- ll[,2]
      rm(ll); gc()
      clim_ll <- stack(clim, lon, lat) %>%
            writeRaster("../climate_ll.tif")
      
      
      
      ##### spatial stats #####
      
      # r-squared of climate variables within landscape
      collinearity <- function(x, ...){
            m <- matrix(x, ncol=2)
            cor(m[,1], m[,2], use="pairwise.complete.obs")^2
      }
      
      # mean climate devaitions of smaller neighborhoods within landscape
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
      
      # mean of univariate standard deviations within landscape
      variation <- function(x, ...){
            m <- matrix(x, ncol=2)
            mean(apply(m, 2, sd, ...))
      }
      
      # proportion of landscape with climate data
      coverage <- function(x, ...){
            length(na.omit(x))/length(x)
      }
      
      # spline of pairwise climate difference ~ geographic distance within landscape
      distances <- function(x, ...){
            n <- 21 # number of points at which to predict spline
            mx <- .41 # max pairwise distance for this neighborhood size, within circle
            if(length(x) != 10000) return(rep(NA, n)) # edge of domain
            
            m <- matrix(x, ncol=4) %>% na.omit()
            if(nrow(m) < 2500) return(rep(NA, n))
            
            gd <- dist(m[,3:4])
            cd <- dist(m[,1:2])
            fit <- data.frame(geo_dist=gd[upper.tri(gd)],
                              clim_dist=cd[upper.tri(cd)]) %>%
                  na.omit() %>%
                  filter(geo_dist <= mx) %>%
                  mutate(gd=ntile(geo_dist, 100)) %>%
                  group_by(gd) %>%
                  sample_n(100) %>%
                  gam(clim_dist ~ s(geo_dist, bs="ps"), data=.)
            
            pred <- data.frame(geo_dist=seq(0, mx, length.out=21))
            y <- as.vector(predict(fit, pred))
            if(length(y) != n) return(rep(NA, n))
            return(y)
      }
      
      
      col <- aggregate(clim, fact=c(50, 50, 2), fun=collinearity)
      het <- aggregate(clim, fact=c(50, 50, 2), fun=heterogeneity)
      var <- aggregate(clim, fact=c(50, 50, 2), fun=variation)
      cvg <- aggregate(clim, fact=c(50, 50, 2), fun=coverage)
      dst <- agg(clim_ll, fact=c(50, 50, 4), fun=distances) # slow -- 24 hrs
      
      s <- stack(col, het, var, cvg)
      writeRaster(s, "data/collinearity_heterogeneity.tif", overwrite=T)
      writeRaster(dst, "data/distance_splines.tif", overwrite=T)
}





################## plot results ####################

s <- stack("data/collinearity_heterogeneity.tif")
splines <- stack("data/distance_splines.tif")
clim <- stack("../climate_het_col.tif")
names(clim) <- c("temp", "ppt")

id <- s[[1]]
id[] <- 1:ncell(id)
s <- stack(s, id)

names(s) <- c("col", "het", "var", "cvg", "id")
d <- s %>%
      rasterToPoints() %>%
      as.data.frame() %>%
      na.omit() %>%
      mutate(het=var, 
             quadrant=paste(col>mean(col), het>mean(het)))

pal <- c("gray", "#00d4ff", "#ff3700", "darkorchid")



### maps and scatterplots for exemplar landscapes ###

# function extracts high-res cliamte data for a given landscape
ls_data <- function(data, climate){
      ext <- extent(data$x-1, data$x+1, data$y-1, data$y+1)
      pts <- crop(climate, ext) %>%
            rasterToPoints() %>%
            as.data.frame() %>%
            na.omit()
      coordinates(pts) <- c("x", "y")
      crs(pts) <- crs(climate)
      pts$id <- extract(id, pts)
      pts <- as.data.frame(pts) %>%
            filter(id==data$id) %>%
            dplyr::rename(lsx=x, lsy=y) %>%
            left_join(data, .)
      return(pts)
}

# categorize landscapes by binary het-col combo
ls <- d %>%
      filter(cvg > .5) %>%
      mutate(heterogeneity = ifelse(het > mean(het), "high", "low"),
             collinearity = ifelse(col > mean(col), "high", "low"),
             # flag corner cases as canditates for examples:
             cl=rescale(col),
             ht=rescale(het),
             pos = rescale(cl+ht),
             neg = rescale(cl-ht),
             #corner = abs(pos-.5) > .45 | abs(neg-.5) > .45) %>%
             corner = (col<.1 | col>.8) & (het<.03 | het>.3)) %>%
      mutate(heterogeneity=factor(heterogeneity, levels=c("low", "high")),
             collinearity=factor(collinearity, levels=c("low", "high")))
#ggplot(ls, aes(col, het, color=corner)) + geom_point()

# select some example landscapes and grab their high-res climates
exemplar <- c(a=40289, b=81455, c=192946, d=47305) 
e <- ls %>%
      #filter(corner) %>%
      #group_by(heterogeneity, collinearity) %>%
      #sample_n(1) %>%
      filter(id %in% exemplar) %>%
      write_csv("figures/het_col/ls_metadata.csv") %>%
      group_by(heterogeneity, collinearity) %>%
      group_map(~ ls_data(.x, clim)) %>%
      group_by(heterogeneity, collinearity) %>%
      mutate(t=rescale(temp),
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
      do.call("rbind", .)# %>%
#group_by(heterogeneity, collinearity) %>%
#mutate(x=rescale(x), y=rescale(y)) %>%
#ungroup()

roman <- c("i", "ii", "iii", "iv")
ROMAN <- to_upper_ascii(roman)

label_letters <- function(x){
      y <- label_both(x, multi_line=F)
      y[[1]] <- roman #letters[1:4]
      return(y)
}

maps <- ggplot(e, aes(lsx, lsy)) +
      geom_raster(fill=e$hex) +
      facet_wrap(heterogeneity~collinearity,
                 scales="free", as.table=F, nrow=1,
                 labeller=label_letters) +
      theme_minimal() +
      theme(legend.position="none",
            legend.text=element_text(color="white"),
            strip.text.x=element_text(face="bold")) +
      labs(x="longitude", y="latitude")

climates <- ggplot(e, aes(temp, ppt)) +
      geom_point(color=e$hex) +
      theme_minimal() +
      facet_wrap(heterogeneity~collinearity,
                 scales="free", as.table=F, nrow=1,
                 labeller=label_letters) +
      theme(legend.position="none",
            strip.text.x=element_text(face="bold")) +
      labs(x="temperature", y="precipitation")



### distance splines ###

add_names <- function(x, nm){names(x) <- nm; return(x)}

ls_cat <- ls %>% select(id, heterogeneity, collinearity)

dd <- splines %>%
      stack(id) %>%
      add_names(c(paste0("v", 1:21), "id")) %>%
      rasterToPoints() %>%
      as.data.frame() %>%
      left_join(ls_cat) %>%
      na.omit() %>%
      gather(geo, clim, v1:v21) %>%
      mutate(geo=sub("v", "", geo) %>% as.integer(),
             geo=rescale(geo, c(0,.41)))

dde <- filter(dd, id %in% e$id)

dd <- dd %>%
      group_by(heterogeneity, collinearity, geo) %>%
      summarize(p01=quantile(clim, .005),
                p05=quantile(clim, .05),
                p25=quantile(clim, .25),
                p50=quantile(clim, .50),
                p75=quantile(clim, .75),
                p95=quantile(clim, .95),
                p99=quantile(clim, .995)) %>%
      mutate(group=paste(collinearity!="low", heterogeneity!="low"))

dde_txt <- dde %>%
      group_by(heterogeneity, collinearity) %>%
      arrange(geo) %>%
      slice(16) %>%
      ungroup() %>%
      mutate(txt=paste0(roman, "\n"))

freqs <- d %>% 
      count(quadrant) %>% 
      mutate(percent = round(n/sum(n)*100),
             quadrant = ROMAN[c(1,3,2,4)]) %>%
      arrange(quadrant)

lb <- function(x){
      y <- label_both(x)
      y[[1]] <- paste0(freqs$quadrant, "\n", y[[1]], "\n", y[[2]], "\n(", freqs$percent, "% of globe)")
      return(y[1])
}

distances <- ggplot() +
      geom_ribbon(data=dd, aes(x=geo, ymin=p01, ymax=p99, fill=group), alpha=.25) +
      geom_ribbon(data=dd, aes(x=geo, ymin=p05, ymax=p95, fill=group), alpha=.25) +
      geom_ribbon(data=dd, aes(x=geo, ymin=p25, ymax=p75, fill=group), alpha=.25) +
      geom_line(data=dd, aes(x=geo, y=p50, color=group), size=1) +
      geom_line(data=dde, aes(x=geo, y=clim), size=1, color="black", linetype=3) +
      geom_text(data=dde_txt, aes(x=geo, y=clim, label=txt)) +
      scale_fill_manual(values=pal) +
      scale_color_manual(values=pal) +
      theme_minimal() +
      facet_wrap(heterogeneity~collinearity, labeller=lb, #scales="free",
                 as.table=F, nrow=1) +
      theme(legend.position="none",
            strip.text.x=element_text(face="bold")) +
      labs(x="pairwise geographic distance (degrees)",
           y="climatic difference")



### global map and scatterplot ###

#d$hex <- colors2d(dplyr::select(d, col, het),
#                  c("yellow", "red", "black", "green"))

ltrs <- expand.grid(sx=mean(d$col) + c(-5, 5) * diff(range(d$col))/15,
                    sy=mean(d$het) + c(-1, 1) * diff(range(d$het))/15) %>%
      mutate(lab = ROMAN,
             heterogeneity=c("low", "low", "high", "high"),
             collinearity=c("low", "high", "low", "high")) %>%
      left_join(dde %>% select(-geo, -clim) %>% distinct()) %>%
      left_join(d)

global_scatter <- ggplot(d, aes(col, het,
                                color=quadrant)) +
      geom_point(size=1) +
      #geom_vline(xintercept=mean(d$col), color="white") +
      #geom_hline(yintercept=mean(d$het), color="white") +
      geom_text(data=ltrs, aes(sx, sy, label=lab), size=4, color="black") +
      geom_text(data=ltrs, aes(col, het, label=to_lower_ascii(lab)), size=5, color="black") +
      theme_minimal() +
      scale_color_manual(values=pal) +
      #theme(text=element_text(size=25)) +
      theme(legend.position="none") +
      scale_x_continuous(breaks=0:1) +
      scale_y_continuous(breaks=0:1) +
      coord_fixed(ratio=diff(range(d$col))/diff(range(d$het))) +
      labs(x="collinearity",
           y="heterogeneity")


global_map <- ggplot() +
      #geom_raster(fill=d$hex) +
      geom_raster(data=d, aes(x, y, fill=paste(col>mean(col), het>mean(het)))) +
      geom_point(data=ltrs, aes(x, y)) +
      geom_text(data=ltrs, aes(x, y, label=to_lower_ascii(lab)), size=5, nudge_x=8) +
      scale_fill_manual(values=pal) +
      theme_void() +
      theme(legend.position="none")



### assemble multipanel figure ###

hts <- c(1.5, 1.2, 1.1, 1.1)
p <- arrangeGrob(global_map, global_scatter, nrow=1, widths=c(3, 1))
p <- arrangeGrob(p, distances, climates, maps, ncol=1, heights=hts)
#ggsave("figures/het_col/het_col.png", p, width=8, height=10)

# modified version of ggsave
ggs <- function (filename, plot = last_plot(), device = NULL, path = NULL, 
          scale = 1, width = NA, height = NA, units = c("in", "cm", 
                                                        "mm"), dpi = 300, limitsize = TRUE, ...) 
{
      source("https://raw.githubusercontent.com/tidyverse/ggplot2/master/R/save.r")
      source("https://raw.githubusercontent.com/tidyverse/ggplot2/master/R/utilities.r")
      dpi <- parse_dpi(dpi)
      dev <- plot_dev(device, filename, dpi = dpi)
      dim <- plot_dim(c(width, height), scale = scale, units = units, 
                      limitsize = limitsize)
      if (!is.null(path)) {
            filename <- file.path(path, filename)
      }
      old_dev <- grDevices::dev.cur()
      dev(filename = filename, width = dim[1], height = dim[2], 
          ...)
      on.exit(utils::capture.output({
            grDevices::dev.off()
            if (old_dev > 1) grDevices::dev.set(old_dev)
      }))
      grid.draw(plot)
      grid.text(LETTERS[1:4], 
                x=.02, 
                y=c(.98, .66, .42, .20),
                gp=gpar(fontsize=20, fontface="bold", col="black"))
      invisible()
}
ggs("figures/het_col/het_col.png", p, width=8, height=10)








