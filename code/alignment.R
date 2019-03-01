
# how does spatial climate correlation across a landscape compare
# to the climate change vector in climate space?


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



##### sign agreement between temp and precip deltas #####

futures <- list.files("F:/chelsa/cmip5_bio19", full.names=T)

if(F){
      for(i in seq(1, length(futures), 2)){
            
            message(i)
            
            clim <- stack(c("F:/chelsa/bio19/CHELSA_bio10_1.tif",
                            futures[i],
                            "F:/chelsa/bio19/CHELSA_bio10_12.tif",
                            futures[i+1]))
            names(clim) <- c("temp0", "temp1", "ppt0", "ppt1")
            
            # exclude polar regions, where square landscapes are hardly square
            clim <- clim %>%
                  crop(extent(-180, 180, -66, 66)) %>%
                  reclassify(c(-Inf, -2000, NA))
            
            clim[[3]] <- log10(clim[[3]])
            clim[[4]] <- log10(clim[[4]])
            
            writeRaster(clim, "../climate_change.tif", overwrite=T)
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
            
            writeRaster(deltas, paste0("../climate_deltas/", i, ".tif"), overwrite=T)
      }
}



##### sign of spatial correlation between temp and precip #####

# 1 if positive correlation, 0 if negative
if(F){
      correlation <- function(x, ...){
            m <- matrix(x, ncol=2)
            r <- cor(m[,1], m[,2], use="pairwise.complete.obs")
            as.integer(sign(r) == 1)
      }
      clim <- stack("../climate_change.tif")
      names(clim) <- c("temp0", "temp1", "ppt0", "ppt1")
      space <- aggregate(clim[[c(1, 3)]], fact=c(50, 50, 2), fun=correlation)
      writeRaster(space, "../correlations.tif", overwrite=T)
}

##### compared #####

# proportion of GCMs with spatial-temporal match
m <- list.files("../climate_deltas/", full.names=T) %>%
      c("../correlations.tif") %>%
      stack() %>%
      rasterToPoints() %>%
      as.data.frame() %>%
      gather(gcm, delta, -x, -y, -correlations) %>%
      mutate(alignment = correlations == delta) %>%
      group_by(x, y) %>%
      summarize(alignment = mean(alignment),
                delta = mean(delta),
                spatial = mean(correlations))


p <- ggplot(m, aes(x, y, fill=alignment)) +
      geom_raster() +
      #scale_fill_viridis_c(direction=-1) +
      scale_fill_gradientn(colours=c("darkorchid", "khaki", "forestgreen")) +
      theme_void() +
      theme(legend.position=c(.1, .2))
ggsave("figures/divergence/divergence.png", p, width=8, height=5)

gm <- m %>%
      gather(stat, value, -x, -y) %>%
      mutate(stat=factor(stat, levels=c("spatial", "delta", "alignment")))

p <- ggplot(gm, aes(x, y, fill=value)) +
      geom_raster() +
      scale_fill_gradientn(colours=c("darkorchid", "khaki", "forestgreen")) +
      facet_grid(stat ~ .) +
      theme_void() +
      theme(strip.text.y=element_text(angle=-90),
            legend.position="bottom") +
      labs(fill=NULL)
ggsave("figures/divergence/divergence_components.png", p, width=8, height=12)


stop()


############ mahalanobis #########

# mahalanobis distance to future mean climate
if(F){
      md <- function(x, n=50, ...){
            
            if(length(x) != n^2 * 4) return(NA) # edge of domain
            if(all(is.na(x))) return(NA)
            
            m <- matrix(x, ncol=4) %>% na.omit()
            if(nrow(m) < n^2) return(NA)
            
            # standardize based on historic period
            t <- (m[,1:2] - mean(m[,1])) / sd(m[,1])
            p <- (m[,3:4] - mean(m[,3])) / sd(m[,3])
            
            m <- cbind(t[,1], p[,1])
            
            d <- try(mahalanobis(matrix(c(mean(t[,2]), mean(p[,2])), 1),
                                 colMeans(m),
                                 cov(m)))
            
            if(class(d) == "try-error") return(NA)
            return(d)
      }
      #mahal500 <- aggregate(clim, fact=c(500, 500, 4), fun=function(x, ...) md(x, n=500))
      mahal <- aggregate(clim, fact=c(50, 50, 4), fun=function(x, ...) md(x, n=50))
      writeRaster(space, "../mahalanobis.tif", overwrite=T)
}
md <- stack("../mahalanobis.tif")


dt <- (means[[2]] - means[[1]]) / sd(means[[1]][], na.rm=T)
dp <- (means[[4]] - means[[3]]) / sd(means[[3]][], na.rm=T)

s <- stack("data/collinearity_heterogeneity.tif") %>%
      stack(mahal) %>%
      stack(dt) %>% stack(dp)
names(s) <- c("col", "het", "var", "cvg", "md", "dt", "dp")
s <- s %>%
      rasterToPoints() %>%
      as.data.frame() %>%
      na.omit()


p <- ggplot(s, aes(x, y, fill=md)) +
      geom_raster() +
      scale_fill_viridis_c(trans="log10") +
      theme_void() +
      theme(legend.position=c(.1, .5))
ggsave("figures/divergence/mahalanobis.png", p, width=8, height=5)



#### is high collinearity associated with high MD,
#### controlling for heterogeneity and magnitude of change?

s$delta <- sqrt(s$dt^2 + s$dp^2)

fit <- lm(log10(md) ~ col + log10(het) + delta, data=s)
summary(fit)
aov(fit)

# partial residual plot
fit <- lm(log10(md) ~ log10(het) + delta, data=s)
s$md_resid <- fit$residuals
fit <- lm(col ~ log10(het) + delta, data=s)
s$col_resid <- fit$residuals
ggplot(sample_n(s, 10000), aes(col_resid, md_resid)) +
      geom_point() +
      geom_smooth(method="lm")

# semipartial residual plot
fit <- lm(log10(md) ~ log10(het) + delta, data=s)
s$md_resid <- fit$residuals
ggplot(sample_n(s, 10000), aes(col, md_resid)) +
      geom_point() +
      geom_smooth(method="lm")



# is there more alignment in higly collinear areas?
# this confounder could explain low r2
