
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

futures <- list.files("F:/chelsa/cmip5_bio19", full.names=T)

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

      writeRaster(deltas, paste0("../climate_deltas/", i, ".tif"))
}

####### workpoint -- modify code below to incorporate multiple gcms

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





############ mahalanobis #########

# mahalanobis distance to future mean climate
md <- function(x, n, ...){

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
mahal500 <- aggregate(clim, fact=c(500, 500, 4), fun=function(x, ...) md(x, n=500))
mahal <- aggregate(clim, fact=c(50, 50, 4), fun=function(x, ...) md(x, n=50))


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

