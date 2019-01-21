


# load libraries
library(raster)
library(dplyr)
library(vegan)
library(colormap)
library(alphashape3d)
library(geometry)
library(alphahull)
library(rgeos)
library(gridExtra)
library(grid)

getData <- raster::getData

# working directory
#setwd("E:/climate_classification")

# climate data
s <- getData("worldclim", var="bio", res=10)

# cali border
b <- getData("GADM", country="USA", level=1)
b <- b[b$NAME_1=="California",]
b <- spTransform(b, crs(s))

# munge
r <- mask(crop(s, extent(b)), b)
s <- mask(crop(s, extent(b)), b)
s[12:19] <- log(s[12:19])
d <- as.data.frame(rasterToPoints(r))
climate <- scale(d[,3:21])

# ordinate
pca <- prcomp(climate)$x[,1:3]


# continuous colors
d[,c("red", "green", "blue")] <- pca %>%
      colors3d(trans="ecdf") %>%
      col2rgb() %>%
      t()
d$hex <- rgb(d$red, d$green, d$blue, maxColorValue=255)

ggplot(d, aes(x, y)) + geom_raster(fill=d$hex)


# 3d

ah <- ashape3d(apply(pca[,1:3], 2, scale), alpha=3)
#plot(ah)





# 2d
pca <- prcomp(climate)$x[,1:2]
#tp <- cbind(climate[,1], climate[,12])
#dupe <- duplicated(tp)

ah <- ahull(pca, alpha=3)

s <- ah2sp(ah) # function below
p <- as.data.frame(pca)
names(p) <- c("tmean", "ppt")
coordinates(p) <- c("tmean", "ppt")


dst <- as.vector(gDistance(p, as(s, "SpatialLines"), byid=T))
plot(s); plot(p, add=T, col=dst)
plot(p[ah$ashape.obj$alpha.extremes,], add=T, col="yellow")

f <- cbind(d, pca, dst)

f$onedge <- F; f$onedge[ah$ashape.obj$alpha.extremes] <- T
f$dst_rcl <- f$dst
f$dst_rcl[f$onedge] <- min(f$dst[!f$onedge])

# continuous colors
f$hex <- colors2d(pca) %>%
      col2rgb() %>%
      t() %>%
      rgb(maxColorValue=255)

sd <- fortify(s)

geo <- ggplot() + 
      geom_raster(data=f, 
                  aes(x, y),
                  fill=f$hex) + 
      geom_raster(data=f, 
                  aes(x, y, alpha=dst_rcl^.01), 
                  fill="gray80") +
      viridis::scale_fill_viridis() +
      theme_minimal() +
      theme(legend.position="none",
            axis.text=element_blank(),
            axis.ticks=element_blank(),
            axis.title=element_blank())

clim <- ggplot() + 
      geom_polygon(data=sd, aes(long, lat),
                   fill=NA, color="black") +
      geom_point(data=f, 
                  aes(PC1, PC2),
                  color=f$hex, size=3) + 
      geom_point(data=f, 
                  aes(PC1, PC2, alpha=dst_rcl^.01), 
                  color="gray80", size=4) +
      viridis::scale_fill_viridis() +
      labs(x="PC1", y="PC2") +
      theme_minimal() +
      theme(legend.position="none",
            axis.text=element_blank(),
            axis.ticks=element_blank())

p <- arrangeGrob(textGrob(label="the edges of extant climate space"), 
                 arrangeGrob(geo, clim, nrow=1), 
                 ncol=1, heights=c(1, 10))
pdf("E:/edges/edge_of_existence.pdf", width=10, height=6)
grid.draw(p)
dev.off()


#####

# http://r-sig-geo.2731867.n2.nabble.com/alpha-hull-ahull-to-polygon-shapefile-td7342734.html

ah2sp <- function(x, increment=360, rnd=10, proj4string=CRS(as.character(NA))){ 
      require(alphahull) 
      require(maptools) 
      if (class(x) != "ahull"){ 
            stop("x needs to be an ahull class object") 
      } 
      # Extract the edges from the ahull object as a dataframe 
      xdf <- as.data.frame(x$arcs) 
      # Remove all cases where the coordinates are all the same       
      xdf <- subset(xdf,xdf$r > 0) 
      res <- NULL 
      if (nrow(xdf) > 0){ 
            # Convert each arc to a line segment 
            linesj <- list() 
            prevx<-NULL 
            prevy<-NULL 
            j<-1 
            for(i in 1:nrow(xdf)){ 
                  rowi <- xdf[i,] 
                  v <- c(rowi$v.x, rowi$v.y) 
                  theta <- rowi$theta 
                  r <- rowi$r 
                  cc <- c(rowi$c1, rowi$c2) 
                  # Arcs need to be redefined as strings of points. Work out the number of points to allocate in this arc segment. 
                  ipoints <- 2 + round(increment * (rowi$theta / 2),0) 
                  # Calculate coordinates from arc() description for ipoints along the arc. 
                  angles <- anglesArc(v, theta) 
                  seqang <- seq(angles[1], angles[2], length = ipoints) 
                  x <- round(cc[1] + r * cos(seqang),rnd) 
                  y <- round(cc[2] + r * sin(seqang),rnd) 
                  # Check for line segments that should be joined up and combine their coordinates 
                  if (is.null(prevx)){ 
                        prevx<-x 
                        prevy<-y 
                  } else if (x[1] == round(prevx[length(prevx)],rnd) && y[1] == round(prevy[length(prevy)],rnd)){ 
                        if (i == nrow(xdf)){ 
                              #We have got to the end of the dataset 
                              prevx<-append(prevx,x[2:ipoints]) 
                              prevy<-append(prevy,y[2:ipoints]) 
                              prevx[length(prevx)]<-prevx[1] 
                              prevy[length(prevy)]<-prevy[1] 
                              coordsj<-cbind(prevx,prevy) 
                              colnames(coordsj)<-NULL 
                              # Build as Line and then Lines class 
                              linej <- Line(coordsj) 
                              linesj[[j]] <- Lines(linej, ID = as.character(j)) 
                        } else { 
                              prevx<-append(prevx,x[2:ipoints]) 
                              prevy<-append(prevy,y[2:ipoints]) 
                        } 
                  } else { 
                        # We have got to the end of a set of lines, and there are several such sets, so convert the whole of this one to a line segment and reset. 
                        prevx[length(prevx)]<-prevx[1] 
                        prevy[length(prevy)]<-prevy[1] 
                        coordsj<-cbind(prevx,prevy) 
                        colnames(coordsj)<-NULL 
                        # Build as Line and then Lines class 
                        linej <- Line(coordsj) 
                        linesj[[j]] <- Lines(linej, ID = as.character(j)) 
                        j<-j+1 
                        prevx<-NULL 
                        prevy<-NULL 
                  } 
            } 
            # Promote to SpatialLines 
            lspl <- SpatialLines(linesj) 
            # Convert lines to polygons 
            # Pull out Lines slot and check which lines have start and end points that are the same 
            lns <- slot(lspl, "lines") 
            polys <- sapply(lns, function(x) { 
                  crds <- slot(slot(x, "Lines")[[1]], "coords") 
                  identical(crds[1, ], crds[nrow(crds), ]) 
            }) 
            # Select those that do and convert to SpatialPolygons 
            polyssl <- lspl[polys] 
            list_of_Lines <- slot(polyssl, "lines") 
            sppolys <- SpatialPolygons(list(Polygons(lapply(list_of_Lines, function(x) { Polygon(slot(slot(x, "Lines")[[1]], "coords")) }), ID = "1")), proj4string=proj4string) 
            # Create a set of ids in a dataframe, then promote to SpatialPolygonsDataFrame 
            hid <- sapply(slot(sppolys, "polygons"), function(x) slot(x, "ID")) 
            areas <- sapply(slot(sppolys, "polygons"), function(x) slot(x, "area")) 
            df <- data.frame(hid,areas) 
            names(df) <- c("HID","Area") 
            rownames(df) <- df$HID 
            res <- SpatialPolygonsDataFrame(sppolys, data=df) 
            res <- res[which(res@data$Area > 0),] 
      }   
      return(res) 
} 


######################## 

#So for example: 
library(alphahull) 
library(maptools) 
n <- 300 
theta <- runif(n,0,2*pi) 
r <- sqrt(runif(n,0.25^2,0.5^2)) 
x <- cbind(0.5+r*cos(theta),0.5+r*sin(theta)) 
a <- ahull(x, alpha = 0.1) 
plot(a) 

#Convert 
b <- ah2sp(a) 
plot(b, col="lightgrey", pbg="white") 


