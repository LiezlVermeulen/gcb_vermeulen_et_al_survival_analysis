## ----------------------------------------------------
## Script for calculating the annual drought and wetness severity per pixel
## ----------------------------------------------------
# Load required packages
library(SPEI)
library (raster)
library(rgdal)
library (rts)
library(terra)
library(gtools)
library(rasterVis)
library(viridis)
library(sp)
library(RColorBrewer)

wd <- 'C:/Users/u0142455/Documents/PhD/Processing/ch1/SA/data/chirps'
setwd(wd)

# Load datasets 
aoi <- shapefile("C:/Users/u0142455/Documents/PhD/Processing/ch1/SA/data/aoi/savanna_studyarea.shp") # study area
in.files <- list.files('D:/chirps/tif', "tif$", full.names = TRUE) # CHIRPS rainfall data
chirps.stack <- stack(in.files)
chirps.reprj <- projectRaster(chirps.stack, TPs.savanna, method="bilinear")

# SPI
# stack years on chirps
chirps.mask <- mask(chirps.reprj, aoi)
spi.brick <- brick(chirps.mask )

# Calculate annual drought severity
xspisev <- function(data) { 
  data <- as.vector(data[1:420])
  chirps.ts <- ts(data, start=1981, end=c(2015, 12), frequency=12)
  spi.12 <- spi(chirps.ts, scale=12, na.rm=TRUE)
  spi.12.vals <- as.vector(as.numeric(spi.12$fitted))
  spi.12.vals[spi.12.vals > -1] <- 0 # values higher than drought threshold
  spi.12.vals[is.na(spi.12.vals)] <- 0
  spi.12.ts <- as.zoo(ts(spi.12.vals, frequency=12, start=c(1981, 1), end=c(2015, 12))) # Change years according to temproal period
  
  # Yearly sum
  as.year <- function(x) as.numeric(floor(as.yearmon(x)))
  drought.sev <- as.ts(aggregate(spi.12.ts, as.year, sum))
  drought.sev.v <- as.numeric(drought.sev)
  
  return(drought.sev.v)
}

spi.sev.output <- calc(spi.brick, fun=xspisev)
plot(spi.sev.output)
sum.spi.sev <- sum(spi.sev.output)

# Plot the drought severity for certain years
r.names <- as.character(seq(1981, 2015, by=1))
names(spi.sev.output) <- r.names
plot(subset(abs(spi.sev.output), 18:20))

plot.subset <- crop(subset(abs(spi.sev.output), 17:22), north)
colr <- colorRampPalette(brewer.pal(9, 'YlOrRd'))
p1 <- levelplot(plot.subset, 
                margin=FALSE, 
                colorkey=list(
                  space='bottom',                   
                  labels=list(at=seq(0, 12, 3), font=4),
                  axis.line=list(col='black'),
                  width=0.75
                ), 
                par.settings=list(
                  strip.border=list(col='transparent'),
                  strip.background=list(col='transparent'),
                  axis.line = list(col = "black")
                ),
                scales=list(draw=FALSE),            
                col.regions=colr,                   
                at=c(-Inf, seq(0, 12), Inf),
                names.attr=c( "1997", "1998", "1999","2000", "2001", "2002"))
p1 <- p1 + latticeExtra::layer(sp.polygons(aoi))
p1

# Calculate annual wetness severity
xspiwetsev <- function(data) { 
  data=as.vector(data[1:420])
  chirps.ts = ts(data, start=1981, end=c(2015, 12), frequency=12)
  spi.12 <- spi(chirps.ts, scale=12, na.rm=TRUE)
  spi.12.vals <- as.vector(as.numeric(spi.12$fitted))
  spi.12.vals[spi.12.vals < 1] <- 0 # values lower than extreme wet threshold
  spi.12.vals[is.na(spi.12.vals)] <- 0
  spi.12.ts <- as.zoo(ts(spi.12.vals, frequency=12, start=c(1981, 1), end=c(2015, 12)))
  
  # Yearly sum
  as.year <- function(x) as.numeric(floor(as.yearmon(x)))
  wetness.sev <- as.ts(aggregate(spi.12.ts, as.year, sum))
  wetness.sev.v <- as.numeric(wetness.sev)
  
  return(wetness.sev.v)
}

spi.wet.sev.output <- calc(spi.brick, fun=xspiwetsev)
spi.wet.sev.output
plot(spi.wet.output)

# Plot the drought severity for certain years
r.names <- as.character(seq(1981, 2015, by=1))
names(spi.wet.output.mask) <- r.names

plot.subset <- crop(subset(spi.wet.output.mask, 15:20), clust.b)
colr <- colorRampPalette(brewer.pal(9, 'YlGnBu'))
p1 <- levelplot(plot.subset, 
          margin=FALSE, 
          colorkey=list(
            space='bottom',                   
            labels=list(at=seq(0, 25, 5), font=4),
            axis.line=list(col='black'),
            width=0.75
          ), 
          par.settings=list(
            strip.border=list(col='transparent'),
            strip.background=list(col='transparent'),
            axis.line = list(col = "black")
          ),
          scales=list(draw=FALSE),            
          col.regions=colr,                   
          at=c(-Inf, seq(0, 25), Inf),
          names.attr=c( "1995", "1996", "1997", "1998", "1999","2000"))
p1 <- p1 + latticeExtra::layer(sp.polygons(aoi))
p1