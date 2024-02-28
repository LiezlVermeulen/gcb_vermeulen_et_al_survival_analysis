# Load required packages
library(lubridate)
library(gtsummary)
library(dplyr)
library (raster)
library(rgdal)
library(terra)
library(gtools)
library(sf)
library(MASS)
library(ggplot2)
library(tidyr)
library(reshape2)
library(lattice)
library(sampling)
library(jtools)
library(effects)
library(ggeffects)
library(ggplot2)
library(RColorBrewer)
library(tidyverse)
library(spaMM)

# Set-up ------------------------------------------------------------------
# Set directories for datasets: ecosystem chacateristics, climatic drivers and breakpoints
wd.eco <- ("C:/Users/u0142455/Documents/PhD/Processing/ch1/SA/data/ecosystem_char/")
wd.drivers <- ("C:/Users/u0142455/Documents/PhD/Processing/ch1/SA/data/drivers/")
wd.tps <- ("C:/Users/u0142455/Documents/PhD/Processing/ch1/SA/data/tps_sem/")

# Load BPs: 1 - BP, 0 = no BP
tps.coord <- shapefile(file.path(wd.tps,"model_input.shp"))
tps.r <- rast(file.path(wd.tps,"model_input.tif"))
tps.coord$Year = as.numeric(as.character(tps.coord$Year)) 
tps.year <- terra::rasterize(vect(tps.coord), rast(tps.r), field="Year")

# Load ecosystem characteristics
tempWet.r <- rast(file.path(wd.eco, "bioclim", "bio8_reproj.tif"))/10
precDry.r <- rast(file.path(wd.eco, "bioclim", "bio14_reproj.tif"))
precSeas.r <- rast(file.path(wd.eco, "bioclim", "bio15_reproj.tif"))
slope.r <- rast(file.path(wd.eco, "terrain", "slope_reproj.tif"))
soils.r = rast(file.path(wd.eco, "soil", "soil_type_reproj.tif"))
soilpH.r <- rast(file.path(wd.eco, "soil", "soilpH_reproj.tif"))/10
soilSand.r <- rast(file.path(wd.eco, "soil", "soilSand_reproj.tif"))
vegType.r = rast(file.path(wd.eco, "bioregions", "vegType_reproj.tif"))

# Load climatic drivers
# 1) Drought
# Long-term drought severity
drought.long.sev <- rast(file.path(wd.drivers, "drought", "spi/spi_12", "spi12_sev_all_na_final_v2.tif"))
drought.long.sev.subset <- subset(drought.long.sev, 2:30)
names(drought.long.sev.subset) <- c('sev.1', 'sev.2', 'sev.3', 'sev.4', 'sev.5', 'sev.6', 'sev.7', 'sev.8', 'sev.9', 'sev.10',
                              'sev.11', 'sev.12', 'sev.13', 'sev.14', 'sev.15', 'sev.16', 'sev.17', 'sev.18', 'sev.19', 'sev.20',
                              'sev.21', 'sev.22', 'sev.23', 'sev.24', 'sev.25', 'sev.26', 'sev.27', 'sev.28', 'sev.29')

# 2) Wetness
# Long-term drought severity
wet.long.sev <- rast(file.path(wd.drivers, "drought", "spi/spi_12", "spi12_wetness_sev_all_na_final_v2.tif"))
wet.long.sev.subset <- subset(wet.long.sev, 2:30)
names(wet.long.sev.subset) <- c('sev.1', 'sev.2', 'sev.3', 'sev.4', 'sev.5', 'sev.6', 'sev.7', 'sev.8', 'sev.9', 'sev.10',
                                    'sev.11', 'sev.12', 'sev.13', 'sev.14', 'sev.15', 'sev.16', 'sev.17', 'sev.18', 'sev.19', 'sev.20',
                                    'sev.21', 'sev.22', 'sev.23', 'sev.24', 'sev.25', 'sev.26', 'sev.27', 'sev.28', 'sev.29')

# Plot
tps.pts <- raster::rasterToPoints(raster(tps.r), spatial=TRUE)
tps.true <- subset(tps.pts, tps.pts$TP == 1)

r.names <- c("TP", "Year",
             "Temp. of Wettest Quarter (Degrees C)",
             "Prec. of Driest Month (mm)", "Prec Seasonality (mm)",
             "Slope (degrees)",
             "Soil Type","Soil pH", "Soil Sand Content",
             "Veg. Type",
             "Drought Frequency (10 year)", "Drought Severity (10 year)")

r.list <- c(tps.r, tps.year,
            tempWet.r,
            precDry.r, precSeas.r, 
            slope.r, twi.r, 
            soils.r, soilpH.r, soilSand.r,
            vegType.r)
names(r.list) <- r.names

# Plot all layers
plotAll<- lapply(r.list, function(x){
  plot(crop(x, tps.true), main=names(x))
  plot(tps.true, add=TRUE)
})

# Create data matrix
model.stack <- r.list
tps.df <- terra::extract(model.stack, vect(tps.coord), xy=TRUE)
df.names <- c("ID", "TP", "Year", "tempWet",
              "precDry", "precSeas",
              "slope", "twi",
              "soils", "soilpH", "soilSand",
              "vegType",
              "x", "y")


# Filter Vegtype according to only select arid savanna bioregions
tps.df2 <- tps.df
names(tps.df2) <- df.names
tps.df2 <- tps.df2[!(tps.df2$vegType == "NKu" | tps.df2$vegType == "FOz" | tps.df2$vegType == "Gh" | tps.df2$vegType == "Gm" | tps.df2$vegType == "AZi" | tps.df2$vegType == "AZf" | tps.df2$vegType == "AZa" | tps.df2$vegType == "NKb"),]
tps.df2$vegType <- droplevels(tps.df2$vegType)
tps.df2$soils <- droplevels(as.factor(tps.df2$soils))
tps.df2 <- na.omit(tps.df2)

tps.df2$time <- tps.df2$Year - 1981

tps.droughtsev <- terra::extract(drought.long.sev.subset, vect(tps.coord))
tps.droughtsev <- reshape(tps.droughtsev, varying=c('sev.1', 'sev.2', 'sev.3', 'sev.4', 'sev.5', 'sev.6', 'sev.7', 'sev.8', 'sev.9', 'sev.10',
                                                    'sev.11', 'sev.12', 'sev.13', 'sev.14', 'sev.15', 'sev.16', 'sev.17', 'sev.18', 'sev.19', 'sev.20',
                                                    'sev.21', 'sev.22', 'sev.23', 'sev.24', 'sev.25', 'sev.26', 'sev.27', 'sev.28', 'sev.29'), 
                          idvar="ID", direction="long")
df.droughtsev <- na.omit(tps.droughtsev)
df.droughtsev <- abs(df.droughtsev)

tps.wetsev <- terra::extract(wet.long.sev.subset, vect(tps.coord))
tps.wetsev <- reshape(tps.wetsev, varying=c('sev.1', 'sev.2', 'sev.3', 'sev.4', 'sev.5', 'sev.6', 'sev.7', 'sev.8', 'sev.9', 'sev.10',
                                                    'sev.11', 'sev.12', 'sev.13', 'sev.14', 'sev.15', 'sev.16', 'sev.17', 'sev.18', 'sev.19', 'sev.20',
                                                    'sev.21', 'sev.22', 'sev.23', 'sev.24', 'sev.25', 'sev.26', 'sev.27', 'sev.28', 'sev.29'), 
                          idvar="ID", direction="long")
df.wetsev <- na.omit(tps.wetsev)
df.wetsev <- abs(df.wetsev)


# Create Longitudinal Dataset ---------------------------------------------
df.ind <-
  tmerge(data1=tps.df2,
         data2=tps.df2,
         id=ID,
         event=event(time, TP))

df.final <-
  tmerge(data1=df.final,
         data2=df.droughtsev,
         id=ID,
         drought.sev=cumtdc(time, sev))

df.final <-
  tmerge(data1=df.final,
         data2=df.wetsev,
         id=ID,
         wet.sev=cumtdc(time, sev))

df.final <- na.omit(df.final)
df.final$event <- as.numeric(as.character(df.final$event))

# Rescale
df.final$soilpH2 <- scale(df.final$soilpH, center = TRUE, scale = TRUE)
df.final$soilSand2 <- scale(df.final$soilSand, center = TRUE, scale = TRUE)
df.final$tempWet2 <- scale(df.final$tempWet, center = TRUE, scale = TRUE)
df.final$precSeas2 <- scale(df.final$precSeas, center = TRUE, scale = TRUE)
df.final$precDry2 <- scale(df.final$precDry, center = TRUE, scale = TRUE)
df.final$slope2 <- scale(df.final$slope, center = TRUE, scale = TRUE)
df.final$drought.sev2 <- scale(df.final$drought.sev, center = TRUE, scale = TRUE)
df.final$wet.sev2 <- scale(df.final$wet.sev, center = TRUE, scale = TRUE)

# 1.2 Discrete Survival Analysis ----------------------------------------------
# 1.2.1 Model fitting -----------------------------------------------------

df.final %>%
  group_by(tstop) %>%
  summarise(event = sum(event),
            total = n()) %>%
  mutate(hazard = event/total) %>%
  ggplot(aes(x = tstop, y = log(-log(1-hazard)))) +
  xlab("Time (Years)") +
  ylab("log(-log(1 - Hazard))") +
  geom_point() +
  geom_smooth()

df.final %>%
  group_by(tstop, vegType) %>%
  summarise(event = sum(event),
            total = n()) %>%
  mutate(hazard = event/total) %>%
  ggplot(aes(x = tstop, 
             y = log(-log(1-hazard)),
             col = vegType)) +
  geom_point() +
  geom_smooth()

df.final %>%
  mutate(drought.sev = cut_interval(drought.sev, 10, labels = F)) %>%
  group_by(drought.sev) %>%
  summarise(event = sum(event),
            total = n()) %>%
  mutate(hazard = event/total) %>%
  ggplot(aes(x = drought.sev, y = log(-log(1-hazard)))) +
  geom_point() +
  geom_smooth(method = "lm")

# Create adjacency matrix for spatial autocorrelation
df.final <- na.omit(df.final)
xy.df <- distinct(as.data.frame(cbind(df.final$ID, df.final$x, df.final$y)))
names(xy.df) <- c("ID", "x", "y")
coordinates(xy.df) <- ~ x + y
knea <- knearneigh(coordinates(xy.df), longlat = TRUE)
W.nb <- knn2nb(knea, sym=TRUE)
W <- nb2mat(W.nb, style="B", zero.policy=FALSE)
rownames(W) <- unique(df.final.5y$ID)

# Model with NO spatial autocorrelation incorporation
model_full_SA_original <- fitme(event ~ tstop + soilpH2 + soilSand2 + slope2 + precSeas2 + precDry2 + tempWet2 + drought.sev2 + wet.sev2 + (1|vegType/ID),
                                adjMatrix=W,
                                family = binomial(link = "cloglog"),
                                data = df.final.train)


# Model with spatial autocorrelation incorporation
model_full_SA <- fitme(event ~ tstop + soilpH2 + soilSand2 + slope2 + precSeas2 + precDry2 + tempWet2 + drought.sev2 + wet.sev2 + (1|vegType/ID) + adjacency(1 | ID),
                             adjMatrix=W,
                             family = binomial(link = "cloglog"),
                             data = df.final)


summary(model_full_SA)
exp(fixef(model_full_SA))
extractAIC(model_full_SA)
summary(model_full_SA_original)
extractAIC(model_full_SA_original)
anova(model_full_SA_original, model_full_SA)

plot_effects(model_full_SA, focal_var = "soilSand2", ylim=c(0, 0.3))
plot_effects(model_full_SA, focal_var = "precDry2", ylim=c(0, 0.3))
plot_effects(model_full_SA, focal_var = "tempWet2", ylim=c(0, 0.3))
plot_effects(model_full_SA, focal_var = "slope2", ylim=c(0, 0.3))