# library("dplyr")
# library(IndexNumR)
# library(imputeTS)
library(mgcViz)
library(lubridate)
library(terra)

# ------------------------------------------------------------------------------
# Load SLA, CHL and SSTA data sets. 
# With this method I prefer to download the extent and time span of interest. 

sla = rast("../../../data/sla/sla_2018_l4_4k.nc")
chl = rast("../../../data/chl/chl_2018_daily_multi_l3_4k.nc")
sst = rast("../../../data/sst/ssta_2018_l4_4k.nc")

# subset spacial polygone like a data frame. 
eez = terra::vect("../../../data/eez/USMaritimeLimitsNBoundaries.shp")
idx = which(eez$REGION == 'Hawaiian Islands')
hawaii = eez[idx,]

# ------------------------------------------------------------------------------
# removing coastal effects

# Select the regions of interest for subseting. 
idx = hawaii$CZ + hawaii$TS
idx = which(idx == 1)
hawaii = hawaii[idx,]

# ------------------------------------------------------------------------------
# Resize to SLA since it is the lowest resolution. 

sst = resample(sst, sla)
chl = resample(chl, sla)
sla = resample(sla, sla)

# So close. Close enough for now.  
sst = terra::mask(sst, hawaii, inverse = TRUE)
chl = terra::mask(chl, hawaii, inverse = TRUE)
sla = terra::mask(sla, hawaii, inverse = TRUE)

# ------------------------------------------------------------------------------
# Create table

# welp this used to be really hard before terra. 
chl = terra::as.data.frame(chl, xy = TRUE, wide = FALSE, na.rm=FALSE, time=TRUE)
sla = terra::as.data.frame(sla, xy = TRUE, wide = FALSE, na.rm=FALSE, time=TRUE)
sst = terra::as.data.frame(sst, xy = TRUE, wide = FALSE, na.rm=FALSE, time=TRUE)

# ------------------------------------------------------------------------------
# remove clouds. 

# this seems not to work.
idx = which(!is.na(chl$values))
chl = chl[idx,]
sla = sla[idx,]
sst = sst[idx,]

time = sla$time
time = yday(time)
alldata = data.frame(time = time, 
                     lon = sla$x, 
                     lat = sla$y, 
                     sla = sla$values, 
                     chl = chl$values,
                     sst = sst$values)
# KOA might not like this. 
rm(sla, chl, sst, time)

# ------------------------------------------------------------------------------
# going for it. 

gam_chl = gam(chl ~ sst + sla + s(time) + s(lon, lat), 
              data = alldata,
              family = Gamma(link = "inverse"))

# ------------------------------------------------------------------------------
# plotting the gam output. 

gam_visual <- getViz(gam_chl)
gridPrint(plot(pterm(gam_visual, 1)) + l_ciPoly() + l_fitLine(),
          plot(pterm(gam_visual, 2)) + l_ciPoly() + l_fitLine(),
          ncol = 2)

# ------------------------------------------------------------------------------
# Preparing a data frame for the GAM
# remove effects of clouds
# xyz = xyz[!is.na(xyz$sla),]
# xyz$time = as.Date(xyz$time)
# gc()

# # Converting lat lon to distance
# df = xyz
# xy   = cbind(df$lat, df$lon)
# utms = rgdal::project(xy, "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0") # need to change zones
# df$northing = utms[, 1]/1000
# df$easting  = utms[, 2]/1000
# df$doy = as.numeric(format(df$time, '%j'))
# 
# df = df[, c("chl", "sta", "sla", "doy", "easting", "northing")]
# 
# # removing outliers. May not be necessary. 
# idx = df$slaa < median(df$slaa) + mad(df$slaa)*5
# df = df[idx, ]
# idx = df$sla > median(df$sla) - mad(df$slaa)*5
# df = df[idx, ]
# gc()

# # ------------------------------------------------------------------------------
# # perform GAM 
# gam_chl = gam(chl ~ s(sta) ~ s(sla) + s(doy) + s(easting, northing), 
#               #method = "REML",
#               data = df, 
#               family = Gamma(link = "inverse"))
# 
# gam_chl = gam(chl ~ s(sla) + s(doy) + s(easting, northing), 
#               #method = "REML",
#               data = df, 
#               family = Gamma(link = "inverse"))
# 
# gc()
# 
# gam_ssta = gam(sta ~ sla + s(doy) + s(easting, northing), 
#                method = "REML",
#                data = df, 
#                family = gaussian(link = "identity"))
# 
# gc()

# gam_2018 = gam(chl ~ s(doy, by = regions, k=4) + s(doy) #+ s(easting, northing), 
#                method = "REML",
#                data = df, 
#                family = gaussian(link = "identity"))
# ------------------------------------------------------------------------------
