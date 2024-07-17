library(yaml)
library(mgcViz)
library(lubridate)
library(terra)

# ------------------------------------------------------------------------------
### To do

# 1. Download SSTA for the same resolution as CHL. 
# 2. Resize to CHL for higher fake resolution. Use nearest neighbor. 
# 3. Run the same GAM as before. 
# 4. Save the model output. 
# 5. Re-download SST I deleted it by accident. 
# 6. SLA resize introduces extreme negative number.  
# 7. Fix the way the way SSTA is being saved. It's extent is messed up. 
# 8. Download a different CHL this one is cloudy. 

# ------------------------------------------------------------------------------
# Load SLA, CHL and SSTA data sets. 
# With this method I prefer to download the extent and time span of interest. 

setwd("/home/jamie/projects/climate/code/rcode/casestudy")

sla = rast("../../../data/sla/sla_2018_l4_4k.nc")
chl = rast("../../../data/chl/chl_2018_daily_multi_l3_4k.nc")
sst = rast("../../../data/sst/ssta_l4_2018_4k_20240717.nc")

# subset spacial polygon like a data frame. 
eez = terra::vect("../../../data/eez/USMaritimeLimitsNBoundaries.shp")
idx = which(eez$REGION == 'Hawaiian Islands')
hawaii = eez[idx,]

# ------------------------------------------------------------------------------
# removing coastal effects

idx = hawaii$CZ + hawaii$TS
idx = which(idx == 1)
hawaii = hawaii[idx,]

# ------------------------------------------------------------------------------
# Resize to SLA since it is the lowest resolution. 

sst = resample(sst, sla, method = "near")
chl = resample(chl, sla, method = "near")
#sla = resample(sla, sla, method = "near")

# So close. Close enough for now. Still leaves little blips. 
sst = terra::mask(sst, hawaii, inverse = TRUE)
chl = terra::mask(chl, hawaii, inverse = TRUE)
sla = terra::mask(sla, hawaii, inverse = TRUE)

# ------------------------------------------------------------------------------
# Create table

# Welp this used to be really hard before terra. 
chl = as.data.frame(chl, xy = TRUE, wide = FALSE, na.rm=FALSE, time=TRUE)
sla = as.data.frame(sla, xy = TRUE, wide = FALSE, na.rm=FALSE, time=TRUE)
sst = as.data.frame(sst, xy = TRUE, wide = FALSE, na.rm=FALSE, time=TRUE)

# ------------------------------------------------------------------------------
# remove clouds. 

# This seems not to work. 
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

write_yaml(gam_chl, "chl_gam_20240715.yml")

# ------------------------------------------------------------------------------
# plotting the gam output. 

gam_visual <- getViz(gam_chl)
gridPrint(plot(pterm(gam_visual, 1)) + l_ciPoly() + l_fitLine(),
          plot(pterm(gam_visual, 2)) + l_ciPoly() + l_fitLine(),
          ncol = 2)

# save the model output. 

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
