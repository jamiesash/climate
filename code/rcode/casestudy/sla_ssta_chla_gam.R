source("../functions.R")
# library("dplyr")
# library(IndexNumR)
# library(imputeTS)
# library(mgcViz)

library(terra)

# ------------------------------------------------------------------------------
# Load SLA, CHL and SSTA datasets. 
# With this method I prefer to download the extent and timespan of interest. 

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

# So close. Close enough for now.  
sst = terra::mask(sst, hawaii, inverse = TRUE)
chl = terra::mask(chl, hawaii, inverse = TRUE)
sla = terra::mask(sla, hawaii, inverse = TRUE)

# ------------------------------------------------------------------------------
# Resize to SLA since it is the lowest resolution. 

sst = resample(sst, sla)
chl = resample(chl, sla)
gc() # koa might not like this. 

### I AM HERE. Just got to get the table working. 

test = as.data.frame(sst)
dim(test)

# ------------------------------------------------------------------------------
# Create table
# this is much simpler code but may take a lot of memory
vectorize = function(ras) {
  df = cbind(coordinates(ras), values(ras))
  df = data.frame(df)
  time = getZ(ras)
  time = rep(time, times = nrow(df))
  df$time = time
  colnames(df) = c("lon", "lat", "value", "time")
  df
}

ras2tbl = function(ras, name = "value"){
  ras_l = list()
  for(i in 1:dim(ras)[3]) ras_l[[i]] = ras[[i]]
  ras_df = lapply(ras_l, vectorize)
  rm(ras_l)
  ras_df = do.call(rbind, ras_df)
  ras_df
}

xdf = ras2tbl(x, name = "value")
ydf = ras2tbl(y, name = "value")
zdf = ras2tbl(z,  name = "value")
gc()

# create one data frame
xyz = data.frame(time = xdf$time, 
                 lat  = xdf$lat, 
                 lon  = xdf$lon, 
                 ssta = xdf$value,
                 sla   = sdf$value,
                 chl   = qdf$value
)

rm(xdf, ydf, zdf)
gc() 

# Remove NA vales for just  these variables. 
xyz = xyz[!is.na(xyz$sla),]
xyz = xyz[!is.na(xyz$chl),]
xyz$time = as.Date(xyz$time)
gc()

# ------------------------------------------------------------------------------
# Preparing a data frame for the GAM
# remove effects of clouds
df = xyz

# Converting lat lon to distance
xy   = cbind(df$lat, df$lon)
utms = rgdal::project(xy, "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0") # need to change zones
df$northing = utms[, 1]/1000
df$easting  = utms[, 2]/1000
df$doy = as.numeric(format(df$time, '%j'))

df = df[, c("chl", "sta", "sla", "doy", "easting", "northing")]

# removing outliers. May not be necessary. 
idx = df$slaa < median(df$slaa) + mad(df$slaa)*5
df = df[idx, ]
idx = df$sla > median(df$sla) - mad(df$slaa)*5
df = df[idx, ]
gc()

# ------------------------------------------------------------------------------
# perform GAM 
gam_chl = gam(chl ~ s(sta) ~ s(sla) + s(doy) + s(easting, northing), 
              #method = "REML",
              data = df, 
              family = Gamma(link = "inverse"))

gc()

gam_chl = gam(chl ~ s(sla) + s(doy) + s(easting, northing), 
              #method = "REML",
              data = df, 
              family = Gamma(link = "inverse"))

gc()

gam_ssta = gam(sta ~ sla + s(doy) + s(easting, northing), 
               method = "REML",
               data = df, 
               family = gaussian(link = "identity"))

gc()

# gam_2018 = gam(chl ~ s(doy, by = regions, k=4) + s(doy) #+ s(easting, northing), 
#                method = "REML",
#                data = df, 
#                family = gaussian(link = "identity"))
# ------------------------------------------------------------------------------
