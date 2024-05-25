source("functions.R")
library("dplyr")
library(IndexNumR)
library(imputeTS)
library(terra)
library(rgdal)
library(mgcViz)

# ------------------------------------------------------------------------------
## the data loading function we've all been waiting for
#dap = function(sdate, edate, e, 
#               id = 'jplMURSST41anom1day', 
#               url = "https://coastwatch.pfeg.noaa.gov/erddap/",
#               var = c("sstAnom", "latitude", "longitude", "time")) {
#  # input is a datetime and an extent. It grabs that raster
#  dap_that_ass = function(x, data_info, e){
#    data = griddap(data_info, 
#                   latitude = e[3:4], 
#                   longitude = e[1:2], 
#                   time = c(as.character(x),as.character(x)), 
#                   fields = 'all')
#    data = nc_open(data$summary$filename)
#    ras  = ncvar_get(data, varid = var[1])
#    lats = ncvar_get(data, varid = var[2])
#    lons = ncvar_get(data, varid = var[3])
#    time = ncvar_get(data, varid = var[4])
#    time = as.Date(as.POSIXct(time, origin = "1970-01-01"))
#    nc_close(data)
#    rm(data)
#    
#    ras = raster(ras)
#    extent(ras) = extent(min(lons), max(lons), min(lats), max(lats))
#    ras = setZ(ras, z = time, name = "time")
#    crs(ras) = "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"
#    ras
#  }
#  
#  time = seq(sdate, edate, by = 1)
#  time = as.list(time)
#  
#  data_info = info(id, 
#                   url = url)
#  
#  ras  = lapply(time, FUN = dap_that_ass, data_info = data_info, e = e)
#  time = lapply(ras, getZ)
#  time = unlist(time)
#  time = as.Date(time)
#  ext    = lapply(ras, extent)
#  ras  = stack(ras)
#  extent(ras) = ext[[1]]
#  ras = setZ(ras, z= time, name = "time")
#  ras
#}
#

# ------------------------------------------------------------------------------
# Set some variables
lons = c(-175, -125)
lats = c(18,   40)
e = extent(lons, lats)
sdate = as.Date("2018-07-05")
edate = as.Date("2018-10-15")

# ------------------------------------------------------------------------------
# Load SLA, CHL and SSTA datasets. 
# With this method I prefer to download the extent and timespan of interest. 

sla = rast("../../data/ssta_2018_2022_l4_4k.nc")
chl = rast("../../data/ssta_2018_2022_l4_4k.nc")
sta = rast("../../data/ssta_2018_2022_l4_4k.nc")

# ------------------------------------------------------------------------------
# Make them the same time domain. Make FSLE same temporal resolution as CHL
# This drops missing values and such. 
match = function(x, y){
  ty = getZ(y)
  ty = as.numeric(ty)
  tx = getZ(x)
  tx = as.numeric(tx)
  idx = which(is.element(tx, ty))
  tx = subset(tx, is.element(tx, ty))
  x = raster::subset(x, idx)
  x = setZ(x, z = as.Date(tx), name = "time")
  x
}

sla = match(x = sla, y = chl)
sta = match(x = sta, y = chl)
gc()

# ------------------------------------------------------------------------------
# Crop all to the same extent and time range. 

# Everything after here used to be the table function. 

resize(ras, ras2){
    ras = crops(ras, ras2)
    tmp = resample(ras, ras2)
    ras["values"] = tmp["values"]
    ras
}

sla = resize(sla, chl)
sta = resize(sta, chl)
gc()

# ------------------------------------------------------------------------------
# removing coastal effects
sta = bufcoast(sta)
chl = bufcoast(chl)
sla = bufcoast(sla)
gc()

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

# removing outliers
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
