
source("functions.R")
source("libraries.R")
source("blooms.R")
library(animation)
library("rerddap")
library("akima")
library("dplyr")
library("ggplot2")
library("mapdata")
library("ncdf4")
library("plot3D")
library(IndexNumR)
library(imputeTS)
library(sp)
#library(rgdal)
library(raster)
library("plotrix") 
rasterOptions(maxmemory = 120e+10, memfrac = 0.9)

# ------------------------------------------------------------------------------

lons = c(-170, -125)
lats = c( 17,   36)
e = extent(lons, lats)
sdate = as.Date("2018-06-18")
edate = as.Date("2018-11-15")

# ------------------------------------------------------------------------------
# loading chl data
# nceiPH53sstd1day
chlInfo = info('noaacwNPPN20S3ASCIDINEOFDaily',
               url = "https://coastwatch.noaa.gov/erddap/")

data = griddap(chlInfo, 
               latitude = lats, 
               longitude = lons, 
               time = c(as.character(sdate),as.character(edate)), 
               fields = 'all')

data = nc_open(data$summary$filename)
ras  = ncvar_get(data)
lats = ncvar_get(data, varid = "latitude")
lons = ncvar_get(data, varid = "longitude")
time = ncvar_get(data, varid = "time")
time = as.Date(as.POSIXct(time, origin = "1970-01-01"))

nc_close(data)
rm(data)

ras = brick(ras) 
ras = setZ(ras, z = time, name = "time")
extent(ras) = extent(min(lons), max(lons), min(lats), max(lats))
crs(ras) = "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"
ras  = oreant(ras, t1 = TRUE)

chl = ras

# ------------------------------------------------------------------------------
# loading fsle
# load fsle from disk
fsle = load_nc(path = "C:\\Users\\james\\Desktop\\jamieslife\\data\\infiles\\case_study\\fsle\\",
               patt = "fsle_po_2018.nc")
extent(fsle) = extent(fsle) - c(360, 360, 0 ,0)
gc()

# ------------------------------------------------------------------------------
# load sla from disk
sla = load_nc(path = "C:\\Users\\james\\Desktop\\jamieslife\\data\\infiles\\case_study\\sla\\",
              patt = "2010_2021_sla.nc",
              vars = c("latitude", "longitude", "sla", "time"))

# ------------------------------------------------------------------------------
# fix oreantation
fsle = oreant(fsle, flip = "y", t1 = TRUE)
gc()
sla  = oreant(sla, flip = "y", t1 = TRUE)
gc()

# ------------------------------------------------------------------------------

# Crop all to the same extent

fsle = raster::crop(fsle, e)
sla  = raster::crop(sla,  e)
chl  = raster::crop(chl,  e)
gc()


# ------------------------------------------------------------------------------
# Calculate a sla anomaly
# slaa = anomalize(ras = sla, detrend = FALSE)
slaa = anom(sla)
rm(sla)
gc()

# ------------------------------------------------------------------------------

# Cut to timespan of bloom from table 1
# Not completly nessisary now will see for latter
chl  = timesnip(chl,   sdate, edate)
fsle = timesnip(fsle,  sdate, edate)
slaa = timesnip(slaa,  sdate, edate)
gc()

# ------------------------------------------------------------------------------
# Remove coastal influence: SOmething is wrong here

chl  = bufcoast(chl, 
                region = "Hawaiian Islands", 
                path = "../data/USMaritimeLimitsAndBoundariesSHP")
fsle = bufcoast(fsle, 
                region = "Hawaiian Islands", 
                path = "../data/USMaritimeLimitsAndBoundariesSHP")
slaa = bufcoast(slaa, 
                region = "Hawaiian Islands", 
                path = "../data/USMaritimeLimitsAndBoundariesSHP")
gc()

# ------------------------------------------------------------------------------

# Match FSLE to SLAA because SLAA lost leap days. Make FSLE and CHL them the same 
# time domain. Make FSLE same temporal resolution as CHL. Not needed for daily 

match = function(x, y){
  ty = getZ(y)
  ty = as.numeric(ty)
  tx = getZ(x)
  tx = as.numeric(tx)
  
  idx  = which(is.element(tx, ty))
  tx   = subset(tx, is.element(tx, ty))
  x = raster::subset(x, idx)
  x = setZ(x, z = as.Date(tx), name = "time")
  x
}

fsle = match(x = fsle, y = chl)
chl  = match(x = chl, y = fsle)
slaa = match(x = slaa, y = fsle)

gc()

# ------------------------------------------------------------------------------

# Resample to the size of fsle

resample = function(ras, to, method = "bilinear"){
  raster::crs(ras) = "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0" 
  time = getZ(ras)
  e = extent(to)
  extent(ras) = extent(e)
  ras = raster::resample(ras,  to, method = method)
  extent(ras) = e
  ras = setZ(ras, z = as.Date(time), name = "time") 
  ras
}

# this can be nearest neighbor if needed
slaa  = resample(slaa,  to = fsle, method = "bilinear")
gc()
# this can be nearest neighbor if needed
fsle  = resample(fsle,  to = fsle, method = "bilinear")
gc()

# ------------------------------------------------------------------------------

# this might not be necissary

fronts = function(ras,
                  downsize = 10,
                  I = 2){
  e    = extent(ras)
  ras  = t(ras)
  ras  = raster::flip(ras, direction = "x")
  extent(ras) = e
  
  thresh = median(ras, na.rm = TRUE) + mad(ras, na.rm = TRUE)*I
  idx  = which(values(ras) > thresh)
  
  e    = extent(ras)
  s    = dim(ras)
  lons = seq(from = e[1], to = e[2], length = s[1]) 
  lats = seq(from = e[3], to = e[4], length = s[2]) 
  grid = pracma::meshgrid(lons, lats)
  fronts = data.frame(lons = grid$X[idx], lats = grid$Y[idx], value = ras[idx])
  idx = seq(1, nrow(fronts), downsize)
  fronts = fronts[idx, ]
  fronts
}

chl_l  = list()
for(i in 1:dim(chl)[3]) chl_l[[i]] = chl[[i]]

chl_points  = lapply(chl_l, fronts, I = 2, downsize = 1)

pretty = function(x, thresh = 0.05, clmp = 0.15) {
  x = subset(x, !is.na(value))
  x = subset(x, value > thresh)
  
  idx = which(x$value > clmp)
  # if there are no points it gives an error
  if(length(idx)>1) x[idx,]$value = clmp
  
  x
}

chl_points = lapply(chl_points, FUN = pretty, thresh = 0.05, clmp = 0.25)

# ------------------------------------------------------------------------------
# Create table
# this is much simpler code but may take a lot of memory
vectorize = function (ras) {
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

slaa_df = ras2tbl(slaa, name = "value")
gc()
fsle_df = ras2tbl(fsle, name = "value")
gc()

# create one data frame
xyz = data.frame(time = slaa_df$time, 
                 lat  = slaa_df$lat, 
                 lon  = slaa_df$lon, 
                 fsle = fsle_df$value,
                 slaa = slaa_df$value)

rm(slaa_df, fsle_df)
gc() 

xyz = xyz[!is.na(xyz$slaa),]
xyz = xyz[!is.na(xyz$fsle),]
xyz$fsle = xyz$fsle * -1

gc()

# ------------------------------------------------------------------------------
# saving the data to harddisck
write.csv(xyz, file = paste("data/collated/xyz", Sys.Date(), ".csv", sep = ""))







