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
library(rasterVis)
rasterOptions(maxmemory = 120e+10, memfrac = 0.9)
knitr::opts_chunk$set(echo = TRUE)

# ------------------------------------------------------------------------------
# the data loading function we've all been waiting for
dap = function(sdate, edate, e, 
               id = 'jplMURSST41anom1day', 
               url = "https://coastwatch.pfeg.noaa.gov/erddap/",
               var = c("sstAnom", "latitude", "longitude", "time")) {
  # input is a datetime and an extent. It grabs that raster
  dap_that_ass = function(x, data_info, e){
    data = griddap(data_info, 
                   latitude = e[3:4], 
                   longitude = e[1:2], 
                   time = c(as.character(x),as.character(x)), 
                   fields = 'all')
    data = nc_open(data$summary$filename)
    ras  = ncvar_get(data, varid = var[1])
    lats = ncvar_get(data, varid = var[2])
    lons = ncvar_get(data, varid = var[3])
    time = ncvar_get(data, varid = var[4])
    time = as.Date(as.POSIXct(time, origin = "1970-01-01"))
    nc_close(data)
    rm(data)
    
    ras = raster(ras)
    extent(ras) = extent(min(lons), max(lons), min(lats), max(lats))
    ras = setZ(ras, z = time, name = "time")
    crs(ras) = "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"
    ras
  }
  
  time = seq(sdate, edate, by = 1)
  time = as.list(time)
  
  data_info = info(id, 
                   url = url)
  
  ras  = lapply(time, FUN = dap_that_ass, data_info = data_info, e = e)
  time = lapply(ras, getZ)
  time = unlist(time)
  time = as.Date(time)
  ext    = lapply(ras, extent)
  ras  = stack(ras)
  extent(ras) = ext[[1]]
  ras = setZ(ras, z= time, name = "time")
  ras
}

# ------------------------------------------------------------------------------
# load ssta
lons = c(-175, -125)
lats = c(18,   40)
e = extent(lons, lats)
sdate = as.Date("2018-07-05")
edate = as.Date("2018-10-15")

ssta = dap(sdate = sdate, 
           edate = edate, 
           e = e,
           id = 'jplMURSST41anom1day',
           var = c("sstAnom", "latitude", "longitude", "time"))
ssta = oreant(ssta, flip = "y", t1 = TRUE)

# ------------------------------------------------------------------------------
# now for chl
sdate = as.Date("2018-06-01")
edate = as.Date("2022-12-30")

chl = dap(sdate = sdate, 
          edate = edate, 
          e = e,
          var = c("chlor_a", "latitude", "longitude", "time"),
          id = 'nesdisVHNnoaaSNPPnoaa20chlaGapfilledDaily')
chl = oreant(chl, t1 = TRUE)
gc()

# ------------------------------------------------------------------------------
# now for sla
sdate = as.Date("2018-06-01")
edate = as.Date("2022-12-30")

sla = dap(sdate = sdate, 
          edate = edate, 
          e = e,
          id = 'nesdisSSH1day',
          url = "https://coastwatch.pfeg.noaa.gov/erddap/",
          var = c("sla", "latitude", "longitude", "time"))

gc()

# ------------------------------------------------------------------------------
# fix sla oreantation
sla  = oreant(sla, flip = "y", t1 = TRUE)
gc()

# ------------------------------------------------------------------------------
# Crop all to the same extent. Not nessisary?
e = extent(ssta)
sla  = raster::crop(sla,  e)
gc()

# ------------------------------------------------------------------------------
# Calculate CHL anomaly and SLA anomaly
# I receive an warning when detrending for cells of all NA values (hawaii)
# chl anom. calc. needs at least one full year to work correctly
chla = anomalize(chl, detrend = FALSE)
slaa = anomalize(sla, detrend = FALSE)
gc(verbose = FALSE, full = TRUE)

# ------------------------------------------------------------------------------
# probs not needed
# by time
sdate = as.Date("2018-07-05")
edate = as.Date("2018-10-15")
chl  = timesnip(chl,  sdate, edate)
chla = timesnip(chla,  sdate, edate)
slaa = timesnip(slaa,  sdate, edate)
sla  = timesnip(sla,  sdate, edate)
ssta  = timesnip(ssta,  sdate, edate)
gc()

# ------------------------------------------------------------------------------
# Make them the same time domain. Make FSLE same temporal resolution as CHL
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

slaa = match(x = slaa, y = chla)
sla = match(x = sla, y = chla)
gc()

# ------------------------------------------------------------------------------
# removing coastal effects
ssta = bufcoast(ssta)
gc()
chla = bufcoast(chla)
chl = bufcoast(chl)
gc()
slaa = bufcoast(slaa)
sla = bufcoast(sla)
gc()

# ------------------------------------------------------------------------------
# cutting down the ssta data for statistics not mapping
# by extent
lons = c(-155, -130)
lats = c(22,   34)
e = extent(lons, lats)
x = raster::crop(ssta, e)
y = raster::crop(chla, e)
z = raster::crop(slaa, e)
s = raster::crop(sla, e)
q = raster::crop(chl, e)

# by time
sdate = as.Date("2018-07-01")
edate = as.Date("2018-10-07")
x = timesnip(x,  sdate, edate)
y = timesnip(y,  sdate, edate)
z = timesnip(z,  sdate, edate)
s = timesnip(s,  sdate, edate)
q = timesnip(q,  sdate, edate)
gc()

# ------------------------------------------------------------------------------
# making same size for statistics
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

z = resample(z, to = y, method = "ngb")
gc()

x = resample(x, to = y, method = "ngb")
gc()

q = resample(q, to = y, method = "ngb")
gc()

s = resample(s, to = y, method = "ngb")
gc()

extent(x) = extent(y)
extent(z) = extent(y)
extent(q) = extent(y)
extent(s) = extent(y)

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

xdf = ras2tbl(x, name = "value")
gc()

ydf = ras2tbl(y, name = "value")
gc()

zdf  = ras2tbl(z,  name = "value")
gc()

qdf  = ras2tbl(q,  name = "value")
gc()

sdf  = ras2tbl(s,  name = "value")
gc()

# create one data frame
xyz = data.frame(time = xdf$time, 
                 lat  = xdf$lat, 
                 lon  = xdf$lon, 
                 ssta = xdf$value,
                 slaa  = zdf$value,
                 sla   = sdf$value,
                 chla  = ydf$value,
                 chl   = qdf$value
)

rm(xdf, ydf, zdf, qdf, sdf)
gc() 

xyz = xyz[!is.na(xyz$sla),]
gc()

# ------------------------------------------------------------------------------
# Preparing a data frame for the GAM
# remove effects of clouds
df = xyz

# Cleaning up the data set
df  = xyz[!is.na(df$chla),]
df$time = as.Date(df$time)
gc()

# Converting lat lon to distance
xy   = cbind(df$lat, df$lon)
utms = rgdal::project(xy, "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0") # need to change zones
df$northing = utms[, 1]/1000
df$easting  = utms[, 2]/1000
df$doy = as.numeric(format(df$time, '%j'))

df = df[, c("chla", "chl", "ssta", "sla", "slaa", "doy", "easting", "northing")]

# cutting the end date back
# idx = df$doy < 309
# df = df[idx, ]

# removing outliers

idx = df$slaa < median(df$slaa) + mad(df$slaa)*5
df = df[idx, ]
idx = df$sla > median(df$sla) - mad(df$slaa)*5
df = df[idx, ]
gc()

# ------------------------------------------------------------------------------
# perform GAM 
library(rgdal)
library(mgcViz)

# Individualy and not anomaly 

gam_chl = gam(chl ~ s(slaa) + s(doy) + s(easting, northing), 
              #method = "REML",
              data = df, 
              family = Gamma(link = "inverse"))

gc()

gam_ssta = gam(ssta ~ s(slaa) + s(doy) + s(easting, northing), 
               #method = "REML",
               data = df, 
               family = gaussian(link = "identity"))

gc()
# gam_2018 = gam(chl ~ s(slaa) + s(fsle) + s(doy) + s(easting, northing), 
#                method = "REML",
#                data = df, 
#                family = gaussian(link = "identity"))
# 
# gam_2018 = gam(chl ~ s(doy, by = regions, k=4) + s(doy), #+ s(easting, northing), 
#                # method = "REML",
#                data = df, 
#                family = gaussian(link = "identity"))

# ------------------------------------------------------------------------------
# Initial plot
b = plot(gam_chla, all = TRUE)
s = plot(gam_ssta, all = TRUE)

# ------------------------------------------------------------------------------
# Visualizing the gam
title = c("a)", "b)", "Time")

x1 = as.data.frame(cbind(b[[1]]$x, b[[1]]$fit, b[[1]]$se))
y1 = as.data.frame(cbind(b[[2]]$x, b[[2]]$fit, b[[2]]$se))
z1 = as.data.frame(cbind(b[[3]]$x, b[[3]]$fit, b[[3]]$se))

# x1 = subset(x1, V1 > -0.2)
# x1 = subset(x1, V1 < 0.2)

x2 = as.data.frame(cbind(s[[1]]$x, s[[1]]$fit, s[[1]]$se))
y2 = as.data.frame(cbind(s[[2]]$x, s[[2]]$fit, s[[2]]$se))
z2 = as.data.frame(cbind(s[[3]]$x, s[[3]]$fit, s[[3]]$se))
# ------------------------------------------------------------------------------
# plotting GAM results
shadedline(x    = x1$V1, 
           ylim = c(-0.017, 0.017),
           xlim = c(-0.12, 0.11),
           y0 = x1$V2,
           y1 = x1$V2 + x1$V3*10, 
           y2 = x1$V2 - x1$V3*10,
           title = "",
           line = -1.5,
           adj  = 0.02,
           main = paste(title[1]),
           cex.main = 1.25,
           ylab = "SSTA Effect [C]",
           col  = "grey69",
           xlab = "SLA",
           xaxis = FALSE,
           #labels = format(as.Date(x$V1), "%Y-%m-%d"),
           mar = c(5, 5, 1, 1)
)
axis(side = 1, 
     las = 1, 
     lwd = 2, 
     mgp = c(2, 1, 0), 
     cex.axis = 1.15)


shadedline(x    = x2$V1, 
           ylim = c(-0.28, 0.12),
           xlim = c(-0.12, 0.11),
           y0 = x2$V2,
           y1 = x2$V2 + x2$V3*10, 
           y2 = x2$V2 - x2$V3*10,
           title = "",
           line = -1.5,
           adj  = 0.02,
           main = paste(title[2]),
           cex.main = 1.25,
           ylab = "SSTA Effect [C]",
           col  = "grey69",
           xlab = "SLA",
           xaxis = FALSE,
           #labels = format(as.Date(x$V1), "%Y-%m-%d"),
           mar = c(5, 5, 1, 1)
)
axis(side = 1, 
     las = 1, 
     lwd = 2, 
     mgp = c(2, 1, 0), 
     cex.axis = 1.15)
# ------------------------------------------------------------------------------
# plotting
layout.matrix = matrix(c(1, 
                         2), 
                       nrow = 2,
                       ncol = 2,
                       byrow = TRUE)
nf = graphics::layout(layout.matrix, 
                      heights = c(2, 2))


shadedline(x    = x1$V1, 
           ylim = c(-0.015, 0.015),
           xlim = c(-0.1, 0.1),
           y0 = x1$V2,
           y1 = x1$V2 + x1$V3*10, 
           y2 = x1$V2 - x1$V3*10,
           title = "",
           line = -1.5,
           adj  = 0.02,
           main = paste(title[1]),
           cex.main = 1.25,
           ylab = "SSTA Effect [C]",
           col  = "grey69",
           xlab = "SLA",
           xaxis = FALSE,
           #labels = format(as.Date(x$V1), "%Y-%m-%d"),
           mar = c(5, 5, 1, 1)
)
axis(side = 1, 
     las = 1, 
     lwd = 2, 
     mgp = c(2, 1, 0), 
     cex.axis = 1.15)

shadedline(x = y1$V1, 
           #ylim = c(-0.0025, 0.0015),
           #xlim = c(0, 0.15),
           xaxis = FALSE,
           y0 = y1$V2,
           y1 = y1$V2 + y1$V3*10, 
           y2 = y1$V2 - y1$V3*10,
           title = "",
           line = -1.5,
           adj = 0.02,
           main = paste(title[2]),
           cex.main = 1.25,
           ylab = "",
           col = "grey69",
           xlab = "FSLE",
           #labels = format(as.Date(y$V1), "%Y-%m-%d"),
           mar = c(5, 3, 1, 3)
)
axis(side = 1, 
     las = 1, 
     lwd = 2, 
     mgp = c(2, 1, 0), 
     cex.axis = 1.15)

shadedline(x = z1$V1, 
           xaxis = FALSE,
           y0 = z1$V2,
           y1 = z1$V2 + z1$V3 * 10, 
           y2 = z1$V2 - z1$V3 * 10,
           title = "",
           line = -1.5,
           adj = 0.02,
           main = paste(title[3]),
           cex.main = 1.25,
           ylab = "CHL resp",
           col = "grey69",
           xlab = "DOY",
           mar = c(5, 3, 1, 3)
)
axis(side = 1, 
     las = 1, 
     lwd = 2, 
     mgp = c(2, 1, 0), 
     cex.axis = 1.15)
# ------------------------------------------------------------------------------

# ------------------------------------------------------------------------------

# ------------------------------------------------------------------------------

# ------------------------------------------------------------------------------
