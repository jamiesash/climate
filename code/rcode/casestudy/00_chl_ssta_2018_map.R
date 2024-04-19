
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

chla = anomalize(chl, detrend = FALSE)
gc()
# ------------------------------------------------------------------------------

# probs not needed
# by time
sdate = as.Date("2018-07-05")
edate = as.Date("2018-10-15")
chl  = timesnip(chl,  sdate, edate)
chla = timesnip(chla,  sdate, edate)
ssta  = timesnip(ssta,  sdate, edate)
gc()


# ------------------------------------------------------------------------------
# removing coastal effects
ssta = bufcoast(ssta)
gc()
chla = bufcoast(chla)
chl = bufcoast(chl)
gc()

# ------------------------------------------------------------------------------
# function for plotting
map = function(ras, colmap, zlim,
               main = "Sea surface temperature anomaly"){
  
  ras   <- raster::clamp(ras, zlim[1], zlim[2])
  wdmap <- getMap(resolution = "high")
  e = extent(ras)
  
  image(ras,
        xlim = e[1:2],
        ylim = e[3:4],
        zlim = zlim,
        xlab = "",
        ylab = "",
        col = colmap,
        axes = FALSE)
  
  plot(wdmap,
       xlim = e[1:2],
       ylim = e[3:4],
       asp = 1,
       bg = "black",
       border = "black",
       col = "black",
       #wrap = c(-180, 180),
       add = TRUE)
  
  axis(side = 2, 
       las = 2, 
       lwd = 2, 
       at = c(18, 22, 26, 30, 34, 38),
       mgp = c(1, 0.75, 0), 
       cex.axis = 1.25)
  
  axis(side = 1, 
       las = 1, 
       lwd = 2,
       mgp = c(2, 1, 0),    
       at = c(-170, -160, -150, -140, -130),
       cex.axis = 1.25)
  
  #title(main = list(main), line = line,  adj = adj, cex.main = cex.main)
  title(main = main, cex.main = 1.75, line = 0.25, adj = 0)
  
  title(ylab = "Latitude", cex.lab = 1.5, line = 2.5)
  
  title(xlab = "Longitude", cex.lab = 1.5, line = 2.5)
  
  box(which = "plot",
      lty = "solid",
      lwd = 3,
      col = colvect("grey12", alpha = 0.9))
}

# ------------------------------------------------------------------------------
# idk
fronts = ssta[[75]]

# making negative values positive
idx = fronts < 0
fronts[idx] = fronts[idx] * -1

# making 0 the greatest value

test = scale(values(fronts), from = 1, to = 0)

vals = scale(values(fronts), from = 1, to = 0)
values(fronts) = vals

plot(fronts)

# ------------------------------------------------------------------------------
# well I guess I also map chl and ssta here thats not good codeing
drawPalette(zlim = c(-0.05, 0.15), 
            col = cmocean("rain")(20),
            plot = TRUE,
            pos  = 4,
            zlab = "Chlorophyll Anomaly")
par(mar = c(5, 4, 2, 5))
map(chla[[75]],
    zlim = c(-0.05, 0.15),
    colmap = cmocean("rain")(50),
    main = paste("Chlorophyll Anomaly ", getZ(chla[[75]]), sep = ""))

drawPalette(zlim = c(-1.75, 3.57), 
            col = cmocean("balance")(20),
            plot = TRUE,
            pos  = 4,
            zlab = "Temperature Anomaly")
par(mar = c(5, 4, 2, 5))
map(ssta[[75]], 
    colmap = cmocean("balance")(50),
    zlim = c(-1.75, 3.75),
    main = paste("Sea Surface Temperature Anomaly ", getZ(ssta[[75]]), sep = ""))
