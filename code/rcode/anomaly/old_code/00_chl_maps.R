# Loading packages and all prewriten functions
# setwd("C:\\Users\\james\\Desktop\\jamieslife\\analysis")
source("code\\libraries.R")
source("code\\functions.R")
rasterOptions(maxmemory = 100e+10)

# ------------------------------------------------------------------------------
# using a two degree box around aloha
lons = c(-175, -125)
lats = c(17, 45)
e = extent(lons, lats)

# ------------------------------------------------------------------------------
# need to fix my opendap function
# Set variables
lat_varid = "lat"
lon_varid = "lon"
sdate = as.Date("2002-01-01")
edate = as.Date("2022-11-01")
var = "CHL"
url = "https://jash:5.Pellegrino@my.cmems-du.eu/thredds/dodsC/cmems_obs-oc_glo_bgc-plankton_my_l4-multi-4km_P1M?"
origin = "1900-01-01"

data = nc_open(url, verbose = FALSE, write = FALSE)
lat  = ncvar_get(data, varid = lat_varid)
lon  = ncvar_get(data, varid = lon_varid)
time = ncvar_get(data, varid = "time")
time = as.Date(time, origin = origin)

idx_lat <- which(lat > lats[1] & lat < lats[2])
idx_lon <- which(lon > lons[1] & lon < lons[2])
idx_time <- which(time >= sdate & time <= edate)


idx_ras <- paste("CHL",
                 paste("[", range(idx_time)[1], ":1:", range(idx_time)[2], "]", sep = ""),
                 paste("[", range(idx_lat)[1],  ":1:", range(idx_lat)[2],  "]", sep = ""),
                 paste("[", range(idx_lon)[1],  ":1:", range(idx_lon)[2],  "]", sep = ""),
                 sep = "")

idx_time <- paste("time", paste("[", range(idx_time)[1], ":1:",range(idx_time)[2], "]", sep = ""), sep = "")
idx_lat  <- paste(lat_varid, paste("[", range(idx_lat)[1], ":1:",range(idx_lat)[2],  "]", sep = ""), sep = "")
idx_lon  <- paste(lon_varid, paste("[", range(idx_lon)[1], ":1:",range(idx_lon)[2],  "]", sep = ""), sep = "")
idx <- paste(idx_lat, idx_lon, idx_time, idx_ras, sep = ",")

url <- paste(url, idx, sep = "")

nc_close(data)
rm(data)

data = nc_open(url, verbose = FALSE, write = FALSE)

lat  = ncvar_get(data, varid = lat_varid)
lon  = ncvar_get(data, varid = lon_varid)
time = ncvar_get(data, varid = "time")
time = as.Date(time, origin = origin)
ras  = ncvar_get(data)
nc_close(data)

s = dim(ras)
ras = raster::brick(ras)
ras = t(ras)
ras = setZ(ras, z = as.Date(time, origin = org), name = "time")
extent(ras) = extent(min(lons), max(lons), min(lats), max(lats))

chl = ras
rm(ras)
gc()

# ------------------------------------------------------------------------------
# Calculate CHL anomaly
# I receive an warning when detrending for cells of all NA values (hawaii)
# chl anom. calc. needs at least one full year to work correctly
chla = anomalize(chl, detrend = FALSE)
rm(chl)
gc()

# ------------------------------------------------------------------------------
# subset the summer months
chl  = subsum(chl)
chla = subsum(chla)
gc()

# ------------------------------------------------------------------------------
# Subset summer and average fro climatology figure of chl anomaly
camap   = calc(chla, fun = mean, na.rm = TRUE)
zlim = c(0, 0.012)
camap = raster::clamp(camap, zlim[1], zlim[2])

# Subset summer and average fro climatology figure
cmap   = calc(chl, fun = mean, na.rm = TRUE)

zlim = c(0.05, 0.10)
cmap = raster::clamp(cmap, zlim[1], zlim[2])

# ------------------------------------------------------------------------------
# 2018 bloom plots
sdate = as.Date("2018-06-01")
edate = as.Date("2018-11-15")
chl_2018 = timesnip(chl,  sdate, edate)
chla_2018 = timesnip(chla,  sdate, edate)
chla_2018 = bufcoast(chla_2018)

# ------------------------------------------------------------------------------
# For CHL anomaly as well?
zlim = c(0.02, 0.15)
cmap_2018   = calc(chl_2018, fun = mean, na.rm = TRUE)
cmap_2018 = raster::clamp(cmap_2018, lower = zlim[1], upper=zlim[2])

zlim = c(-0.005, 0.025)
cmap_2018_anom   = calc(chla_2018, fun = mean, na.rm = TRUE)
cmap_2018_anom = raster::clamp(cmap_2018_anom, lower=zlim[1], upper=zlim[2])
# ------------------------------------------------------------------------------

# Plotting CHL anomaly official plot
wdmap = getMap(resolution = "high")
e = extent(camap)
colmap = cmocean("rain")(50)

png(filename= paste("chla_map", Sys.Date(), ".png", sep = ""),
    width = 3.5,
    height = 6,
    units = "in",
    res = 300,
    pointsize = 10)

lay = matrix(c(1, 2), 
             nrow = 1, 
             byrow = TRUE)

nf = graphics::layout(lay, 
                      widths = c(14, 1))

par(mar = c(4,4,2,1))
image(camap,
      ylim = c(18, 35),
      xlim = c(-170, -130),
      xlab = "",
      ylab = "",
      zlim = zlim,
      col = colmap,
      axes = FALSE)
plot(wdmap,
     xlim = e[3:4],
     ylim = e[1:2],
     asp = 1,
     bg = "black",
     border = "black",
     col = "black",
     add = TRUE,
     lwd = 4)
axis(side = 2, 
     las = 2, 
     lwd = 2, 
     at = c(18, 22, 26, 30, 34),
     mgp = c(1, 0.75, 0), 
     cex.axis = 1.25)
axis(side = 1, 
     las = 1, 
     lwd = 2,
     mgp = c(2, 1, 0),    
     at = c(-170, -160, -150, -140, -130),
     cex.axis = 1.25)
#title(main = list(main), line = line,  adj = adj, cex.main = cex.main)
title(xlab = "Longitude", line = 2.6, cex.lab = 1.5)
title(ylab = "Latitude", line = 2.6, cex.lab = 1.5)
box(which = "plot",
    lty = "solid",
    lwd = 4,
    col = colvect("grey12", alpha = 0.9))

drawPalette(zlim = zlim,
            col  = colmap,
            plot = TRUE,
            pos  = 4,
            zlab = "",
            cex  = 1.25,
            fullpage = TRUE)

dev.off()

# ------------------------------------------------------------------------------
# Plotting CHL raw official plot
zlim = c(0.05, 0.10)

wdmap = getMap(resolution = "high")
e = extent(cmap)
colmap = cmocean("rain")(50)

png(filename= paste("chl_map", Sys.Date(), ".png", sep = ""),
    width = 3.5,
    height = 6,
    units = "in",
    res = 300,
    pointsize = 10)

lay = matrix(c(1, 2), 
             nrow = 1, 
             byrow = TRUE)

nf = graphics::layout(lay, 
                      widths = c(14, 1))

par(mar = c(4,4,2,1))
image(cmap,
      ylim = c(18, 35),
      xlim = c(-170, -130),
      xlab = "",
      ylab = "",
      zlim = zlim,
      col = colmap,
      axes = FALSE)
plot(wdmap,
     xlim = c(-170, -130),
     ylim = c(17, 35),
     asp = 1,
     bg = "black",
     border = "black",
     col = "black",
     add = TRUE,
     lwd = 4)
axis(side = 2, 
     las = 2, 
     lwd = 2, 
     at = c(18, 22, 26, 30, 34),
     mgp = c(1, 0.75, 0), 
     cex.axis = 1.25)
axis(side = 1, 
     las = 1, 
     lwd = 2,
     mgp = c(2, 1, 0),    
     at = c(-170, -160, -150, -140, -130),
     cex.axis = 1.25)
#title(main = list(main), line = line,  adj = adj, cex.main = cex.main)
title(xlab = "Longitude", line = 2.6, cex.lab = 1.5)
title(ylab = "Latitude", line = 2.6, cex.lab = 1.5)
box(which = "plot",
    lty = "solid",
    lwd = 4,
    col = colvect("grey12", alpha = 0.9))

drawPalette(zlim = c(0.03, 0.12),
            col  = colmap,
            plot = TRUE,
            pos  = 4,
            zlab = "",
            cex  = 1.25,
            fullpage = TRUE)
dev.off()

# ------------------------------------------------------------------------------
# Plotting CHL 2018 official plot
wdmap = getMap(resolution = "high")
e = extent(cmap_2018)
zlim = c(0.03, 0.15)
colmap = cmocean("rain")(50)

png(filename= paste("chl2018_map", Sys.Date(), ".png", sep = ""),
    width = 3.5,
    height = 6,
    units = "in",
    res = 300,
    pointsize = 10)

lay = matrix(c(1, 2), 
             nrow = 1, 
             byrow = TRUE)

nf = graphics::layout(lay, 
                      widths = c(14, 1))

par(mar = c(4,4,2,1))
image(cmap_2018,
      ylim = c(18, 35),
      xlim = c(-170, -130),
      xlab = "",
      ylab = "",
      zlim = zlim,
      col = colmap,
      axes = FALSE)
plot(wdmap,
     xlim = c(-170, -130),
     ylim = c(17, 35),
     asp = 1,
     bg = "black",
     border = "black",
     col = "black",
     add = TRUE,
     lwd = 4)
axis(side = 2, 
     las = 2, 
     lwd = 2, 
     at = c(18, 22, 26, 30, 34),
     mgp = c(1, 0.75, 0), 
     cex.axis = 1.25)
axis(side = 1, 
     las = 1, 
     lwd = 2,
     mgp = c(2, 1, 0),    
     at = c(-170, -160, -150, -140, -130),
     cex.axis = 1.25)
#title(main = list(main), line = line,  adj = adj, cex.main = cex.main)
title(xlab = "Longitude", line = 2.6, cex.lab = 1.5)
title(ylab = "Latitude", line = 2.6, cex.lab = 1.5)
box(which = "plot",
    lty = "solid",
    lwd = 4,
    col = colvect("grey12", alpha = 0.9))

drawPalette(zlim = zlim,
            col  = colmap,
            plot = TRUE,
            pos  = 4,
            zlab = "",
            cex  = 1.25,
            fullpage = TRUE)
dev.off()
