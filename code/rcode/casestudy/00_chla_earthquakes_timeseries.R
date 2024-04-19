
source("code/functions.R")
source("code/libraries.R")
library(imputeTS)
rasterOptions(maxmemory = 120e+10, memfrac = 0.9)

quak <- read.csv("data/earthquakes/quakes_1997_20230601.csv")
quak$time <- substr(quak$time, start = 1, stop = 10)
quak$time <- as.Date(quak$time)

gry <- colvect("grey22", alpha = 0.3)

# ------------------------------------------------------------------------------
# using a two degree box around aloha
lons = c(-158, -144)
lats = c(22, 26)
e = extent(lons, lats)

# same as calibration
# lons = c(-159, -157)
# #lons = c(-158.05, -157.95)
# lats = c( 22,   23)
# e = extent(lons, lats)
# loading chl the grossway -----------------------------------------------------

# I have a function in functions.R that does this but do not want to take the 
# time to update the code

# Set variables
lat_varid = "lat"
lon_varid = "lon"
sdate = as.Date("1997-01-01")
edate = as.Date("2023-05-01")
var = "CHL"
#url = "https://jash:5.Pellegrino@my.cmems-du.eu/thredds/dodsC/cmems_obs-oc_glo_bgc-plankton_my_l4-multi-4km_P1M?"
# daily and same as calibration
url = "https://jash:5.Pellegrino@my.cmems-du.eu/thredds/dodsC/cmems_obs-oc_glo_bgc-plankton_my_l3-multi-4km_P1D?"
origin = "1900-01-01"

data = nc_open(url, verbose = FALSE, write = FALSE)
lat  = ncvar_get(data, varid = lat_varid)
lon  = ncvar_get(data, varid = lon_varid)
time = ncvar_get(data, varid = "time")
time = as.Date(time, origin = origin)

idx_lat = which(lat > lats[1] & lat < lats[2])
idx_lon = which(lon > lons[1] & lon < lons[2])
idx_time = which(time >= sdate & time <= edate)


idx_ras = paste("CHL",
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

# lons = c(-158.1, -157.95)
# lats = c( 22.65,   22.8)
# e = extent(lons, lats)
# chl = raster::crop(chl, e)
# dim(chl)

# ------------------------------------------------------------------------------
# remove coastal effects
chl = bufcoast(chl, 
               region = "Hawaiian Islands", 
               path = "data/USMaritimeLimitsAndBoundariesSHP")

# ------------------------------------------------------------------------------
# calculate the chl anomaly 
chla = anomalize(chl, detrend = TRUE)
# chla = anomalize(chl, detrend = FALSE)
gc()

# ------------------------------------------------------------------------------
# stats to be plotted
tc   <- cellStats(chl, stat="mean", na.rm = TRUE)
tc   <- as.vector(tc)
# tc[is.nan(tc)] <- NA
time <- getZ(chl)

tca   <- cellStats(chla, stat="mean", na.rm = TRUE)
tca   <- as.vector(tca)

u = mad(tca, na.rm = TRUE)
o = median(tca, na.rm = TRUE)

# ------------------------------------------------------------------------------
# plotting just chl timeseries

dt = gsub("-", "", as.character(Sys.Date()))
pdf(paste("figures\\blooms2_", dt, ".pdf", sep = ""),   # The directory you want to save the file in
    width = 6, # The width of the plot in inches
    height = 4,
    pointsize = 10) # The height of the plot in inches

par(mar = c(4, 4, 4, 1))
plot(time, tca, 
     pch = 4, 
     cex = 1, 
     xlab = "",
     ylab = "",
     col = gry, 
     axes = FALSE,
     xlim = c(as.Date("1997-10-01"), as.Date("2023-01-01"))
     )
abline(v = as.Date("2018-04-25"), col = "firebrick", lwd = 1, lty = 2)
# abline(h = u, 
#        col="red", lwd = 1, lty = 1)
axis(side = 2,
     las = 2, 
     lwd = 1, 
     mgp = c(1, 0.75, 0), 
     cex.axis = 1)
axis(side = 1,
     labels = as.character(year(seq(min(time), max(time), length = 10))),
     at =  as.numeric(seq(min(time), max(time), length = 10)), 
     las = 1, 
     lwd = 1, 
     mgp = c(2, 1, 0), 
     cex.axis = 1)
#title(main = list("satellite chlorophyll anomaly at Stn. ALOHA", cex = 1.5),
#      line = -1.5,
#      adj  = 0.05)
title(ylab    = expression(Mean ~ CHL ~ (mg ~ m^{-3})), 
      cex.lab = 1,
      line    = 2.5)
title(xlab    = "Year", 
      cex.lab = 1,
      line    = 2.5)
box(which = "plot", lty = "solid", lwd = 1.5, col = "black")
dev.off()
# ------------------------------------------------------------------------------
# plotting just earthquaes

dt = gsub("-", "", as.character(Sys.Date()))
pdf(paste("..\\figures\\quakes_", dt, ".pdf", sep = ""),   # The directory you want to save the file in
    width = 6, # The width of the plot in inches
    height = 4,
    pointsize = 10) # The height of the plot in inches

par(mar = c(4, 4, 4, 1))
plot(quak$time, quak$mag, 
     pch = 4, 
     cex = 1, 
     xlab = "",
     ylab = "",
     col = gry, 
     axes = FALSE,
     xlim = c(as.Date("1997-10-01"), as.Date("2023-01-01"))
     )
abline(v = as.Date("2018-04-25"), col = "firebrick", lwd = 1, lty = 2)
axis(side = 1,
     labels = as.character(year(seq(min(time), max(time), length = 10))),
     at =  as.numeric(seq(min(time), max(time), length = 10)), 
     las = 1, 
     lwd = 1, 
     mgp = c(2, 1, 0), 
     cex.axis = 1)
axis(side = 2,
     las = 2, 
     lwd = 1, 
     mgp = c(1, 0.75, 0), 
     cex.axis = 1)
# title(main = list("Earthquakes along the Hawaiian Islands", cex = 1),
#       line = -1.5,
#       adj  = 0.05)
title(ylab    = "Magnitude", 
      cex.lab = 1,
      line    = 2.5)
title(xlab    = "Year", 
      cex.lab = 1,
      line    = 2.5)
box(which = "plot", lty = "solid", lwd = 1.5, col = "black")
dev.off()

# ------------------------------------------------------------------------------
# both in the same figure

dt = gsub("-", "", as.character(Sys.Date()))
pdf(paste("figures\\blooms_quakes_", dt, ".pdf", sep = ""),   # The directory you want to save the file in
    width = 8, # The width of the plot in inches
    height = 8,
    pointsize = 10) # The height of the plot in inches

layout( matrix(c(1,
                 2),
               ncol = 1,
               nrow = 2, 
               byrow = TRUE),
        heights = c(1, 1))

par(mar = c(0.25, 4, 4, 1))
plot(time, tca, 
     pch = 4, 
     cex = 1, 
     xlab = "",
     ylab = "",
     col = gry, 
     axes = FALSE,
     xlim = c(as.Date("1997-01-01"), as.Date("2023-05-01"))
)
abline(v = as.Date("2018-04-25"), col = "firebrick", lwd = 1, lty = 2)
axis(side = 2,
     las = 2, 
     lwd = 1, 
     mgp = c(1, 0.75, 0), 
     cex.axis = 1)
axis(side = 1,
     #labels =  as.character(year(seq(min(time), max(time), length = 10))),
     labels =  rep("", length = 10),
     at =  as.numeric(seq(min(time), max(time), length = 10)), 
     las = 1, 
     lwd = 1, 
     mgp = c(2, 1, 0), 
     cex.axis = 1)
#title(main = list("satellite chlorophyll anomaly at Stn. ALOHA", cex = 1.5),
#      line = -1.5,
#      adj  = 0.05)
title(ylab    = expression(Mean ~ CHL ~ Anomaly ~ (mg ~ m^{-3})), 
      cex.lab = 1,
      line    = 2.5)
box(which = "plot", lty = "solid", lwd = 1.5, col = "black")

par(mar = c(4, 4, 0.25, 1))
plot(quak$time, quak$mag, 
     pch = 4, 
     cex = 1, 
     xlab = "",
     ylab = "",
     col = gry, 
     axes = FALSE,
     xlim = c(as.Date("1997-01-01"), as.Date("2023-05-01"))
)
abline(v = as.Date("2018-04-25"), col = "firebrick", lwd = 1, lty = 2)
axis(side = 1,
     labels = as.character(format(seq(min(time), max(time), length = 10), "%Y-%m")),
     at =  as.numeric(seq(min(time), max(time), length = 10)), 
     las = 1, 
     lwd = 1, 
     mgp = c(2, 1, 0), 
     cex.axis = 1)
axis(side = 2,
     las = 2, 
     lwd = 1, 
     mgp = c(1, 0.75, 0), 
     cex.axis = 1)
# title(main = list("Earthquakes along the Hawaiian Islands", cex = 1),
#       line = -1.5,
#       adj  = 0.05)
title(ylab    = "Earthquake Magnitude", 
      cex.lab = 1,
      line    = 2.5)
title(xlab    = "Year Month", 
      cex.lab = 1,
      line    = 2.5)
box(which = "plot", lty = "solid", lwd = 1.5, col = "black")
dev.off()

# ------------------------------------------------------------------------------
# I think chl anomaly again?

# png(filename= paste("chl_anom", Sys.Date(), ".png", sep = ""),
#     width = 3.5,
#     height = 6,
#     units = "in",
#     res = 300,
#     pointsize = 10)
# 
# par(mar = c(5, 5, 5, 1))
# plot(time, tca, 
#      pch = 19, 
#      cex = 0.75, 
#      xlab = "",
#      ylab = "",
#      col = "grey22", 
#      axes = FALSE,
#      xlim = c(as.Date("1997-10-01"), as.Date("2023-01-01"))
# )
# abline(h = median(tca) + mad(tca), 
#        col="red", lwd = 1, lty = 1)
# axis(side = 2,
#      las = 2, 
#      lwd = 2, 
#      mgp = c(1, 0.75, 0), 
#      cex.axis = 1.5)
# axis(side = 1,
#      labels = as.character(year(seq(min(time), max(time), length = 10))),
#      at =  as.numeric(seq(min(time), max(time), length = 10)), 
#      las = 1, 
#      lwd = 2, 
#      mgp = c(2, 1, 0), 
#      cex.axis = 1.5)
# title(main = list("Satellite chlorophyll anomaly at Stn. ALOHA", cex = 1.5),
#       line = -1.5,
#       adj  = 0.05)
# title(ylab    = "CHL [mg/m^3]", 
#       cex.lab = 1.5,
#       line    = 3.5)
# title(xlab    = "Year", 
#       cex.lab = 1.5,
#       line    = 2.5)
# box(which = "plot", lty = "solid", lwd = 3, col = "grey12")
# dev.off()
# 
# # ------------------------------------------------------------------------------
# # plotting raw chl timeseries
# png(filename= paste("earthquakes", Sys.Date(), ".png", sep = ""),
#     width = 3.5,
#     height = 6,
#     units = "in",
#     res = 300,
#     pointsize = 10)
# 
# par(mar = c(5, 5, 5, 1))
# plot(time, tc, 
#      pch = 19, 
#      cex = 0.75, 
#      xlab = "",
#      ylab = "",
#      col = "grey22", 
#      axes = FALSE,
#      xlim = c(as.Date("1997-10-01"), as.Date("2023-01-01"))
# )
# abline(h = median(tc) + mad(tc), 
#        col="red", lwd = 1, lty = 1)
# axis(side = 2,
#      las = 2, 
#      lwd = 2, 
#      mgp = c(1, 0.75, 0), 
#      cex.axis = 1.5)
# axis(side = 1,
#      labels = as.character(year(seq(min(time), max(time), length = 10))),
#      at =  as.numeric(seq(min(time), max(time), length = 10)), 
#      las = 1, 
#      lwd = 2, 
#      mgp = c(2, 1, 0), 
#      cex.axis = 1.5)
# title(main = list("Satellite chlorophyll at Stn. ALOHA", cex = 1.5),
#       line = -1.5,
#       adj  = 0.05)
# title(ylab    = "CHL [mg/m^3]", 
#       cex.lab = 1.5,
#       line    = 3.5)
# title(xlab    = "Year", 
#       cex.lab = 1.5,
#       line    = 2.5)
# box(which = "plot", lty = "solid", lwd = 3, col = "grey12")
# dev.off()








