# ------------------------------------------------------------------------------
### Liraries and functions

library(anytime)
library(oce)
library(cmocean)
library(rworldmap)
library(rworldxtra)
library(ncdf4)
library(raster)
library(lubridate)

source("functions.R")

# ------------------------------------------------------------------------------
### Loading the dataset

# using a two degree box around aloha
lons = c(-175, -125)
lats = c(17, 45)

# load_nc
chl = load_nc(path = "/home/jamesash/blooms/data/",
              patt = "chl_1998_2023_l4_month_multi_4k.nc",
              vars = c("latitude", "longitude", "CHL", "time"))
 
chl = oreant(chl, flip = 'y', t1 = TRUE)

# print("Step 1. Loaded and flipped.")

# ------------------------------------------------------------------------------

# Remmove the seasonal climitological trend. 
chla = anomalize(chl, detrend = FALSE)
rm(chl)

# print("Step 2. Anomalized.")

# Remove all but the summer months
chla = subsum(chla, mnth = 6:10)
# print("Step 3. Summer left standing.")

# ------------------------------------------------------------------------------
### Math on dataset

cmap = calc(chla, fun = mean, na.rm = TRUE)
zlim = c(0, 0.012)
cmap = raster::clamp(cmap, zlim[1], zlim[2])

# print("Step 4. Averaged.")

# ------------------------------------------------------------------------------

wdmap <- getMap(resolution = "high")
dt = gsub("-", "", as.character(Sys.Date()))
colmap = cmocean("rain")(50)
e = extent(cmap)

# print("Step 5. Worldmap loaded.")

# ------------------------------------------------------------------------------
### Plotting the data 
# Plotting CHL anomaly official plot

pdf(paste("chla_map", dt, ".pdf", sep = ""),   # The directory you want to save the file in
    width = 6, # The width of the plot in inches
    height = 4,
    pointsize = 10) # The height of the plot in inches

#par(mar = c(4,4,2,1))
image(cmap,
      ylim = c(18, 40),
      xlim = c(-175, -130),
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
     lwd = 1)
axis(side = 2, 
     las = 2, 
     lwd = 1, 
     at = c(18, 22, 26, 30, 34),
     mgp = c(1, 0.75, 0), 
     cex.axis = 1)
axis(side = 1, 
     las = 1, 
     lwd = 1,
     mgp = c(2, 1, 0),    
     at = c(-170, -160, -150, -140, -130),
     cex.axis = 1)
title(xlab = "Longitude", line = 2.6, cex.lab = 1)
title(ylab = "Latitude", line = 2.6, cex.lab = 1)
box(which = "plot",
    lty = "solid",
    lwd = 1.5,
    col = "black")

dev.off()


