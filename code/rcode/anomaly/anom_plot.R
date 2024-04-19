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

# load_nc
# chl = load_nc(path = "/home/jamesash/blooms/data/",
#               patt = "chl_1998_2023_l4_month_multi_4k.nc",
#               vars = c("latitude", "longitude", "CHL", "time"))
 

# print("Step 1. Loaded and flipped.")
# chl = oreant(chl, flip = 'y', t1 = TRUE)

chl = rast("/home/jamesash/blooms/data/chl_2022_2023_l3_daily_multi_4k.nc")

# ------------------------------------------------------------------------------

# Perform the regression right awway fug it. 
t = 1:nlyr(chl)
cli = regress(chl, t, na.rm = TRUE)
# I need to remove the slope somehow. 
# output is a spacial raster.  

# -------------------------------------------------------------------------------
### Or do it ussing matrix multiplication
# This is considered a local regression since it is by grid cell. 
# taken from https://gis.stackexchange.com/questions/403811/linear-regression-analysis-in-r
x = cbind(1, 1:nlyr(chl))
## pre-computing constant part of least squares
invxtx = solve(t(x) %*% x) %*% t(x)
## [2] is to just get the slope
quickfun = function(y) (invxtx %*% y)[2]
slope = app(chl, quickfuni) 

# -------------------------------------------------------------------------------
# using the lm function doirectly for na.rm
fun = function(x) { 
	t = 1:length(x)
	lm(x ~ t, na.exclude = TRUE)$coefficients[2] 
	}

slp = terra::app(chl, fun)

# plot and save teh plot
dt = gsub("-", "", as.character(Sys.Date()))
pdf(paste("cli_map_", dt, ".pdf", sep = ""),   # The directory you want to save the file in
    width = 5, # The width of the plot in inches
    height = 4,
    pointsize = 10) # The height of the plot in inches

plot(slope)

dev.off()



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
    width = 5, # The width of the plot in inches
    height = 4,
    pointsize = 10) # The height of the plot in inches

#par(mar = c(4,4,2,1))
image(cmap,
      ylim = c(17, 36),
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

