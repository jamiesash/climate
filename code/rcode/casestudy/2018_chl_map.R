library(mgcViz)
library(lubridate)
library(terra)
library(cmocean)
library(fields) 

# ------------------------------------------------------------------------------
### Script To Do

# 1. Make those subplots. 
# 2. Save a figure of peak bloom.

# ------------------------------------------------------------------------------
# Load SLA, CHL and SSTA data sets. 

setwd("/home/jamie/projects/climate/code/rcode/casestudy")
chl = rast("../../../data/chl/chl_2018_glob_daily_multi_l4_4k.nc")

### Loading the dataset
states = vect("../../../data/states/s_05mr24.shp")
idx = states[[c(1)]] == "HI"
hawaii = states[idx]

# ------------------------------------------------------------------------------
colmap = cmocean("rain")(100)
zmap = cmocean("rain")(30)

# Truncating lower bounds.
bounds = c(0.02, 0.3)

chl_clamp = clamp(chl, lower=bounds[1], upper=bounds[2], values=TRUE)

par(mar=c(4, 4, 3, 7))
image(chl_clamp[[100]],
      col = colmap,
      ylim = c(16, 39),
      xlim = c(-175, -130),
      axes = FALSE)
plot(hawaii,
     col = "black",
     add = TRUE)
image.plot(zlim = c(0.02, 0.15),
           smallplot= c(0.9, 0.95, 0.11, 0.92),
           axis.args = c(lwd = 0, lwd.ticks = 0.5),
           col = zmap,
           nlevel = 10,
           add = TRUE, 
           legend.only = TRUE)
axis(side = 2, 
     las = 2, 
     lwd = 0,
     lwd.ticks = 1,
     at = c(18, 22, 26, 30, 34, 38),
     mgp = c(1, 0.75, 0), 
     cex.axis = 1)
axis(side = 1, 
     las = 1, 
     lwd = 0,
     lwd.ticks = 1,
     mgp = c(2, 1, 0),    
     at = c(-173, -168, -163, -158, -153, -148, -143, -138, -133),
     cex.axis = 1)
box(lwd = 1, col = "black")
title(expression(GlobColor~CHL~(mg~m^-3~year^-1)))
title(ylab = "Latitude", line = 2)
title(xlab = "Longitude", line = 2)
