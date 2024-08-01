library(mgcViz)
library(lubridate)
library(terra)
library(cmocean)
library(fields) 
library(oce)
library(rasterVis)
# ------------------------------------------------------------------------------
### Script To Do

# 1. Make those subplots. 
# 2. Save a figure of peak bloom.

# ------------------------------------------------------------------------------
# Load SLA, CHL and SSTA data sets. 

setwd("/home/jamie/projects/climate/code/rcode/casestudy")
# chl = rast("../../../data/chl/chl_2018_glob_daily_multi_l4_4k.nc")
chl = rast("../../../data/chl/chl_2018_glob_monthly_multi_l4_4k.nc")

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
image(chl_clamp[[1]],
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

# ------------------------------------------------------------------------------
# Truncating lower bounds.
bounds = c(0.05, 100)
chl_clamp = clamp(chl, lower=bounds[1], upper=bounds[2], values=FALSE)
bounds = c(0, 0.25)
chl_clamp = clamp(chl_clamp, lower=bounds[1], upper=bounds[2], values=TRUE)

colmap = cmocean("algae")(100)

png(filename= paste("../../../figures/chl_2018_panel_map", Sys.Date(), ".png", sep = ""),
    width = 9.7,#
    height = 6,
    units = "in",
    res = 300,
    pointsize = 12)

#pdf(paste("../../../figures/chl_2018_panel_map", Sys.Date(), ".pdf", sep = ""),  
#    width = 9.7, # inches
#    height = 6)
    #pointsize = 7) # inches

lay = matrix(c(c(1, 2, 5), c(3, 4, 5)), 
             nrow = 2,
             ncol = 3,
             byrow = TRUE)
nf = graphics::layout(lay, 
                      widths = c(4, 4, 0.5))

# plot 1
par(mar=c(1, 4, 4, 0))
image(chl_clamp[[2]],
      col = colmap,
      ylim = c(16, 39),
      #asp = 0.83,
      xlim = c(-175, -130),
      axes = FALSE)
plot(hawaii,
     col = "black",
     add = TRUE)
axis(side = 2, 
     las = 2, 
     lwd = 0,
     lwd.ticks = 1,
     at = c(18, 22, 26, 30, 34, 38),
     mgp = c(1, 0.75, 0), 
     cex.axis = 1)
box(lwd = 1.25, col = "black")
title(ylab = "Latitude",  line =2)
title(main = as.character(time(chl_clamp[[2]])), 
      line = -1.2, 
      adj = 0.02,
      cex.main= 1)

# plot 2
par(mar=c(1, 2, 4, 2))
image(chl_clamp[[3]],
      col = colmap,
      #asp = 0.83,
      ylim = c(16, 39),
      xlim = c(-175, -130),
      axes = FALSE)
plot(hawaii,
     col = "black",
     add = TRUE)
box(lwd = 1.25, col = "black")
title(main = as.character(time(chl_clamp[[3]])), 
      line = -1.2, 
      adj = 0.02,
      cex.main= 1)

# plot 3
par(mar=c(4, 4, 1, 0))
image(chl_clamp[[4]],
      col = colmap,
      #asp = 0.83,
      ylim = c(16, 39),
      xlim = c(-175, -130),
      axes = FALSE)
plot(hawaii,
     col = "black",
     add = TRUE)
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
box(lwd = 1.25, col = "black")
title(ylab = "Latitude",  line =2)
title(xlab = "Longitude", line = 2)
title(main = as.character(time(chl_clamp[[4]])), 
      line = -1.2, 
      adj = 0.02,
      cex.main= 1)

# plot 4 
par(mar=c(4, 2, 1, 2))
image(chl_clamp[[5]],
      col = colmap,
      ylim = c(16, 39),
      #asp = 0.83,
      xlim = c(-175, -130),
      axes = FALSE)
plot(hawaii,
     col = "black",
     add = TRUE)
axis(side = 1, 
     las = 1, 
     lwd = 0,
     lwd.ticks = 1,
     mgp = c(2, 1, 0),    
     at = c(-173, -168, -163, -158, -153, -148, -143, -138, -133),
     cex.axis = 1)
box(lwd = 1.25, col = "black")
title(xlab = "Longitude", line = 2)
title(main = as.character(time(chl_clamp[[5]])), 
      line = -1.2, 
      adj = 0.02,
      cex.main= 1)

par(mar=c(4, 0, 4, 15))
drawPalette(zlim = c(0.05, 0.25),
            col  = colmap,
            plot = TRUE,
            pos  = 4,
            zlab = "",
            cex  = 1.25,
            fullpage = TRUE)

mtext(expression(The~2018~CHL~bloom~(mg~m^-3)), side = 3, line = -3, outer = TRUE)
dev.off()

# title(expression(GlobColor~CHL~(mg~m^-3~year^-1)))
# ------------------------------------------------------------------------------

bounds = c(0.05, 100)
chl_clamp = clamp(chl, lower=bounds[1], upper=bounds[2], values=FALSE)
bounds = c(0, 0.25)
chl_clamp = clamp(chl_clamp, lower=bounds[1], upper=bounds[2], values=TRUE)
colmap = cmocean("algae")(100)


pax = list(xat=c(-173, -165, -157, -149, -141, -133),
           yat=c(18, 24, 30, 36),
           cex.axis = 1,
           lab = c())
plg = list(size = c(2, 1.5), 
           tic = "out", 
           cex = 1, 
           at = c(0.05, 0.07, 0.09, 0.11, 0.13, 0.15, 0.17, 0.19))

png(filename= paste("../../../figures/chl2018_panel_map", Sys.Date(), ".png", sep = ""),
    width = 5.8,
    height = 3,
    units = "in",
    res = 300,
    pointsize = 10)

panel(chl_clamp[[c(2, 3, 4, 5)]], 
      fun=\() plot(hawaii, add = TRUE, col = "black"),
      ylim = c(16, 38.5),
      xlim = c(-175, -130),
      pax = pax,
      xlab = "test",
      ylab= "test",
      smooth = TRUE,
      box = TRUE,
      plg = plg,
      col = colmap)

levelplot(chl_clamp[[c(1, 2, 3, 4)]], 
          col.regions = colmap,
          main = "test",
          names.attr = c("", "", "", ""))
#title(main = expression(Excerpts~from~the~2018~CHL~bloom~(mg~m^-3)), line = 3.65, cex.main = 1)
#title(ylab = "Latitude", line = 3)
#title(xlab = "Longitude", line = 4.5)
dev.off()











