# ------------------------------------------------------------------------------
### Liraries and functions
library(cmocean)
library(anytime)
library(terra)
# for image.plot colorbar
library(fields) 

print("2. Librares loaded.")

# ------------------------------------------------------------------------------
### Loading the dataset
states = vect("../../data/states/s_05mr24.shp")
idx = states[[c(1)]] == "HI"
hawaii = states[idx]

cli = rast("../../data/cli_mon_sum_20240118.nc")
cli = cli[[2]]
# per year not month
cli = cli*12

print("3. Data loaded")

# ------------------------------------------------------------------------------
colmap = cmocean("delta")(100)
zmap = cmocean("delta")(30)

# Truncating lower bounds.
bounds = c(-0.02, 0)
lcli = clamp(cli, lower=bounds[1], upper=bounds[2], values=TRUE)

# Truncating upper bounds
bounds = c(0, 0.005)
ucli = clamp(cli, lower=bounds[1], upper=bounds[2], values=FALSE)
bounds = c(0.005, 1)
sucli = clamp(cli, lower=bounds[1], upper=bounds[2], values=FALSE)

# Plotting the function. 
pdf(paste("../../figures/", "cli_mon_", dt, ".pdf", sep = ""),  
    width = 6.5, # inches
    height = 4,
    pointsize = 7) # inches

# Note: asp will fix the white box issue but it comes at a price of 
# messing with the lat lon ratio.
par(mar=c(4, 4, 3, 7))
image(lcli,
      col = colmap[1:50],
      ylim = c(16, 39),
      xlim = c(-175, -132),
      axes = FALSE)
image(slcli,
      col = colmap[1:3],
      ylim = c(16, 40),
      xlim = c(-175, -132),
      add = TRUE,
      axes = FALSE)
image(ucli,
      col = colmap[50:100],
      ylim = c(16, 40),
      xlim = c(-175, -132),
      add = TRUE,
      axes = FALSE)
image(sucli,
      col = colmap[90:100],
      ylim = c(16, 39),
      xlim = c(-175, -130),
      add = TRUE,
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
title(expression(Climate~Signal~of~Summertime~Chlorophyll~from~1997~to~2023~(mg~m^-3~year^-1)))
title(ylab = "Latitude", line = 2)
title(xlab = "Longitude", line = 2)
image.plot(zlim = c(0, 0.005),
           smallplot= c(0.89, 0.93, 0.52, 0.91),
           axis.args = c(lwd = 0, lwd.ticks = 0.5),
           col = zmap[15:30],
           nlevel = 10,
           add = TRUE, 
           legend.only = TRUE)
image.plot(zlim = c(-0.02,0),
           axis.args = c(lwd = 0, lwd.ticks = 0.5),
           smallplot= c(0.89, 0.93, 0.12, 0.51),
           col = zmap[1:15],
           nlevel = 10,
           add = TRUE, 
           legend.only = TRUE)
box(lwd = 1, col = "black")

dev.off()

# getting frustrated and starting my own color bar.
# A = matrix(nrow = 1,
#            ncol = 50,
#            data = 1:50) 
# 
# image(A,
#       col = zmap,
#         axes = FALSE)
# axis(side = 4, 
#      las = 1, 
#      lwd = 0,
#      lwd.ticks = 1,
#      mgp = c(2, 1, 0),    
#      cex.axis = 1)
# box(lwd = 1, col = "black")






