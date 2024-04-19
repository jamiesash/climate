
# ------------------------------------------------------------------------------
### Liraries and functions
library(cmocean)
library(rworldmap)
library(rworldxtra)
library(terra)
library(anytime)

source("functions.R")

# ------------------------------------------------------------------------------
### Loading the dataset
# chl = rast("/home/jamesash/blooms/data/chl_1998_2023_l4_month_multi_4k.nc")
chl = rast("/home/jamesash/blooms/data/chl_1998_2023_l3_multi_4k.nc")

# ------------------------------------------------------------------------------

chla = anomalize(chl)

# ------------------------------------------------------------------------------
# Perform the regression right awway fug it. 
t = 1:nlyr(chla)
# the first layer is the intercept, and I assume the second layer is the slope. 
cli = regress(chla, t, na.rm = TRUE)

# -------------------------------------------------------------------------------
# using the lm function doirectly for na.rm
# t = 1:nlyr(chl)
# locreg = function(x) { 
#	# t = 1:length(x)
#        lm(x ~ t, na.action = na.exclude)$coefficients[2]
#         }
# # apply the linear regression to each grid cell. 
# slp = app(chl, locreg)

# -------------------------------------------------------------------------------
wdmap <- getMap(resolution = "high")
colmap = cmocean("delta")(100)
dt = gsub("-", "", as.character(Sys.Date()))
e = ext(cli)

# Plotting the function. 
pdf(paste("/home/jamesash/blooms/figures/", "cli_", dt, ".pdf", sep = ""),  
    width = 5.5, # inches
    height = 4,
    pointsize = 10) # inches

plot(cli[[2]], 
	# ylim = c(16, 40),
	# xlim = c(-175, -130),
	col = colmap, 
	mar = c(3.1, 3.1, 2.1, 7.1),
	plg = list(size = c(1, 1.25)),
	range = c(-0.001, 0.001),
	ylab = "Latitude",
	xlab = "Longitude")
	#breaks = 100)
plot(wdmap,
     ylim = e[3:4],
     xlim = e[1:2],
     asp = 1,
     bg = "black",
     border = "black",
     col = "black",
     add = TRUE,
     lwd = 1)
dev.off()

# ------------------------------------------------------------------------------

