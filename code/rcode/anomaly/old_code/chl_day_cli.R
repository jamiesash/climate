# functions --------------------------------------------------------------------
# Input is a Rasterstack object with a set Z value as a Datetime variable. 
# Output is a rasterstack with 12 raster layers  one for each month. To be 
# ploted using levelplot or mylevelplot.
# Arguments:
#    x: Raster stack with a Datetime value as the Z value. 
#       Should have at least one value for each month, otherwise an error given
bloomclim = function(x) {
  themonths <- c("January","February", "March", "April", "May","June",  "July",
                 "August", "September", "October", "November","December")
  sdate <- time(x)
  sdate <- anydate(sdate)
  m     <- months(sdate) # getting the months from the Z value
  j    = 0 # start j at 0 to count loops. i is a string
  # placeholder raster
  clim = x[[1]]

  # looping through each month and finding a spacial average using calc()
  for (i in themonths) {
    j = j + 1
    idx = which(m == i)
    z =  x[[idx]]
    z = median(z, na.rm = TRUE)
    clim = c(clim, z)
  }
  # Should add a time stamp and extent to the raster
  clim[[2:13]]
  }

anomalize = function(ras, detrend = FALSE, f = 0.6){
  themonths <- c("January","February", "March", "April", "May","June",  "July",
                 "August", "September", "October", "November","December")

  # find the monthly climotology of the data set 
  ras_clim = bloomclim(ras)
  
  # subtract each month from corresponding daily data set
  ogt = time(ras)
  ogt = anydate(ogt)
  mon_raw  = months(ogt)
  s = dim(ras_clim)
  
  # placeholder raster
  chla  = ras[[1]] 
  j    = 0 # start j at 0 to count loops. i is a string

  for (mon in themonths) {
    j = j + 1
    ind  = which(mon_raw == mon)
    temp = ras[[ind]] - ras_clim[[j]]
    chla = c(chla, temp)
  }
  
  rm(temp)
  chla = chla[[2:nlyr(chla)]]
  idx = order(time(chla))
  chla = chla[[idx]]
  
  if(detrend == TRUE) {
    clim = smooth.time.series(chla, f = f, smooth.data = TRUE)
    chla = chla - clim
    chla = setZ(chla, z = anydate(t), name = "time")
    }

  chla
}

print("1. Functions loaded.")

# ------------------------------------------------------------------------------
### Liraries and functions
library(cmocean)
library(rworldmap)
library(rworldxtra)
library(terra)
library(anytime)

print("2. Librares loaded.")

# ------------------------------------------------------------------------------
### Loading the dataset
chl = rast("/home/jamesash/blooms/data/chl_1998_2023_l3_multi_4k.nc")

print("3. Data loaded")

# ------------------------------------------------------------------------------

chla = anomalize(chl)

print("3. Anomolized")

# ------------------------------------------------------------------------------
# Change to time variable. 
t = 1:nlyr(chla)
# the first layer is the intercept, and I assume the second layer is the slope. 
cli = regress(chla, t, na.rm = TRUE, cores = 7)

print("4. Regression complete.")

# ------------------------------------------------------------------------------
### Save Raster
save = FALSE
if (save == TRUE){
writeCDF(cli,
	filename = paste("/home/jamesash/blooms/data/", "cli_mon_", dt,, sep = ""),
	overwrite = TRUE) 
        # varname = "CHL", 
        # longname="cllimatology of chl from monthly data", 
        # unit="mg/m^3", 
        # split=FALSE)
	}

# ------------------------------------------------------------------------------
### Plotting
# I may need to extract values here. I don't want to lower the suprimum.

o = sd(cli, na.rm = TRUE)
l = min(cli, na.rm = TRUE) 
u = max(cli, na.rm = TRUE) 
infi = l # + o*2
supi = u - o*2

cli = clamp(cli, lower=infi, upper=supi, values=TRUE)
print("5. Clamped.")

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

