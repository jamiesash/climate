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


# The terra version of anomalize. I think. 
anomalize = function(ras){
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
  
  chla
}

subsum    = function(x, mnths = 7:10) {
  sdate     <- time(x)
  themonths <- c("January","February",
                 "March", "April",
                 "May","June",
                 "July","August",
                 "September", "October",
                 "November","December")
  m     <- months(sdate)
  idx_t <- which(is.element(m, themonths[mnths]))
  x     <- subset(x, idx_t)
  sdate <- sdate[idx_t]
  time(x) = sdate
  x
}
# ------------------------------------------------------------------------------
### Liraries and functions
library(anytime)
library(terra)

print("1. Librares loaded.")

# ------------------------------------------------------------------------------
### Loading the dataset
# chl = rast("/home/jamesash/climate/data/chl_1998_2023_l4_month_multi_4k.nc")
chl = rast("/home/jamesash/climate/data/chl/chl_1998_2023_l4_month_multi_4k.nc")

print("2. Data loaded")

# ------------------------------------------------------------------------------
# remove the seasonal climatologic signal. 
chl = anomalize(chl)
gc()

# remove all but the summer months.
chl = subsum(chl, mnths = 6:10)
gc()

print("3. Anamolized and summer subset.")

# ------------------------------------------------------------------------------
# Perform find the climatology. 
cli = app(chl, fun = mean, na.rm = TRUE)

print("4. Averaged.")

# ------------------------------------------------------------------------------

### Save Raster
writeCDF(cli, 
	filename = paste("/mnt/lustre/koa/scratch/jamesash/", "climatology_mon_sum", ".nc", sep = ""), 
	overwrite = TRUE,
	varname = "CHL")

