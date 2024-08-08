# functions --------------------------------------------------------------------

# input is a raster output is a boolian raster of bloom/not bloom
bool = function(x){
  u = app(x, fun = median, na.rm = TRUE)
  o = app(x, fun = sd, na.rm = TRUE)
  boo = x > (u + o)
  ext(boo) = ext(x)
  time(boo) = time(x)
  boo
}

# ------------------------------------------------------------------------------
### Liraries and functions
library(anytime)
library(terra)

print("1. Librares loaded.")

# ------------------------------------------------------------------------------
### Loading the dataset
chl = rast("/home/jamesash/climate/data/chl/chla_day_l3_2017_2023_20240804.nc")

print("2. Data loaded")

# ------------------------------------------------------------------------------

# remove the seasonal climatologic signal. 
chl = bool(chl)
gc()

print("3. Blooms found")

# ------------------------------------------------------------------------------

### Save Raster
dt = gsub("-", "", as.character(Sys.Date()))

writeCDF(chl, 
	filename = paste("/home/jamesash/climate/data/chl/", "blooms_day_l3_2017_2023_", dt, ".nc",sep = ""), 
	overwrite = TRUE,
	varname = "blooms")

