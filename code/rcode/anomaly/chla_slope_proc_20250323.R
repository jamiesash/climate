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

lm_ras = function(ras) {
  # Load metadata from the raster
  ras_time <- terra::time(ras)
  ras_ext <- terra::ext(ras)
  ras_crs <- terra::crs(ras)
  ras_names <- terra::names(ras)
  s <- dim(ras)
  print("Got all the metadata")
  # Time information
  start <- c(year(ras_time[1]), month(ras_time[1]))
  end <- c(year(ras_time[length(ras_time)]), month(ras_time[length(ras_time)]))
  print("Initial setup complete.")
  slope_matrix = matrix(NA, nrow = s[1], ncol = s[2]) 
  pvalue_matrix = matrix(NA, nrow = s[1], ncol = s[2]) 
  months = month(ras_time)
  month_idx = months %in% c(6,7,8,9,10)
  lm_time = 1:s[3]
  # Process each slice in a loop, keeping raster format as long as possible
  for (i in 1:s[1]) {
    slice = ras[i,, drop = FALSE]
    slice = as.matrix(slice)  # Convert only 1 slice at a time
    print(paste("Processing slice:", i))
    gc()
    for (j in 1:s[2]) {
      tryCatch({
        pix = slice[j, ]  # Extract time series
        # print("sliced")
        pix = unname(unlist(pix))
        if (all(is.na(pix))) {
          next  # Skip empty time series
        }
        na_mask = is.na(pix)  # Store NA locations
        pix = na.approx(pix, rule = 2, na.rm = FALSE)  # Fill NAs

        # Perform linear regression on the residuals
        model = lm(pix[month_idx] ~ lm_time[month_idx]) 
        model_summary = summary(model)   # Get the summary of the model
        # Extract the p-value for the slope coefficient (second coefficient)
        p_value = coef(model_summary)[2, 4]  # p-value is in the 4th column
        slope = coef(model)[2]  
        # Store the slope in the corresponding grid cell
        slope_matrix[i, j] = slope
        pvalue_matrix[i, j] = p_value      
      }, error = function(e) {
        print(paste("Skipped pixel:", i, j))
      })
    }
  }
  print("stl filtered.")
  slope_ras = rast(slope_matrix)
  pvalue_ras = rast(pvalue_matrix)
  crs(slope_ras) = ras_crs
  ext(slope_ras)  = ras_ext
  ext(pvalue_ras)  = ras_ext
  crs(pvalue_ras) = ras_crs
  c(slope_ras, pvalue_ras)
} 

# ------------------------------------------------------------------------------
### Liraries and functions
library(anytime)

library(zoo)
library(lubridate)
library(terra)
library(ncdf4)

print("2. Librares loaded.")

# ------------------------------------------------------------------------------
### Loading the dataset
chl = rast("/home/jamesash/climate/data/chl/chl_1998_2023_l3_multi_4k.nc")
# chl = rast("/home/jamesash/climate/data/chl/chl_1999_2023_day_small_l3.nc")

print("3. Data loaded")

# Get the number of time steps
ntime = nlyr(chl)
nlat  = nrow(chl)
nlon  = ncol(chl)

# Calculate half the number of time steps
half_time = ceiling(ntime / 3)
half_lat  = ceiling(nlat / 3)
half_lon  = ceiling(nlon / 3)

# Subset the raster to the first half of the time dimension
chl = chl[1:half_lat, 1:half_lon, 1:half_time]

print("3. Cut data down")

# ------------------------------------------------------------------------------

# remove the seasonal climatologic signal. 
chl = anomalize(chl)
gc()

chl = anomalize(chl)
slopes = lm_ras(chl)
print("4. anomalized")

# ------------------------------------------------------------------------------

### Save Raster
dt = gsub("-", "", as.character(Sys.Date()))

writeCDF(chl, 
	 zname="time",
	 filename = paste("/home/jamesash/koa_scratch/", "chla_day_l3_1998_2023_", dt, ".nc",sep = ""), 
	 overwrite = TRUE)

