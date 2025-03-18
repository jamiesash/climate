# ------------------------------------------------------------------------------
### Liraries and functions
library(lubridate)
library(stlplus)
library(zoo)
library(terra)
terraOptions(memfrac=.90)
print("1. Functions loaded.")

# functions --------------------------------------------------------------------

stlfilter = function(ras) {
  # Load metadata from the raster
  ras_time <- terra::time(ras)
  ras_ext <- terra::ext(ras)
  ras_crs <- terra::crs(ras)
  ras_names <- terra::names(ras)

  print("Got all the metadata")

  s <- dim(ras)
  print("Dimension time")

  # Create an empty terra raster for output
  # y_ras <- rast(ras)  # Copy input raster structure
  # values(y_ras) <- NA  # Set all values to NA

  print("Initialized output raster")

  # Time information
  start <- c(year(ras_time[1]), month(ras_time[1]))
  end <- c(year(ras_time[length(ras_time)]), month(ras_time[length(ras_time)]))

  print("Initial setup complete.")
  
  temp_list = list()
  # Process each slice in a loop, keeping raster format as long as possible
  for (i in 1:s[1]) {
    slice = ras[i,, drop = FALSE]
    slice = as.matrix(slice)  # Convert only 1 slice at a time
    # Allocate a temporary matrix for storing output before writing
    # temp_result = matrix(NA, nrow = s[2], ncol = s[3])
    temp_result = array(NA, c(s[2], s[3]))
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
        # print("if all na's passed")
        na_mask = is.na(pix)  # Store NA locations
        # print("na indexed")  
        pix = na.approx(pix, rule = 2, na.rm = FALSE)  # Fill NAs
        # print("approxed")
        pix_ts = ts(data = pix, start = start, end = end, frequency = 365.24)
        # print("timeseried")
        peaces = stlplus(pix_ts, s.window = "periodic", s.jump = 20, l.jump = 20, t.jump = 20)
        # print("stled")
        anom = peaces$data$remainder
        # print("anoms taken out")
        anom[na_mask] = NA  # Restore original NAs
        # print("na's back") 
        # y_ras[i,j,] = anom
        temp_result[j, ] = as.numeric(anom)  # Store in temporary matrix
        # print("anom's saved")
      }, error = function(e) {
        print(paste("Skipped pixel:", i, j))
      })
    }
    temp_ras = rast(temp_result)
    temp_list[[i]] = temp_ras
    # y_ras[[i]] = temp_ras 
    rm(temp_result, temp_ras)
    gc()
  }
  
  y_ras = subset(ras, 1)
  one_band = subset(ras, 1)
  values(one_band) = NA
  values(y_ras) = NA 
  for (i in 1:s[3]){
    for (j in 1:s[1]){
      temp_ras = temp_list[[j]]
      one_band[j,,] = temp_ras[, i, drop = TRUE]
      }
    y_ras = c(y_ras, one_band)
  }
  print("jankily way to put the rasters into the correct format.")

  # y_ras = terra::rast(temp_list)
  print("stl filtered.")
  # crs(y_ras) = ras_crs
  # names(y_ras) = ras_names
  # terra::time(y_ras) = ras_time
  # ext(y_ras)  = ras_ext
  y_ras  
} 

print("2. Funcitons loaded.")

# ------------------------------------------------------------------------------
### Loading the dataset

# chl = rast("/home/jamesash/climate/data/chl/chl_1999_2023_day_small_l3.nc")
chl = rast("/home/jamesash/climate/data/chl/chl_1998_2023_l3_multi_4k.nc")

print("3. Data loaded")

# ------------------------------------------------------------------------------
# apply the stl filter. 
clim = stlfilter(chl)
gc()
print("5. Decomposed")

# ------------------------------------------------------------------------------
### the data

dt = gsub("-", "", as.character(Sys.Date()))
writeCDF(clim, 
         filename = paste("/home/jamesash/koa_scratch/", "chla_stl_day_lage_", dt, ".nc", sep = ""), 
         overwrite = TRUE,
         varname = "CHL")

print("6. Saved")
