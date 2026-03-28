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
  s <- dim(ras)
  print("Got all the metadata")
  # Time information
  start <- c(year(ras_time[1]), month(ras_time[1]))
  end <- c(year(ras_time[length(ras_time)]), month(ras_time[length(ras_time)]))
  print("Initial setup complete.")
 
  u_matrix = matrix(NA, nrow = s[1], ncol = s[2]) 
  months = month(ras_time)
  month_idx = months %in% c(6,7,8,9,10)
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
        # print("if all na's passed")
        na_mask = is.na(pix)  # Store NA locations
        # print("na indexed")  
        pix = na.approx(pix, rule = 2, na.rm = FALSE)  # Fill NAs
        # print("approxed")
        pix_ts = ts(data = pix, start = start, end = end, frequency = 365.24)
        # print("timeseried")
        peaces = stlplus(pix_ts, s.window = "periodic", s.jump = 3, l.jump = 3, t.jump = 3)
        # print("stled")
        trend = peaces$data$trend

        trend[na_mask] = NA  # Restore original NAs
        trend = trend[month_idx]
        u = mean(trend, na.rm = TRUE)
        
        # Store the slope in the corresponding grid cell
        u_matrix[i, j] = u

      }, error = function(e) {
        print(paste("Skipped pixel:", i, j))
      })
    }
  }

  print("stl filtered.")
  u_ras = rast(u_matrix)
  crs(u_ras) = ras_crs
  ext(u_ras)  = ras_ext

  u_ras
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
print("5. Decomposed")

# ------------------------------------------------------------------------------
### the data

dt = gsub("-", "", as.character(Sys.Date()))
writeCDF(clim, 
         filename = paste("/home/jamesash/koa_scratch/", "clim_sum_stl_mean_", dt, ".nc", sep = ""), 
         overwrite = TRUE,
         varname = "CHL")

print("6. Saved")
