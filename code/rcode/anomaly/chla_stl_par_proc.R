# ------------------------------------------------------------------------------
### Libraries and functions
library(zoo)
library(terra)
library(lubridate)
library(stlplus)
library(doParallel)
library(foreach)
terraOptions(memmax = 125)  # mmemfrac=.98)

print("1. Functions and libraries loaded.")

# Parallel STL filter function
stlfilter_parallel = function(ras) {
  
  ## Initialize variables
  s = dim(ras) 
  ras_time = time(ras)  
  ras_ext = ext(ras) 
  start = c(year(ras_time[1]), month(ras_time[1]))
  end = c(year(ras_time[length(ras_time)]), month(ras_time[length(ras_time)]))
  
  ras = array(ras, s) # DoParrellel does not like rasters
  gc(verbose = TRUE)
  
  print("converted to array.")
  
  # Store NA indices and interpolate only missing values before the loop
  na_masks = array(FALSE, dim(ras))  # Track original NA locations
  
  for (i in 1:s[1]) {
    for (j in 1:s[2]) {
      pix = ras[i, j, ]
      na_masks[i, j, ] = is.na(pix)  # Store NA locations
      ras[i, j, ] = na.approx(pix, rule = 2, na.rm = FALSE)  # Interpolate
    }
  }
  gc(verbose = TRUE)

  print("Nan's filled, and nanmask found.")
  
  ## Set up parallel backend
  print("number of cores is...")
  detectCores()
  cl = makeCluster(detectCores())
  registerDoParallel(cl)
  
  print("clustered.")

  ## Perform STL decomposition in parallel
  results <- foreach(slice = iter(ras, by = "row"), .packages = c('terra', 'lubridate', 'stlplus', 'zoo')) %dopar% {
    y = array(NA, c(s[2], s[3]))  # Temporary array for this row
    
    print("inside the clustered.")
    for (j in 1:s[2]) {
      tryCatch({
        pix = unlist(unname(slice[j, ])) 
	pix_ts = ts(data = pix, start = start, end = end, frequency = 365.24)
        # use s.jump, t.jump, l.jump for speed. will let the window jump rather than one every step. 
        peaces = stlplus(pix_ts, s.window = "periodic", s.jump = 20, l.jump = 20, t.jump = 20)
        y[j, ] = peaces$data$remainder  # Keep the remainder
        }, error=function(e){
          #skips to next iteration if there's an error  
          }) 
        }
    y
  }
    
  print("loop ended.")

  # Convert list of results back to array
  output_array = array(NA, c(s[1], s[2], s[3]))
  for (i in 1:s[1]) {
    output_array[i, , ] = results[[i]]
  }
   
  ## Stop the parallel backend
  stopCluster(cl)
  gc(verbose = TRUE)  
  
  print("cluster stopped.")
  
  # Restore NA values to original positions
  output_array[na_masks] = NA
  # get back the extent and time
  output_array = rast(output_array)
  time(output_array) = ras_time
  ext(output_array) = ras_ext
  gc()

  output_array
}

print("2. Functions loaded.")

# ------------------------------------------------------------------------------
### Loading the dataset
chl <- rast("/home/jamesash/climate/data/chl/chl_1999_2023_day_small_l3.nc")
# chl <- rast("/home/jamesash/climate/data/chl/chl_1998_2023_l3_multi_4k.nc")
# chl <- rast("/home/jamesash/climate/data/chl/chl_1998_2023_l4_month_multi_4k.nc")
print("3. Data loaded")

# ------------------------------------------------------------------------------
# Apply the STL filter in parallel
clim <- stlfilter_parallel(chl)

print("5. Decomposed")

# ------------------------------------------------------------------------------
### Save the data
dt <- gsub("-", "", as.character(Sys.Date()))
writeCDF(clim, 
         filename = paste("/home/jamesash/koa_scratch/", "chla_stl_day_small", dt, ".nc", sep = ""), 
         overwrite = TRUE,
         varname = "CHL")

print("5. Saved")
