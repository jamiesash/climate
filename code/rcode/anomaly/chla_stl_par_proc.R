# ------------------------------------------------------------------------------
### Libraries and functions
library(zoo)
library(terra)
library(lubridate)
library(stlplus)
library(doParallel)
library(foreach)
terraOptions(memfrac=0.97)  # mmemfrac=.98)

print("1. Functions and libraries loaded.")

# Parallel STL filter function
stlfilter_parallel = function(ras) {
  
  ## Initialize variables
  s = dim(ras) 
  ras_time = time(ras)  
  ras_ext = ext(ras)
  ras_crs = crs(ras) 
  ras_names = names(ras)
  start = c(year(ras_time[1]), month(ras_time[1]))
  end = c(year(ras_time[length(ras_time)]), month(ras_time[length(ras_time)]))
  
  ras = array(ras, s) # DoParrellel does not like rasters
  slices = list()
  for (i in 1:s[1]) slices[[i]] = ras[i,,]
  rm(ras)
  gc()
  print("converted to array and sliced.")

  ## Set up parallel backend
  print("number of cores is...")
  cores = detectCores()
  print(as.character(cores))
  
  print("requesting 10 cores")
  cl = makeCluster(10)
  registerDoParallel(cl)
  
  print("clustered.")

  ## Perform STL decomposition in parallel
  results <- foreach(slice = slices, .packages = c('lubridate', 'stlplus', 'zoo')) %dopar% {
    y = array(NA, c(s[2], s[3]))  # Temporary array for this row
    for (j in 1:s[2]) {
      tryCatch({
        pix = slice[j, ] 
        na_mask = is.na(pix)
        pix = na.approx(pix, rule = 2, na.rm = FALSE)
	pix_ts = ts(data = pix, start = start, end = end, frequency = 365.24)
        # use s.jump, t.jump, l.jump for speed. will let the window jump rather than one every step. 
        peaces = stlplus(pix_ts, s.window = "periodic", s.jump = 30, l.jump = 30, t.jump = 30)
        anom = peaces$data$remainder
        anom[na_mask] = NA
        y[j, ] = anom  
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
  # get back the extent and time
  output_array = rast(output_array)
  terra::time(output_array) = ras_time
  crs(output_array) = ras_crs
  names(output_array) = ras_names
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
