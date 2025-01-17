# ------------------------------------------------------------------------------
### Libraries and functions
library(terra)
library(lubridate)
library(stlplus)
library(doParallel)
library(foreach)

print("1. Functions and libraries loaded.")

# Parallel STL filter function
stlfilter_parallel <- function(ras) {
  ## Initialize variables
  s <- dim(ras)  ## Dimensions
  t <- time(ras)  ## Fake datetime vector for STL input
  e <- ext(ras)  ## Raster extent
  start <- c(year(t[1]), month(t[1]))
  end <- c(year(t[length(t)]), month(t[length(t)]))
  
  ## Set up parallel backend
  num_cores <- detectCores() - 1  # Leave one core free
  cl <- makeCluster(num_cores)
  registerDoParallel(cl)
  
  ## Perform STL decomposition in parallel
  results <- foreach(i = 1:s[1], .combine = 'list', .packages = c('terra', 'lubridate', 'stlplus')) %dopar% {
    y <- array(NA, c(s[2], s[3]))  # Temporary array for this row
    
    for (j in 1:s[2]) {
      tryCatch({
        pix <- unlist(unname(ras[i, j, ]))
        pix_ts <- ts(data = pix, start = start, end = end, frequency = 12)
        peaces <- stlplus(pix_ts, s.window = "periodic")
        y[j, ] <- peaces$data$remainder  # Keep the remainder
      }, error = function(e) {
        # Skip iteration if an error occurs
      })
    }
    return(y)
  }
  
  ## Combine results into a 3D array
  y <- array(NA, c(s[1], s[2], s[3]))
  for (i in 1:s[1]) {
    y[i, , ] <- results[[i]]
  }
  
  ## Finalize raster
  y_rast <- rast(y)
  time(y_rast) <- t
  ext(y_rast) <- e
  
  ## Stop the parallel backend
  stopCluster(cl)
  
  return(y_rast)
}

print("2. Functions loaded.")

# ------------------------------------------------------------------------------
### Loading the dataset
chl <- rast("/home/jamie/projects/climate/data/chl/chl_1999_2023_day_small_l3.nc")

print("3. Data loaded")

# ------------------------------------------------------------------------------
# Apply the STL filter in parallel
clim <- stlfilter_parallel(chl)

print("4. Decomposed")

# ------------------------------------------------------------------------------
### Save the data
dt <- gsub("-", "", as.character(Sys.Date()))
writeCDF(clim, 
         filename = paste("/home/jamesash/koa_scratch/", "chla_stl_day_small", dt, ".nc", sep = ""), 
         overwrite = TRUE,
         varname = "CHL")

print("5. Saved")
