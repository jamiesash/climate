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
  ras_time <- time(ras)  ## Fake datetime vector for STL input
  ras_ext <- ext(ras)  ## Raster extent
  start <- c(year(ras_time[1]), month(ras_time[1]))
  end <- c(year(ras_time[length(ras_time)]), month(ras_time[length(ras_time)]))
  
  ras = focal(ras, w=21, fun=mean, na.policy="only", na.rm=T)
  ras = array(ras, s) # convert to array. DoParrellel does not like
  gc()
  
  ## Set up parallel backend
  num_cores <- detectCores() - 6  # Leave one core free
  cl <- makeCluster(num_cores)
  registerDoParallel(cl)
  
  ## Perform STL decomposition in parallel
  results <- foreach(i = 1:s[1], .packages = c('terra', 'lubridate', 'stlplus')) %dopar% {
    y <- array(NA, c(s[2], s[3]))  # Temporary array for this row
    
    for (j in 1:s[2]) {
      tryCatch({
        pix <- unlist(unname(ras[i, j, ])) 
        pix_ts <- ts(data = pix, start = start, end = end, frequency = 365.24)
        # setting the window to 13 should speed this up. Right now it estimates the window. 
        # t.window: harder to know. it is calculated as nextodd(ceiling((1.5*period) / (1-(1.5/s.window))))
        # use s.jump, t.jump, l.jump for speed. will let the window jump rather than one every step. 
        peaces <- stlplus(pix_ts, s.window = "periodic", s.jump = 20, t.jump = 20)
        y[j, ] <- peaces$data$remainder  # Keep the remainder
        }, error=function(e){
          #skips to next iteration if there's an error  
          }) 
        }
    y
  }
  
  ## Combine results into a 3D array
  y <- array(NA, c(s[1], s[2], s[3]))
  for (i in 1:s[1]) {
    y[i, , ] <- results[[i]]
  }
  
  ## Finalize raster
  y_rast <- rast(y)
  time(y_rast) <- ras_time
  ext(y_rast) <- ras_ext
  
  ## Stop the parallel backend
  stopCluster(cl)
  
  y_rast
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