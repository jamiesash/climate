# ------------------------------------------------------------------------------
### Liraries and functions
library(terra)
library(lubridate)
library(stlplus)

print("1. Functions loaded.")

# functions --------------------------------------------------------------------

stlfilter = function(ras) {
  ## initialize an empty vector of the same size.
  s = dim(ras) ## dimensions
  t = time(ras) ## fake datetime vector for stl input
  e = ext(ras)
  y = array(rep(NA, s[1]*s[2]*s[3]), c(s[1], s[2],s[3])) ## the empty array
  gc()
  
  start = c(year(t[1]), month(t[1]))
  end = c(year(t[length(t)]), month(t[length(t)]))
  
  for (i in 1:s[1]) {
    for (j in 1:s[2]) {
      tryCatch({
        pix = unlist(unname(ras[i,j,]))
        pix = ts(data = pix,
                 start = start,
                 end = end, 
                 frequency = 12)
        peaces = stlplus(pix, s.window = "periodic")
        y[i,j,] = peaces$data$remainder  #keep the seasonal
      }, error=function(e){
        #skips to next iteration if there's an error  
      }) 
    } 
    if (i == max(s[1])) print("Finish!")
  }
  gc()
  y = rast(y)
  # y[na_id]  <- NA #X and Y are not same dim()
  time(y) = t
  ext(y)  = e
  y  
} 

print("2. Funcitons loaded.")

# ------------------------------------------------------------------------------
### Loading the dataset
# setwd("/home/jamie/projects/climate/code/rcode/anomaly")
chl = rast("/home/jamesash/climate/data/chl/chl_1999_2024_small_daily_multi_l3_4k.nc")

print("3. Data loaded")

# ------------------------------------------------------------------------------
# apply the stl filter. 

clim = stlfilter(chl)

print("4. Decomposed")

# ------------------------------------------------------------------------------
### the data

dt = gsub("-", "", as.character(Sys.Date()))
writeCDF(clim, 
         filename = paste("/home/jamesash/koa_scratch/", "chla_stl_day_small", dt, ".nc",sep = ""), 
         overwrite = TRUE,
         varname = "CHL")

print("5. Saved")
