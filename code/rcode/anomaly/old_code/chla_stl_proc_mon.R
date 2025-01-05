# ------------------------------------------------------------------------------
### Liraries and functions
library(terra)
library(lubridate)
library(stlplus)

print("1. Functions loaded.")

# functions --------------------------------------------------------------------

stlfilter = function(ras, trend = FALSE, anom = FALSE, seas = FALSE) {
  ## initialize an empty vector of the same size.
  s = dim(ras) ## dimensions
  t = time(ras) ## fake datetime vector for stl input
  e = ext(ras)
  y = array(rep(NA, s[1]*s[2]*s[3]), c(s[1], s[2],s[3])) ## the empty array
  gc()
  
  start = c(year(t[1]), month(t[1]))
  end = c(year(t[length(t)]), month(t[length(t)]))
  
  for (i in 1:s[1]) {
    if (i == round((s[1]/5)))   print("So much more to go")
    if (i == round(s[1]/2))     print("You're half way there")
    if (i == round((s[1]/5)*4)) print("Last leg")
    for (j in 1:s[2]) {
      tryCatch({
        pix = unlist(unname(ras[i,j,]))
        pix = ts(data = pix,
                 start = start,
                 end = end, 
                 frequency = 12)
        peaces = stlplus(pix, s.window = "periodic")
        if (anom  == TRUE) y[i,j,] = peaces$data$remainder  #keep the seasonal
        if (seas  == TRUE) y[i,j,] = peaces$data$seasonal  #keep the seasonal
        if (trend == TRUE) y[i,j,] = peaces$data$trend #keep the residual
      }, error=function(e){
        #skips to next iteration if there's an error  
      }) 
    } 
    if (i == max(s[1])) print("Finish!")
  }
  gc()
  y = rast(y)
  gc()
  
  # y[na_id]  <- NA #X and Y are not same dim()
  time(y) = t
  ext(y)  = e
  y  
} 

print("2. Funcitons loaded.")

# ------------------------------------------------------------------------------
### Loading the dataset
# setwd("/home/jamie/projects/climate/code/rcode/anomaly")
chl = rast("/home/jamesash/climate/data/chl/chl_1998_2023_l4_month_multi_4k.nc")

print("3. Data loaded")

# ------------------------------------------------------------------------------
# removing coastal effects

# subset spacial polygon like a data frame. 
eez = terra::vect("/home/jamesash/climate/data/eez/USMaritimeLimitsNBoundaries.shp")
idx = which(eez$REGION == 'Hawaiian Islands')
hawaii = eez[idx,]
idx = hawaii$CZ + hawaii$TS
idx = which(idx == 1)
hawaii = hawaii[idx,]

# mask land.
# chl = terra::mask(chl, hawaii, inverse = TRUE)
# gc()

print("3. Masked")

# ------------------------------------------------------------------------------
# apply the stl filter. 

clim = stlfilter(chl, anom = TRUE)

print("4. Decomposed.")

# ------------------------------------------------------------------------------
### the data
dt = gsub("-", "", as.character(Sys.Date()))
writeCDF(clim, 
         filename = paste("/home/jamesash/koa_scratch/", "chla_stl_mon_", dt, ".nc",sep = ""), 
         overwrite = TRUE,
         varname = "CHL")

