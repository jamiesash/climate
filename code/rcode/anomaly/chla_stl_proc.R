# ------------------------------------------------------------------------------
### Liraries and functions
library(lubridate)
library(stlplus)
library(zoo)
library(terra)
terraOptions(memmax = 125)
print("1. Functions loaded.")

# functions --------------------------------------------------------------------

stlfilter = function(ras) {
  # initialize an empty vector of the same size.
  s = dim(ras) 
  t = time(ras) 
  e = ext(ras)
  y = array(rep(NA, s[1]*s[2]*s[3]), c(s[1], s[2],s[3])) # the empty array
  gc()
  
  start = c(year(t[1]), month(t[1]))
  end = c(year(t[length(t)]), month(t[length(t)]))
  
  # should retrieve the rasters mask.
  na_mask = is.na(ras)
  # same approximate from stat library wrapped in terra function. 
  ras = approximate(ras, rule=2)

  gc(verbose = TRUE)

  for (i in 1:s[1]) {
    gc(verbose = TRUE)
    for (j in 1:s[2]) {
      tryCatch({
        pix = unlist(unname(ras[i,j,]))
        pix = ts(data = pix,
                 start = start,
                 end = end, 
                 frequency = 365.24)
        peaces = stlplus(pix, s.window = "periodic", s.jump = 30, l.jump = 30, t.jump = 30)
        y[i,j,] = peaces$data$remainder  #keep the seasonal
      }, error=function(e){
        #skips to next iteration if there's an error  
      }) 
    } 
  }
  gc()

  y = rast(y)
  values(y)[values(na_mask)] = NA
  crs(y) = crs(ras)
  names(y) = names(ras)
  time(y) = t
  ext(y)  = e
  y  
} 

print("2. Funcitons loaded.")

# ------------------------------------------------------------------------------
### Loading the dataset

chl = rast("/home/jamesash/climate/data/chl/chl_1999_2023_day_small_l3.nc")
#chl = rast("/home/jamesash/climate/data/chl/chl_1998_2023_l3_multi_4k.nc")

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
         filename = paste("/home/jamesash/koa_scratch/", "chla_stl_day_small", dt, ".nc",sep = ""), 
         overwrite = TRUE,
         varname = "CHL")

print("6. Saved")
