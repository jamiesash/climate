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
        y[i,j,] = peaces$data$seasonal  #keep the seasonal
      }, error=function(e){
        #skips to next iteration if there's an error  
      }) 
    } 
  }
  gc()
  y = rast(y)
  gc()
  time(y) = t
  ext(y)  = e
  y  
} 

print("2. Funcitons loaded.")

# ------------------------------------------------------------------------------
### Loading the dataset
chl = rast("/home/jamesash/climate/data/chl/chl_1998_2023_l4_month_multi_4k.nc")
print("3. Data loaded")

# ------------------------------------------------------------------------------
# apply the stl filter. 

clim = stlfilter(chl)

print("4. Decomposed.")

# ------------------------------------------------------------------------------
### the data
dt = gsub("-", "", as.character(Sys.Date()))
writeCDF(clim, 
         filename = paste("/home/jamesash/koa_scratch/", "seas_stl_mon_", dt, ".nc",sep = ""), 
         overwrite = TRUE,
         varname = "CHL")

