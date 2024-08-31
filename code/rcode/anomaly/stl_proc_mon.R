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
        gc()
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

print("1. Functions loaded.")

# -----------------------------------------------------------------------------
### testing
# t = time(ras)
# chlu = global(ras, "mean", na.rm=TRUE)
# chlu = unlist(unname(chlu))
# start = c(year(t[1]), month(t[1]))
# end = c(year(t[length(t)]), month(t[length(t)]))
# 
# chl_dec = ts(data = chlu,
#              start = start,
#              end = end,
#              frequency = 12)
# peaces = stl(chl_dec, s.window = "periodic", s.degree = 0, t.degree = 0, l.degree=0, robust=TRUE)
# gc()
# 
# test = chl_dec - peaces$time.series[,2]  
# plot(test)

# ------------------------------------------------------------------------------
### Liraries and functions
library(anytime)
library(terra)
library(lubridate)
library(stlplus)

print("2. Librares loaded.")

# ------------------------------------------------------------------------------
### Loading the dataset
setwd("/home/jamie/projects/climate/code/rcode/anomaly")
chl = rast("../../../data/chl/chl_1998_2023_l4_month_multi_4k.nc")

print("3. Data loaded")

# Subset smaller area.

# ------------------------------------------------------------------------------
# removing coastal effects

# subset spacial polygon like a data frame. 
eez = terra::vect("../../../data/eez/USMaritimeLimitsNBoundaries.shp")
idx = which(eez$REGION == 'Hawaiian Islands')
hawaii = eez[idx,]
idx = hawaii$CZ + hawaii$TS
idx = which(idx == 1)
hawaii = hawaii[idx,]

# mask land.
chl = terra::mask(chl, hawaii, inverse = TRUE)
gc()

# ------------------------------------------------------------------------------
# remove the seasonal climatologist signal. 

# calculate the climate trend
clim = stlfilter(chl, trend = TRUE)

clim_map = mean(clim, na.rm = TRUE)
clim_ts = global(clim, "mean", na.rm = TRUE)
clim_ts = unlist(unname(clim_ts))
gc()

print("4. Decomposed.")

t = time(chl)

df = data.frame(t, clim_ts)
colnames(df) = c("time", "clim")

dt = gsub("-", "", as.character(Sys.Date()))

# ------------------------------------------------------------------------------
### the data

write.csv(df, paste("/home/jamie/projects/climate/data/chl/", "stl_ts_mon_", dt, ".csv", sep = ""), row.names = FALSE)

writeCDF(clim, 
         filename = paste("/home/jamie/projects/climate/data/chl/", "stl_mon_", dt, ".nc",sep = ""), 
         overwrite = TRUE,
         varname = "CHL")

