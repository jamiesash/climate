# this script is really old ad could use some work. left a few legacy functions
# at teh top so that it will still run.

source("code\\libraries.R")
source("code\\functions.R")
source("code\\blooms.R")
rasterOptions(maxmemory = 120e+10, memfrac = 0.9)
knitr::opts_chunk$set(echo = TRUE)

# using a two degree box around aloha
lons = c(-158, -130)
lats = c(22, 35)
e = extent(lons, lats)

# low functions ----------------------------------------------------------------
# input is a raster output is a boolian raster of bloom/not bloom
bool = function(x){
  u <- calc(x, fun = median, na.rm = TRUE)
  o <- calc(x, fun = mad, na.rm = TRUE)
  boo <- x > (u + o)
  extent(boo) <- extent(x)
  boo <- setZ(boo, z = getZ(x), name = "time")
  boo
}

# Input is 3D chla raster output is time series vector of bloom area
barea = function(b){
  t   = getZ(b)
  #boo = bool(x)
  rr  = reclassify(b, cbind(0,NA))
  a   = raster::area(rr)
  #aa  = a * boo
  temp <- raster::mask(a, rr)
  km2 = cellStats(temp, stat = "sum", na.rm = TRUE)
  data.frame(area = as.numeric(km2), time= t)
}

# input x is a raster of raw chl
# inpur x is a raster of 0/1 bloom non-bloom
bmag <- function(x, b){
  t = getZ(x)
  # create a raster of bloom 0, 1's as b
  # boo = bool(x)
  # temp = x * boo
  
  # Area time series for duration 
  rr <- reclassify(b, cbind(0,NA))
  # find area of each cell + mask cells that matter
  temp <- raster::mask(x, rr)
  
  mag  = cellStats(temp, stat ="max", na.rm = TRUE)
  data.frame(mag = as.numeric(mag), time= t)
}  

vectorize  <- function(x) {
  sdate <- getZ(x)
  # x = raster::flip(x, direction = "x")
  # x = raster::flip(x, direction = "y")
  # 
  # e = extent(x)
  # x = t(x)
  # extent(x) = e
  # 
  x     <- rasterToPoints(x)
  x     <- data.frame(x)
  colnames(x) <- c("lon", "lat", as.character(sdate))
  x     <- reshape2::melt(x, id.vars = c("lat", "lon"))
  colnames(x) <- c("lats", "lons", "time", "val")
  x
}

fronts = function(ras,
                  downsize = 10,
                  I = 2){
  # ras  = raster::flip(ras, direction = "x")
  thresh = median(ras, na.rm = TRUE) + mad(ras, na.rm = TRUE)*I
  idx  = which(values(ras) > thresh)
  e    = extent(ras)
  s    = dim(ras)
  lons = seq(from = e[1], to = e[2], length = s[1]) 
  lats = seq(from = e[3], to = e[4], length = s[2]) 
  grid = pracma::meshgrid(lons, lats)
  fronts = data.frame(lons = grid$X[idx], lats = grid$Y[idx], value = ras[idx])
  idx = seq(1, nrow(fronts), downsize)
  fronts = fronts[idx, ]
  fronts
}

bcent <- function(ras){
  # create a raster of bloom 0, 1's as b
  pix  <- cellStats(ras, stat = mean, na.rm = TRUE)
  idx  <- as.numeric(which(pix == max(pix, na.rm = TRUE)))
  # temp = raster::mask(subset(ras, ind), subset(b, ind))
  
  # chla layer with largest bloom area
  temp = subset(ras, idx)
  # smooth that layer to reduce costal effects
  temp = raster.gaussian.smooth(x = temp, n = 21, sigma = 5)
  #temp = oreant(temp, t1 = TRUE)
  
  temp = vectorize(temp)
  idx = which(temp$val == max(temp$val, na.rm = TRUE))
  idx = max(idx)
  temp[idx,]
}

mymax = function(x, par = "area"){
  idx = which(x[, par] == max(x[, par], na.rm = TRUE))
  x = x[idx, ]
  x[1, par]
}

# connects blooms separated by desired duration
fillblips = function(x, dur = 2){
  blooms = rle(x)
  # index runs sorter than the desired duration
  shortblooms = blooms$lengths <= dur
  # index short that are false false runs
  deadair = shortblooms & !blooms$values
  # replace shorts with true
  blooms$values[deadair] = TRUE
  rep(blooms$values, blooms$lengths)
}

# removes all blooms shorter than the desired duration
# determines if element is a member of consecutive TRUE/FALSE
longblooms = function(x, dur = 2) {
  blooms    = rle(x)
  bigblooms = blooms$lengths > dur
  #bigblooms = blooms$values == TRUE
  bigblooms = rep(bigblooms, blooms$lengths)
  # TRUE * FALSE = FALSE so this filters the non blooms
  bigblooms = bigblooms & x
  #bigblooms = bigblooms == 1
  bigblooms
}

# input is a vector of true/false bloom/not bloom values
# output is a data frame of start and end indices of bloom times
bloomtails = function(x){
  # Compute endpoints of run
  blooms = rle(x)
  # subsetting only trues
  #truelengths = blooms$lengths[blooms$values]
  
  # subset trues after the fact
  end = cumsum(blooms$lengths)
  start = c(1, lag(end)[-1] + 1)
  
  data.frame(start, end)[blooms$values, ]
} 

# x is a data frame of bloom start and end indexes from bloom tails function
# ras is a raster brick with a z as the time variable
idxbloom = function(x, ras){
  subset(ras, x[1]:x[2])
}

idxbloom = function(x, ras){
  t = getZ(ras)
  y = subset(ras, x[1]:x[2])
  t = t[x[1]:x[2]]
  setZ(y, z = t, name = "time")
}

# high functions ---------------------------------------------------------------
filt_stl = function(ras = chl, res = 8, sub.start = 2){
  #res = 8 # 8 day resolution
  # needs at least 2 years to run
  ts = cellStats(ras, stat = "mean", na.rm = TRUE)
  ts = as.numeric(ts)
  
  so = stlplus(x   = ts,      # One time series of the raster
               t   = getZ(ras),       # datetime vector
               n.p = floor(365/res),      # give the period of seasonality
               s.window  = floor(30/res), # length of window: about a month 
               sub.start = sub.start, # if data does not start on XXXX-01-01
               outer     = 4) # idk something with the iterations
  
  raw  = so$data$raw
  anom = so$data$remainder
  clim = so$data$trend
  seas = so$data$seasonal
  seas_clim =  seas + clim
  
  data.frame(raw = raw, anom = anom, clim = clim, seas = seas, seas_clim =  seas + clim)
}

bloomspan = function(anom, time = getZ(chl), gaps = 0, dur = 1){
  # apply loess to 
  
  u   = median(anom, na.rm = TRUE)
  o   = mad(anom, na.rm = TRUE)
  boo = anom > u + o/2
  boo[is.na(boo)] = FALSE
  
  # Instead of filling blips, apply loess
  #boo = fillblips(boo, dur = 4)
  boo = fillblips(boo, dur = gaps)
  boo = longblooms(boo, dur = dur)
  
  bloom_idx = bloomtails(boo)
  
  dablooms  = data.frame(sdate = time[bloom_idx$start], 
                         edate = time[bloom_idx$end], 
                         duration = abs(time[bloom_idx$start] - time[bloom_idx$end]),
                         start = bloom_idx$start,
                         end = bloom_idx$end)
  idx_summer = month(dablooms$sdate) > 4 & month(dablooms$sdate) < 10 
  #bloom_idx  = bloom_idx[idx_summer, ]
  dablooms   = dablooms[idx_summer, ]
  dablooms
}

sumstats = function(blooms = chlb, signal1 = chl, signal2 = chla, sdate, edate){
  
  ts = as.Date(sdate)
  te = as.Date(edate)
  tc = getZ(signal1)
  ta = getZ(signal2)
  
  nearest = function(x, y) which.min(abs(y - x))
  
  start = list()
  end   = list()
  for(i in 1:length(ts)) start[[i]] = as.Date(ts[i])
  for(i in 1:length(te)) end[[i]] = as.Date(te[i])
  
  start = lapply(X = start, FUN = nearest, y = tc)
  end   = lapply(X = end, FUN = nearest, y = tc) 
  
  idx = data.frame(start = unlist(start), end = unlist(end))
  #idx = tbl[,c("start", "end")]
  
  cb = apply(idx, MARGIN = 1, FUN = idxbloom, ras = blooms)
  
  # instead of calculating a 0/1 I can vectorize the blooms and turn them into 
  # long-format tables
  
  co = apply(idx, MARGIN = 1, FUN = idxbloom, ras = signal1)
  
  ca = apply(idx, MARGIN = 1, FUN = idxbloom, ras = signal2)
  
  # then apply max function to list of time series to get value
  #cents = mapply(bcent, ras = co, SIMPLIFY = FALSE)
  cents = mapply(bcent, ras = ca, SIMPLIFY = FALSE)
  cents = do.call(rbind, cents)
  
  mags  = mapply(bmag, x = co, b = cb, SIMPLIFY = FALSE)
  areas = lapply(cb,    FUN = barea)
  
  m     = lapply(mags,  FUN = mymax, par = "mag")
  a     = lapply(areas, FUN = mymax, par = "area")
  
  duration = as.Date(edate) - as.Date(sdate)
  
  sumtbl = data.frame(sdate    = as.Date(sdate),
                      edate    = as.Date(edate),
                      duration = duration,
                      mag      = round(unlist(m), 2), 
                      area     = round(unlist(a), 2),
                      lat      = cents[,1],
                      lon      = cents[,2])
  
  # idx = sumtbl$duration > 16
  # sumtbl = sumtbl[idx, ]
  sumtbl
}

# Fast furiuer series decomposiiton summing first n components
nff = function(x = NULL, n = NULL, up = 1L){
  #The direct transformation
  #The first frequency is DC, the rest are duplicated
  dff = fft(x)
  #The time
  t = seq(from = 1, to = length(x))
  #Upsampled time
  nt = seq(from = 1, to = length(x)+1-1/up, by = 1/up)
  #New spectrum
  ndff = array(data = 0, dim = c(length(nt), 1L))
  ndff[1] = dff[1] #Always, it's the DC component
  if(n != 0){
    ndff[2:(n+1)] = dff[2:(n+1)] #The positive frequencies always come first
    #The negative ones are trickier
    ndff[length(ndff):(length(ndff) - n + 1)] = dff[length(x):(length(x) - n + 1)]
  }
  #The inverses
  indff = fft(ndff/73, inverse = TRUE)
  idff = fft(dff/73, inverse = TRUE)
  ret = data.frame(time = nt, signal = Mod(indff))
  return(ret)
}

# LOAD CHL DATA ----------------------------------------------------------------

# Set variables
lat_varid = "lat"
lon_varid = "lon"
sdate = as.Date("1997-01-01")
edate = as.Date("2022-12-31")
var = "CHL"
# smaller data set for testing. 
# url = "https://jash:5.Pellegrino@my.cmems-du.eu/thredds/dodsC/cmems_obs-oc_glo_bgc-plankton_my_l4-multi-4km_P1M?" 
url = "https://jash:5.Pellegrino@my.cmems-du.eu/thredds/dodsC/cmems_obs-oc_glo_bgc-plankton_my_l3-multi-4km_P1D?"
origin = "1900-01-01"

data = nc_open(url, verbose = FALSE, write = FALSE)
lat  = ncvar_get(data, varid = lat_varid)
lon  = ncvar_get(data, varid = lon_varid)
time = ncvar_get(data, varid = "time")
time = as.Date(time, origin = origin)

idx_lat <- which(lat > lats[1] & lat < lats[2])
idx_lon <- which(lon > lons[1] & lon < lons[2])
idx_time <- which(time >= sdate & time <= edate)


idx_ras <- paste("CHL",
                 paste("[", range(idx_time)[1], ":1:", range(idx_time)[2], "]", sep = ""),
                 paste("[", range(idx_lat)[1],  ":1:", range(idx_lat)[2],  "]", sep = ""),
                 paste("[", range(idx_lon)[1],  ":1:", range(idx_lon)[2],  "]", sep = ""),
                 sep = "")

idx_time <- paste("time", paste("[", range(idx_time)[1], ":1:",range(idx_time)[2], "]", sep = ""), sep = "")
idx_lat  <- paste(lat_varid, paste("[", range(idx_lat)[1], ":1:",range(idx_lat)[2],  "]", sep = ""), sep = "")
idx_lon  <- paste(lon_varid, paste("[", range(idx_lon)[1], ":1:",range(idx_lon)[2],  "]", sep = ""), sep = "")
idx <- paste(idx_lat, idx_lon, idx_time, idx_ras, sep = ",")

url <- paste(url, idx, sep = "")

nc_close(data)
rm(data)

data = nc_open(url, verbose = FALSE, write = FALSE)

lat  = ncvar_get(data, varid = lat_varid)
lon  = ncvar_get(data, varid = lon_varid)
time = ncvar_get(data, varid = "time")
time = as.Date(time, origin = origin)
ras  = ncvar_get(data)
nc_close(data)

s = dim(ras)
ras = raster::brick(ras)
ras = t(ras)
ras = setZ(ras, z = as.Date(time, origin = org), name = "time")
extent(ras) = extent(min(lons), max(lons), min(lats), max(lats))

chl = ras
rm(ras)
gc()

# ------------------------------------------------------------------------------
# Remove South of Hawaii and Chlorophyll front 
# I want all that gone
chl = polymask(chl)

gc()

# Calculate CHL anomaly and boolian bloom array
# ------------------------------------------------------------------------------

# I receive an warning when detrending for cells of all NA values (Hawaii)
# chl anom. calc. needs at least one full year to work correctly
chl_a = anomalize(chl, detrend = TRUE)
chl_b = bool(chl_a)
gc()

# ------------------------------------------------------------------------------
# Apply functions and catinate summary table
# I need to find the start and end index of each bloom
sumtbl = sumstats(blooms = chl_b, 
                  signal1 = chl, 
                  signal2 = chl_a, 
                  sdate = blooms$sdate, 
                  edate = blooms$edate)
gc()

# ------------------------------------------------------------------------------
# saving tables as csv
write.csv(sumtbl, 
          file = paste("data/summary_statistics/full_sum_", Sys.Date(), ".csv", sep = ""))

















