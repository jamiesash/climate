
#load functions and libraries
source("code\\functions.R")
source("code\\libraries.R")
source("code\\blooms.R")
rasterOptions(maxmemory = 120e+10, memfrac = 0.9)
knitr::opts_chunk$set(echo = TRUE)

# ------------------------------------------------------------------------------
# I added a few funcions here because its an older script and I changed the 
# functions in my lararies R file

# input is a boolian raster layer of 0/1 
# output is a data frame with only high values included
# Downlsize reduces the size of the data set
# I may like to simultaniously index the chl signal as well
chl_fronts = function(blooms,
                      signal,
                      downsize = 10){
  ras  = raster::flip(blooms, direction = "x")
  
  idx  = which(values(blooms) == TRUE) # == TRUE
  
  e    = extent(blooms)
  s    = dim(blooms)
  lons = seq(from = e[1], to = e[2], length = s[1]) 
  lats = seq(from = e[3], to = e[4], length = s[2]) 
  grid = pracma::meshgrid(lons, lats)
  # using the blooms index to subset the signal values
  fronts = data.frame(lons = grid$X[idx], lats = grid$Y[idx], value = signal[idx])
  idx = seq(1, nrow(fronts), downsize)
  fronts = fronts[idx, ]
  fronts
}

# use the apply function to apply the calc mean function for every 8 layers

avestack = function(ras, by = 8, fun){
  t = getZ(ras)
  s      = dim(ras)
  idx_8d = rep(1:ceiling(s[3]/8), 8)
  idx_8d = sort(idx_8d)
  rem    = abs(s[3] - length(idx_8d))
  idx_8d = idx_8d[1:(length(idx_8d)-rem)]
  idx_og = 1:s[3]
  t_idx  = unique(idx_8d)
  idx    = data.frame(idx_8d = idx_8d, idx_og = idx_og)
  idx$idx_8d = as.factor(idx$idx_8d)
  idx = split(idx, idx_8d)
  
  # subset and average every eight days
  tempset = function(x, ras) ras[[x$idx_og]]
  days = lapply(idx, FUN = tempset, ras = ras)
  days = lapply(days, FUN = calc, fun = fun, na.rm = TRUE)
  days = brick(days)
  days = setZ(days, z = t[t_idx], name = "time")
  days
}

# Writing this as a function to take a date range, a chl anomaly signal, 
# and output an area time-series

# Input is 3D chla raster output is time series vector of bloom area
bloom_area = function(b){
  t    = getZ(b)
  rr   = reclassify(b, cbind(0,NA))
  a    = raster::area(rr)
  temp = raster::mask(a, rr)
  
  km2  = cellStats(temp, stat = "sum", na.rm = TRUE)
  data.frame(area = as.numeric(km2), time= t)
}

# input x is a raster of raw chl
# inpur x is a raster of 0/1 bloom non-bloom
bloom_con = function(x, b, stat = "mean"){
  t = getZ(x)
  # make bloom raster a mask
  rr   = reclassify(b, cbind(0,NA))
  # mask bloom raster with bloom raster
  temp = raster::mask(x, rr)
  mag  = cellStats(temp, stat = stat, na.rm = TRUE)
  data.frame(con = as.numeric(mag), time= t)
}  

basic_theme <- function(x,
                        y, 
                        mar = c(1,1,1,1),
                        #asp = 1,
                        ylim = range(y),
                        xlim = range(x),
                        main = "",
                        dt = TRUE,
                        ylab = "",
                        xlab = ""
){
  if(is.numeric(x)) x <- round(x, 1)
  if(is.numeric(y)) y <- round(y, 1)
  par(mar = mar)
  plot(0,
       ylim = ylim,
       xlim = xlim,
       main = "",
       xlab = "",
       ylab = "",
       axes = FALSE,
  )
  grid(nx = NULL, # X-axis divided in two sections
       ny = NULL, # Y-axis divided in three sections
       lty = 2, 
       col = colvect(c("gray69"), alpha = 0.5), lwd = 1)
  box(which = "plot", 
      lty = "solid", 
      lwd = 3, 
      col = colvect("grey22", alpha = 0.9))
  title(main = main,
        cex.lab = 2,
        line= 0.75,
        adj = 0)
  title(ylab = ylab, cex.lab = 1.5, line = 2.5)
  title(xlab = xlab, cex.lab = 1.5, line = 2.5)
}

bool = function(x){
  u <- calc(x, fun = median, na.rm = TRUE)
  o <- calc(x, fun = mad, na.rm = TRUE)
  boo <- x > (u + o)
  extent(boo) <- extent(x)
  boo <- setZ(boo, z = getZ(x), name = "time")
  boo
}

vectorize  <- function(x) {
  sdate <- getZ(x)
  # x = raster::flip(x, direction = "x")
  # x = raster::flip(x, direction = "y")
  # 
  e = extent(x)
  x = t(x)
  extent(x) = e
  # 
  x     <- rasterToPoints(x)
  x     <- data.frame(x)
  colnames(x) <- c("lon", "lat", as.character(sdate))
  x     <- reshape2::melt(x, id.vars = c("lat", "lon"))
  colnames(x) <- c("lats", "lons", "time", "val")
  x
}

tsect <- function(x, y, z, xreach = 1, yreach = 1, xlen = 120, ylen = 40){
  #make a max min vector of sla and fsle
  coord = c(max(x, na.rm = TRUE), min(x, na.rm = TRUE), 
            max(y, na.rm = TRUE), min(y, na.rm = TRUE))
  
  # use maxmin vector to create fake x, y vectors (downsized) as meshgrid input
  y_vec <- seq(from = coord[4], to = coord[3], length.out = ylen) 
  # may need to invert to be same length 
  x_vec <- seq(from = coord[2], to = coord[1], length.out = xlen)
  xy_grid <- meshgrid(x_vec, y_vec)
  
  #inisilize a matrix of dim fs_grid[1] filled with NA values
  dims <- dim(xy_grid[[1]])
  z_grid <- matrix(data = NA, 
                   nrow = dims[1], 
                   ncol = dims[2], 
                   dimnames = NULL)
  
  for(iy in 1:dims[1]) {
    for(ix in 1:dims[2]){
      # where in the df is the difference greater than the units
      box_x <- which(abs(x - xy_grid[[1]][iy, ix]) <= xreach)
      box_y <- which(abs(y - xy_grid[[2]][iy, ix]) <= yreach)
      # I think the grids are dif sizes and should be subet differently
      #index vector of both cox_sla and box_fsle as one
      #box <- sort(match(box_y, box_x))
      
      # I think this is the correct way to do this
      box <- box_x[box_x %in% box_y]
      
      z_grid[iy, ix] <- mean(z[box], na.rm = TRUE)
    }
  }
  
  list(z_grid, y_vec, x_vec)
}

jamie_theme <- function(x,
                        y, 
                        line = 0.75,
                        adj = 0,
                        #asp = 1,
                        ylim = range(y),
                        xlim = range(x),
                        main = "",
                        labels = round(seq(range(x)[1]+1, range(x)[2]-1, 
                                           length = 6), 1),
                        at =  round(seq(range(x)[1]+1, range(x)[2]-1, 
                                        length = 6), 1),
                        dt = TRUE,
                        yaxes = TRUE,
                        xaxes = TRUE,
                        ylab = "",
                        xlab = "",
                        cex.main = 1.25
){
  
  if(dt == TRUE) labels <- format(as.Date(labels), "%b-%d")
  if(is.numeric(x)) x <- round(x, 1)
  if(is.numeric(y)) y <- round(y, 1)
  
  plot(0,
       ylim = ylim,
       xlim = xlim,
       main = "",
       xlab = "",
       ylab = "",
       axes = FALSE,
       #asp = asp,
  )
  grid(nx = NULL, # X-axis divided in two sections
       ny = NULL, # Y-axis divided in three sections
       lty = 2, 
       col = colvect(c("gray69"), alpha = 0.6), lwd = 1)
  # box(which = "plot", 
  #     lty = "solid", 
  #     lwd = 3, 
  #     col = colvect("grey22", alpha = 0.9))
  if(xaxes){axis(side = 1,
                 at = at,
                 labels = labels,
                 las = 1, 
                 lwd = 1.5, 
                 mgp = c(2, 1, 0), 
                 cex.axis = 1,
                 col = colvect("grey22", alpha = 0.9))}
  if(yaxes){axis(side = 2,
                 las  = 2, 
                 lwd  = 1.5, 
                 mgp  = c(1, 0.75, 0), 
                 cex.axis = 1,
                 col = colvect("grey22", alpha = 0.9))}
  title(main = main,
        cex.lab = cex.main,
        line = line,
        adj = adj)
  title(ylab = ylab, cex.lab = 1, line = 2.5)
  title(xlab = xlab, cex.lab = 1, line = 2.25)
}
# ------------------------------------------------------------------------------
# loading data via opendap had issues with over engineered function
sdate = as.Date("2018-01-01")
edate = as.Date("2019-01-01") 

lons = c(-160, -125)
lats = c( 18,   37)
e = extent(lons, lats)

# cache_list()
cache_delete_all(force = FALSE)

# Set variables
lat_varid = "lat"
lon_varid = "lon"
var = "CHL"
url = "https://jash:5.Pellegrino@my.cmems-du.eu/thredds/dodsC/cmems_obs-oc_glo_bgc-plankton_my_l3-multi-4km_P1D?"

origin = "1900-01-01"

data = nc_open(url, verbose = FALSE, write = FALSE)
lat  = ncvar_get(data, varid = lat_varid)
lon  = ncvar_get(data, varid = lon_varid)
time = ncvar_get(data, varid = "time")
time = as.Date(time, origin = origin)

idx_lat  = which(lat > lats[1] & lat < lats[2])
idx_lon  = which(lon > lons[1] & lon < lons[2])
idx_time = which(time >= sdate & time <= edate)

idx_ras = paste("CHL",
                paste("[", range(idx_time)[1], ":1:", range(idx_time)[2], "]", sep = ""),
                paste("[", range(idx_lat)[1],  ":1:", range(idx_lat)[2],  "]", sep = ""),
                paste("[", range(idx_lon)[1],  ":1:", range(idx_lon)[2],  "]", sep = ""),
                sep = "")

idx_time = paste("time", paste("[", range(idx_time)[1], ":1:",range(idx_time)[2], "]", sep = ""), sep = "")
idx_lat  = paste(lat_varid, paste("[", range(idx_lat)[1], ":1:",range(idx_lat)[2],  "]", sep = ""), sep = "")
idx_lon  = paste(lon_varid, paste("[", range(idx_lon)[1], ":1:",range(idx_lon)[2],  "]", sep = ""), sep = "")
idx = paste(idx_lat, idx_lon, idx_time, idx_ras, sep = ",")

url = paste(url, idx, sep = "")

nc_close(data)
rm(data)

data = nc_open(url, verbose = FALSE, write = FALSE)

lat  = ncvar_get(data, varid = lat_varid)
lon  = ncvar_get(data, varid = lon_varid)
time = ncvar_get(data, varid = "time")
time = as.Date(time, origin = origin)
ras  = ncvar_get(data)
nc_close(data)

s   = dim(ras)
ras = raster::brick(ras)
ras = t(ras)
ras = setZ(ras, z = as.Date(time, origin = org), name = "time")
extent(ras) = extent(min(lons), max(lons), min(lats), max(lats))

chl = ras
rm(ras)
gc() 

# ------------------------------------------------------------------------------
# Loading the float data using ERDDAP

fields =  c("temp",
            "temp_qc",
            "temp_adjusted", 
            "psal", 
            "psal_adjusted", 
            "pres", 
            "pres_adjusted", 
            "time",
            "cycle_number", 
            "latitude", 
            "longitude",
            "platform_number")

float = loadfloat(sdate  = sdate, 
                  edate = edate, 
                  lons  = lons, 
                  lats  = lats,
                  fields = fields,
                  id = "ArgoFloats")
gc()

# ------------------------------------------------------------------------------
# removing coastal effects
chl = bufcoast(chl, 
               region = "Hawaiian Islands", 
               path = "data/USMaritimeLimitsAndBoundariesSHP")
gc()

# ------------------------------------------------------------------------------
# Calculate CHL anomaly
chla = anomalize(chl, detrend = TRUE)
chlb = bool(chla)
gc()

# ------------------------------------------------------------------------------
# checking how well the filter did
# c_time = cellStats(chla, stat='mean', na.rm=TRUE)
# plot(1:length(c_time), c_time)

# ------------------------------------------------------------------------------
# Subset time-span of interest: 2018 summer bloom
# I could also just load using errdap. This is faster tho
#ppp_2018  = timesnip(ppp,  as.Date(sdate), as.Date(edate)
idx = year(blooms$sdate) == 2018
sdate = as.Date(blooms$sdate[idx]) - 16
edate = as.Date(blooms$edate[idx]) + 15

chl_2018  = timesnip(chl,  as.Date(sdate), as.Date(edate))
chla_2018 = timesnip(chla, as.Date(sdate), as.Date(edate))
chlb_2018 = timesnip(chlb, as.Date(sdate), as.Date(edate))

# ------------------------------------------------------------------------------
# Calculate patch area
# integrated = depth integrated
tsa = bloom_area(chlb_2018)
# magnitude might like the raw CHL
tsc = bloom_con(chl_2018, chlb_2018, stat = "mean")

#tsa = area_ts(chla, sdate = sdate, edate = edate, stat = "mean")
bloom_2018 = data.frame(area = tsa$area, 
                        con  = tsc$con, 
                        time = tsa$time)

# cleaning work space
rm(tsa, tsc)

# ------------------------------------------------------------------------------
# Process float data
# time-series of the mixed layer depth for the 
# Converting all columns to numeric
float$time = as.numeric(as.Date(float$time))

#float$float_serial_no = as.factor(float$float_serial_no)
float$float_serial_no = as.factor(float$platform_number)

float = data.frame(lapply(float, function(x) as.numeric(as.character(x))))

# removing bad float data
float = subset(float, 
               psal > 3 
               & temp > 5 
               & temp < 45 
               & pres > 0
               & pres < 800)

#float$float_serial_no <- as.numeric(float$float_serial_no)
float$float_serial_no <- as.numeric(float$platform_number)
float$time            <- as.numeric(float$time)

# Identifying bad floats and removing them
badf  = subset(float, temp < 10 & pres < 100)
#badid = is.element(float$float_serial_no, badf$float_serial_no)
badid = is.element(float$float_serial_no, badf$platform_number)
badcy = is.element(float$cycle_number, badf$cycle_number)
float = subset(float, !badid & !badcy)
rm(badcy, badid, badf)

# Identifying bad floats and removing them
badf = subset(float, temp > 20 & pres > 300)
#badid = is.element(float$float_serial_no, badf$float_serial_no)
badid = is.element(float$float_serial_no, badf$platform_number)
badcy = is.element(float$cycle_number, badf$cycle_number)
float = subset(float, !badid & !badcy)

rm(badcy, badid, badf)

# ------------------------------------------------------------------------------
# Calculate density and mixed layer depth
# I changed the region so floats suck

# Density Calculation: function from the oce package and calc. density
float$rho = swRho(salinity = float$psal,
                  temperature = float$temp,
                  pressure = float$pres,
                  latitude = float$latitude,
                  eos = getOption("oceEOS", default = "gsw"))
# I should make sure float factors are not lost to NA values

#float = subset(float, !is.na(float_serial_no))
float = subset(float, !is.na(platform_number))

# Mixed layer depth per profile
#id = unique(float[, c("float_serial_no", "cycle_number")])
id = unique(float[, c("platform_number", "cycle_number")])

mixed = data.frame(matrix(data = NA, 
                          nrow = nrow(id), 
                          ncol = ncol(float)+1, 
                          dimnames = NULL))
colnames(mixed) = colnames(float)

for(i in 1:nrow(id)){
  oneprof = subset(float, 
                   float_serial_no == id$platform_number[i] & 
                     cycle_number == id$cycle_number[i])
  if(nrow(oneprof) > 4) mixed[i, ] = mldchu(ctd = oneprof, 
                                            x = "rho", 
                                            y = "pres_adjusted", 
                                            n = 10)
}

# ------------------------------------------------------------------------------
# Last minute remove bad data
mix = mixed
mix = subset(mixed, !pres == 8)
mix = subset(mix, !pres > 100)
chl_float = subset(float, pres < 40) 

# ------------------------------------------------------------------------------
# Loess running average
# I may want to change tsa to include magnitude. Maybe keep separate
t     = as.numeric(bloom_2018$time)

curve = loess(mix$pres_adjusted ~ mix$time) 
x_mod = seq(min(t), max(t), length.out = length(t))
mix_mod = predict(curve, newdata =  x_mod, type = "response")

bloom_2018$mld = mix_mod

# ------------------------------------------------------------------------------
# CHL concentration and units conversion
# handling clouds. I hope hawiann islands are not in this extent
s = dim(chl_2018)
n = s[1] * s[2] * s[3]
# total cloud cover percent across all data
clouds = sum(is.na(values(chl_2018))) / n
clouds = unname(clouds)
# True * clouds = measured and moving clouds over. Should be less than 1. 
clouds = 1 - clouds
# adding the missing percent back assuming even cloud cover
bloom_2018$area =  bloom_2018$area/clouds

# ------------------------------------------------------------------------------

# In mg
bloom_2018$bio = (bloom_2018$mld * bloom_2018$con * bloom_2018$area * 1e+06)

# mg to kilotons 
bloom_2018$bio = bloom_2018$bio/(1e+12)
bloom_2018$area = bloom_2018$area/1e+05 #10e05
bloom_2018$con = bloom_2018$con # in mg/m^3
elephants = max(bloom_2018$bio*1000)/ 5.98742 # elephant in t

# ------------------------------------------------------------------------------
# saving the data to harddisck
# write.csv(bloom_2018, 
#           file = paste("data/fig_data/lehan", Sys.Date(), ".csv", sep = ""))

bloom_2018 = read.csv("data/fig_data/lehan2023-05-06.csv")
bloom_2018$time = as.Date(bloom_2018$time)
# ------------------------------------------------------------------------------

# ------------------------------- FIGURES --------------------------------------

# ------------------------------------------------------------------------------

# ------------------------------------------------------------------------------
# subsetiign for plotting
tmp = subset(bloom_2018, time < as.Date("2018-10-20"))
tmp = subset(tmp, time > as.Date("2018-07-15"))
head(tmp)
# ------------------------------------------------------------------------------
# Biomass Timesries 
dt = gsub("-", "", as.character(Sys.Date()))
pdf(paste("figures/biomass_timeseries_", dt, ".pdf", sep = ""),   # The directory you want to save the file in
    width = 5, # The width of the plot in inches
    height = 6, # The height of the plot in inches
    pointsize = 15) 

# saving to an eps file
setEPS()
postscript(paste("figures/biomass_timeseries_", Sys.Date(), ".eps", sep = ""))

layout(matrix(c(1,  
                2,
                3), 
              nrow = 3, 
              ncol = 1, 
              byrow = TRUE),
       widths  = c(1),
       heights = c(1, 1, 1.15))

# Area
par(mar = c(0, 5, 2, 3))
x = as.numeric(tmp$time)
y = tmp$area
jamie_theme(x = x, 
            y = y, 
            dt  =FALSE,
            ylim = c(0, 35),
            main = "", 
            xlab = "",
            ylab = expression(Surface ~ Area ~ (km^2 ~ x10^{5})),
            line = 0.75,
            adj = 0,
            yaxes = TRUE,
            xaxes = FALSE)
points(x, y,  pch = 20, cex = 1) 

# CHL Concentration
par(mar = c(1, 5, 1, 3))
x = tmp$time
y = tmp$bio
jamie_theme(x = x,
            y = y,
            ylim = c(0, 12),
            dt = FALSE,
            line = 0.75,
            adj = 0.0,
            main = "",
            xlab = "",
            yaxes = TRUE,
            ylab = expression(Chlorophyll ~ Biomass ~ (kt)),
            xaxes = FALSE)
points(x, y, pch = 20, cex = 1)

# CHL Summed Magnitude
par(mar = c(4, 5, 0, 3))
y = tmp$con/1000
x = as.Date(tmp$time)
jamie_theme(x = x, 
            y = y, 
            ylim = c(0.06, 0.15),
            dt = TRUE,
            line = -1,
            adj = 0.03,
            main = "", 
            xaxes = TRUE,
            yaxes = TRUE,
            xlab = "Month Day",
            ylab = expression(Chlorophyll ~ (mg ~ m^{-3}))
            )
points(x, y,  pch = 20, cex = 1)

dev.off()

# ------------------------------------------------------------------------------
# png(filename= paste("mixed_layer_", Sys.Date(), ".png", sep = ""),
#     width = 3.5,
#     height = 6,
#     units = "in",
#     res = 300,
#     pointsize = 10)
# 
# # Mixed Layer Depth
# curve = loess(mix$pres_adjusted ~ as.numeric(mix$time)) 
# curve = predict(curve, type = "response")
# 
# jamie_theme(x = mix$time, y = mix$pres_adjusted, 
#             ylim = c(80, 0),
#             main = "2018 Bloom Mixed Layer Depth", 
#             xlab = "Date Time",
#             ylab = "Pressure [dbar]")
# points(mix$time, mix$pres_adjusted, pch = 20, cex = 0.6)
# points(mix$time, curve, pch = 20, cex = 1)
# points(bloom_2018$time, bloom_2018$mld, col = "red", lwd = 2, pch = 20)
# dev.off()
# ------------------------------------------------------------------------------
# Float Tracks

e = extent(-165, -132, 18, 37)
wdmap = getMap(resolution = "high")

#id <- unique(float[, c("longitude", "latitude", "float_serial_no", "time")])
id <- unique(float[, c("longitude", "latitude", "platform_number", "time")])
#id$float_serial_no <- as.factor(id$float_serial_no)
id$float_serial_no <- as.factor(id$platform_number)
#nb.cols <- length(unique(id$float_serial_no))
nb.cols <- length(unique(id$platform_number))
mycolors <- colorRampPalette(brewer.pal(12, "Paired"))(nb.cols)
#levels(id$float_serial_no) <- mycolors
levels(id$platform_number) <- mycolors

png(filename= paste("float_tracks_", Sys.Date(), ".png", sep = ""),
    width = 3.5,
    height = 6,
    units = "in",
    res = 300,
    pointsize = 10)

jamie_theme(x = id$longitude, y = id$latitude,
            dt = FALSE,
            main = "2018 bloom argo-floats tracks",
            xlab = "Longitude",
            ylab = "Latitude",
            xlim = e[1:2],
            ylim = e[3:4])
sz <- scale01(id$time)
points(id$longitude, id$latitude,
       #col = as.character(id$float_serial_no),
       col = as.character(id$platform_number),
       cex = sz*1.5,
       pch = 20)
plot(wdmap,
     xlim = e[1:2],
     ylim = e[3:4],
     asp = 1,
     bg = "black",
     border = "black",
     col = "black",
     #wrap = c(-180, 180),
     add = TRUE)
box(which = "plot", lty = "solid", lwd = 3, col = colvect("grey22", alpha = 0.9))

dev.off()


