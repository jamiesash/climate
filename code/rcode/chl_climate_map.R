# ----------------------------------------------------------------------------
# The purpose of this script is to perform a linear regression on the chl 
# dataset from hopefull all the data 1997-2023, then to plot just the slope 
# of each regression on a map of the study region. 
# ----------------------------------------------------------------------------
# functions
# A function to load the dataset from opendap. 
# This function is copy paste from teh functions script. 
opendap = function(url = "https://jash:5.Pellegrino@my.cmems-du.eu/thredds/dodsC/cmems_obs-oc_glo_bgc-plankton_my_l4-gapfree-multi-4km_P1D?",
                   lons = c(-165, -135),
                   lats = c( 17,   35),
                   sdate = as.Date("2018-07-01"),
                   edate = as.Date("2018-11-30"),
                   lat_varid = "lat",
                   lon_varid = "lon",
                   var = "CHL",
                   origin = "1900-01-01"){
  # delete all built up files in cache
  # cache_list()

  cache_delete_all(force = FALSE)
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
  ras
}

# ------------------------------------------------------------------------------
# employing the function to download the data for the region of interest. 
chl = opendap(lats = c(17,45), 
	      lons = c(-175, -125), 
              var = "CHL", 
              sdate = "1998-01-01", 
	      edate = "2023-12-30")

# ------------------------------------------------------------------------------
# calculate CHL anomaly
# i receive an warning when detrending for cells of all NA values (hawaii)
# chl anom. calc. needs at least one full year to work correctly
chla = anomalize(chl, detrend = FALSE)
rm(chl)
gc()

# -------------------------------------------------------------------------------
# perform the regression. either as a for loop across each indecies or using a 
# prewritten function. 

# -------------------------------------------------------------------------------

