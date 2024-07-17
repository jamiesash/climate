library("rerddap")
library("ncdf4")
library("terra")

# ------------------------------------------------------------------------------
# Reading entire data set at once. 
# This will only work for data sets less than 2Gb. 

sstainfo = info('jplMURSST41anom1day')
murSSTA = griddap(sstainfo, 
                  stride = c(1, 50, 50),
                  latitude = c(16, 40), 
                  longitude = c(-175, -130), 
                  time = c("2018-06-01", "2018-11-01"), 
                  fields = 'sstAnom')

data = nc_open(murSSTA$summary$filename)
ras  = ncvar_get(data, varid = "sstAnom")
lats = ncvar_get(data, varid = "latitude")
lons = ncvar_get(data, varid = "longitude")
time = ncvar_get(data, varid = "time")
time = as.Date(as.POSIXct(time, origin = "1970-01-01"))
nc_close(data)

e = ext(min(lons), max(lons), min(lats), max(lats))
sst = rast(ras, extent = e)
time(sst) = time

writeCDF(sst, 
         filename = paste("/home/jamesash/projects/climate/data/", "ssta_l4_lowres_", dt, ".nc",sep = ""), 
         varname = "ssta",
         overwrite = TRUE)

# ------------------------------------------------------------------------------
# Reading and writing a large data set in chunks. 
dap = function(sdate, edate, x, y, 
               stride = 1,
               id = 'jplMURSST41anom1day', 
               var = c("sstAnom", "latitude", "longitude", "time")) {
    dap_that_ass = function(day, data_info, x, y){
      data = griddap(data_info, 
                     latitude = y[1:2],
                     stride =  c(1, 5, 5),
                     longitude = x[1:2], 
                     time = c(as.character(day), as.character(day)),
                     fmt = "nc",
                     fields = 'all')
      data = nc_open(data$summary$filename)
      ras  = ncvar_get(data, varid = var[1])
      lats = ncvar_get(data, varid = var[2])
      lons = ncvar_get(data, varid = var[3])
      t = ncvar_get(data, varid = var[4])
      t = as.Date(as.POSIXct(t, origin = "1970-01-01"))
      nc_close(data)
      
      # set all the variables. 
      e = ext(min(lons), max(lons), min(lats), max(lats))
      ras = rast(ras, extent = e)
      time(ras) = t
      ras
    }
    
    sdate = as.Date(sdate)
    edate = as.Date(edate)
    t = seq(as.Date(sdate), as.Date(edate), by = 1)
    t = as.list(t)
    data_info = info(id)
    
    # download all rasters as a list.
    ras  = lapply(t, 
                  FUN = dap_that_ass, 
                  data_info = data_info, 
                  x = x, 
                  y = y)
    
    ras  = rast(ras)
    # fix it right. 
    e = ext(ras)
    ras = t(ras)
    ras = flip(ras)
    ext(ras) = e
    ras
}

# Note this fills the hell out of the tmp directory /tmp/RtmpjCHdUF/R/rerddap/
# /tmp/RtmpTTr6r3/R/rerdda
# /tmp/RtmpKdtXge/R/rerddap
# Will fail if there are corrupted files in that directory. 
sst = dap(sdate = "2018-06-01",
          edate = "2018-11-01",
          x = c(-175, -130), 
          y = c(16, 40))

dt = gsub("-", "", as.character(Sys.Date()))
writeCDF(sst, 
         filename = paste("/home/jamie/projects/climate/data/", "ssta_l4_2018_lowres_", dt, ".nc", sep = ""), 
         overwrite = TRUE)





