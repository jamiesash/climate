library(mgcViz)
library(lubridate)
library(terra)

# ------------------------------------------------------------------------------
### Script To Do

# 1. Resize to CHL for higher fake resolution. Use nearest neighbor. 
# 2. Calculate distance not lat/lon for the GAM. 
# 3. Download a different CHL this one is cloudy. 
# 4. Mask does not remove the center of the island. 
# 5. Re-download SST I deleted it by accident. 
# 6. SLA resize introduces extreme negative number.  

# ------------------------------------------------------------------------------
# Load SLA, CHL and SSTA data sets. 
# With this method I prefer to download the extent and time span of interest. 

setwd("/home/jamie/projects/climate/code/rcode/casestudy")

sla = rast("../../../data/sla/sla_2018_l4_4k.nc")
chl = rast("../../../data/chl/chl_2018_daily_multi_l3_4k.nc")
sst = rast("../../../data/sst/ssta_l4_2018_4k_20240717.nc")
chl = rast("../../../data/chl/chl_2018_glob_daily_multi_l3_4k.nc")

# subset spacial polygon like a data frame. 
eez = terra::vect("../../../data/eez/USMaritimeLimitsNBoundaries.shp")
idx = which(eez$REGION == 'Hawaiian Islands')
hawaii = eez[idx,]

# ------------------------------------------------------------------------------
# removing coastal effects

idx = hawaii$CZ + hawaii$TS
idx = which(idx == 1)
hawaii = hawaii[idx,]

# ------------------------------------------------------------------------------
# Resize to SLA since it is the lowest resolution. 

sst = resample(sst, sla, method = "near")
chl = resample(chl, sla, method = "near")
# sla = resample(sla, sla, method = "near")

# So close. Close enough for now. Still leaves little blips. 
sst = terra::mask(sst, hawaii, inverse = TRUE)
chl = terra::mask(chl, hawaii, inverse = TRUE)
sla = terra::mask(sla, hawaii, inverse = TRUE)

# ------------------------------------------------------------------------------
# Create table

# Welp this used to be really hard before terra. 
chl = as.data.frame(chl, xy = TRUE, wide = FALSE, na.rm=FALSE, time=TRUE)
sla = as.data.frame(sla, xy = TRUE, wide = FALSE, na.rm=FALSE, time=TRUE)
sst = as.data.frame(sst, xy = TRUE, wide = FALSE, na.rm=FALSE, time=TRUE)

# ------------------------------------------------------------------------------
# remove clouds. 

# This seems not to work. 
idx = which(!is.na(chl$values))
chl = chl[idx,]
sla = sla[idx,]
sst = sst[idx,]

time = sla$time
time = yday(time)
alldata = data.frame(time = time, 
                     lon = sla$x, 
                     lat = sla$y, 
                     sla = sla$values, 
                     chl = chl$values,
                     sst = sst$values)
# KOA might not like this. 
rm(sla, chl, sst, time)

# ------------------------------------------------------------------------------

# Subset data so I can easily run this on my local computer. 
idx = 1:nrow(alldata)
idx = sample(idx, 100000)
smalldata = alldata[idx,]

# ------------------------------------------------------------------------------
# going for it. 

gam_chl_sla = gam(chl ~ sla + s(time) + s(lon, lat), 
              data = smalldata,
              # method = "REML",
              family = Gamma(link = "inverse"))

gam_sla_sst = gam(sst ~ sla + s(time) + s(lon, lat), 
              # method = "REML",
              data = smalldata)

# save the model output. 
save(gam_chl_sla, file="../../../data/gam/chl_gam.Rdata")
save(gam_sla_sst, file="../../../data/gam/sla_gam.Rdata")

# ------------------------------------------------------------------------------
# plotting the GAM output. 

gam_chl_sla_visual <- getViz(gam_chl_sla)
gam_sla_sst_visual <- getViz(gam_sla_sst)

dt = gsub("-", "", as.character(Sys.Date()))
pdf(paste("../../../figures/chl_sla_gam_", dt, ".pdf", sep = ""),   # The directory you want to save the file in
    width = 8, # The width of the plot in inches
    height = 8,
    pointsize = 10) # The height of the plot in inches
gridPrint(plot(pterm(gam_chl_sla_visual, 1)) + xlab("SLA [m]") + ylab(expression(CHL ~ Effect ~ (mg ~ m^{-3}))) + l_ciPoly() + l_fitLine(),
          ncol = 1)
dev.off()

dt = gsub("-", "", as.character(Sys.Date()))
pdf(paste("../../../figures/sla_sst_gam_", dt, ".pdf", sep = ""),   # The directory you want to save the file in
    width = 8, # The width of the plot in inches
    height = 8,
    pointsize = 10) # The height of the plot in inches
gridPrint(plot(pterm(gam_sla_sst_visual, 1)) + xlab("SSTA [C]") + ylab(expression(SLA ~ Effect ~ (m))) + l_ciPoly() + l_fitLine(),
          ncol = 1)
dev.off()

dt = gsub("-", "", as.character(Sys.Date()))
pdf(paste("../../../figures/gam_subplots_", dt, ".pdf", sep = ""),   # The directory you want to save the file in
    width = 9, # The width of the plot in inches
    height = 4,
    pointsize = 10) # The height of the plot in inches
gridPrint(plot(pterm(gam_sla_sst_visual, 1)) + xlab(expression(SLA ~ (m))) + ylab(expression(SSTA ~ Effect ~ (C))) + l_ciPoly() + l_fitLine(),
          plot(pterm(gam_chl_sla_visual, 1)) + xlab(expression(SLA ~ (m))) + ylab(expression(CHL ~ Effect ~ (mg ~ m^{-3}))) + l_ciPoly() + l_fitLine(),
          ncol = 2)
dev.off()

# ------------------------------------------------------------------------------

# Preparing a data frame for the GAM
# remove effects of clouds
# xyz = xyz[!is.na(xyz$sla),]
# xyz$time = as.Date(xyz$time)
# gc()

# # Converting lat lon to distance
# df = xyz
# xy   = cbind(df$lat, df$lon)
# utms = rgdal::project(xy, "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0") # need to change zones
# df$northing = utms[, 1]/1000
# df$easting  = utms[, 2]/1000
# df$doy = as.numeric(format(df$time, '%j'))
# 
# df = df[, c("chl", "sta", "sla", "doy", "easting", "northing")]
# 
# # removing outliers. May not be necessary. 
# idx = df$slaa < median(df$slaa) + mad(df$slaa)*5
# df = df[idx, ]
# idx = df$sla > median(df$sla) - mad(df$slaa)*5
# df = df[idx, ]
# gc()

# # ------------------------------------------------------------------------------
# # perform GAM 
# gam_chl = gam(chl ~ s(sta) ~ s(sla) + s(doy) + s(easting, northing), 
#               #method = "REML",
#               data = df, 
#               family = Gamma(link = "inverse"))
# 
# gam_chl = gam(chl ~ s(sla) + s(doy) + s(easting, northing), 
#               #method = "REML",
#               data = df, 
#               family = Gamma(link = "inverse"))
# 
# gc()
# 
# gam_ssta = gam(sta ~ sla + s(doy) + s(easting, northing), 
#                method = "REML",
#                data = df, 
#                family = gaussian(link = "identity"))
# 
# gc()

# gam_2018 = gam(chl ~ s(doy, by = regions, k=4) + s(doy) #+ s(easting, northing), 
#                method = "REML",
#                data = df, 
#                family = gaussian(link = "identity"))
# ------------------------------------------------------------------------------
