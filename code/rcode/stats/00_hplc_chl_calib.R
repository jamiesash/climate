# Climotology ------------------------------------------------------------------
# Loading packages and all prewriten functions
# setwd("C:\\Users\\james\\Desktop\\jamieslife\\analysis")
source("code\\libraries.R")
source("code\\functions.R")
library(MASS)
library(ggeffects)
rasterOptions(maxmemory = 100e+10)

lons = c(-158.25, -157.25)
lats = c( 22.7,   22.8)
e = extent(lons, lats)

# ------------------------------------------------------------------------------
# for the box plot
boxit <- function(x1, x2, y, data,
                  xaxes = TRUE, 
                  yaxes = TRUE,
                  col=c("slateblue1" , "tomato"),
                  ylim = c(min(y, na.rm = TRUE), max(y, na.rm = TRUE)),
                  sub = "Car Milage Data",
                  xlab = "", 
                  main = "",
                  ylab = "Miles Per Gallon",
                  legend = c("ALOHA", "30N"),
                  add = FALSE,
                  at = 1:24,
                  labels = as.character(1:12),
                  cex.main = 1.5,
                  cex.lab = 1.5,
                  mar = c(5, 5, 5, 3)) {
  par(mar = mar)
  boxplot(y ~ x1*x2, 
          data = box_df,
          ylim = ylim,
          main = main,
          col=col,
          xlab = "", 
          ylab = "", 
          add = add,
          whisklty = 1,
          medlty = "blank",
          outline = FALSE, 
          axes = FALSE)
  legend(x = 9, 
         y = 0.21, 
         legend = legend, 
         pch = 22, 
         bty ="n",
         pt.bg = col, 
         cex = 1)
  if(yaxes) axis(side = 2,
                 las = 2, 
                 lwd = 1, 
                 mgp = c(1, 0.75, 0), 
                 cex.axis = 1)
  if(xaxes) axis(side = 1, 
                 las = 1, 
                 lwd = 1, 
                 at = at,
                 labels = labels, 
                 mgp = c(2, 1, 0), 
                 cex.axis = 1)
  box(which = "plot", lty = "solid", lwd = 1.5, col = "black")
  title(ylab = ylab, 
        cex.lab = 1,
        line = 2.5)
  title(xlab = xlab, 
        cex.lab = 1,
        line = 2.5)
}

# for the scatter plot
plotit <- function(x, y, 
                   asp = 1,
                   abline = c(m, b, r),
                   tex_x = 0.10,
                   tex_y = 0.05,
                   xaxes = TRUE, 
                   yaxes = TRUE,
                   ylim = c(min(y, na.rm = TRUE), 
                            max(y, na.rm = TRUE)),
                   xlim = c(min(x, na.rm = TRUE), 
                            max(x, na.rm = TRUE)),
                   sub = "Car Milage Data",
                   xlab = "", 
                   main = "",
                   pch = 19,
                   ylab = "Miles Per Gallon",
                   cex.main = 1.5,
                   cex.lab = 1.5,
                   mar = c(5, 5, 5, 3)) {
  
  # This is for a linear regression
  xm <- seq(xlim[1]-xlim[1]*2, xlim[2]+xlim[2], length = 100)
  ym <- xm * abline[1] + abline[2]
  
  par(mar = mar)
  #plot.new()
  plot(1,
       asp = asp,
       ylim = ylim,
       xlim = xlim,
       main = main,
       xlab = "",
       ylab = "",
       axes = FALSE)
  #text(0.19, 0.02, pos = 1,
  #     labels = paste('r =', as.character(abline[3]), sep = " "))
  text(tex_x, tex_y, pos = 1,
       labels = paste('y =', 
                      as.character(abline[1]), "* x",
                      "+",
                      as.character(abline[2]), 
                      sep = " "))
  grid(nx = NULL, # X-axis divided in two sections
       ny = NULL, # Y-axis divided in three sections
       lty = 2, 
       col = colvect(c("gray69"), alpha = 0.6), 
       lwd = 1)
  lines(xm, ym, lwd = 1,  col = "grey12")
  points(x, y, 
         pch = pch, 
         col = colvect(c("gray12"), alpha = 0.8))
  if(yaxes) axis(side = 2,
                 las = 2, 
                 lwd = 1, 
                 mgp = c(1, 0.75, 0), 
                 cex.axis = 1)
  if(xaxes) axis(side = 1, 
                 las = 1, 
                 lwd = 1, 
                 mgp = c(2, 1, 0), 
                 cex.axis = 1)
  box(which = "plot", lty = "solid", lwd = 1.5, col = "black")
  title(ylab = ylab, 
        cex.lab = 1,
        line = 2.5)
  title(xlab = xlab, 
        cex.lab = 1,
        line = 2.5)
}

# ------------------------------------------------------------------------------

lons = c(-158, -157)
#lons = c(-158.05, -157.95)
lats = c( 21.5,   22.5)
e = extent(lons, lats)

# Set variables
lat_varid = "lat"
lon_varid = "lon"
sdate = as.Date("1997-06-01")
edate = as.Date("2022-12-31")
var = "CHL"
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

lons = c(-158.1, -157.95)
lats = c( 22.65,   22.8)
e = extent(lons, lats)
chl = raster::crop(chl, e)
dim(chl)

# ------------------------------------------------------------------------------
# load hplc data
# Load the HPLC Data text file from HOT_DOGS webserver
hplc <- data.frame(fread("data//hotdogs//HPLC_HOT_19881001_20211203.txt", 
                         colClasses = c(rep("character", 6)), 
                         header=TRUE),
                   stringsAsFactors = FALSE)
colnames(hplc) <- c("id", "date", "time", "press", "chl", "NA")

# ------------------------------------------------------------------------------
# Could try apply a gassion running smoother like I did for the magnitude
# satelite data is gotten noisy in the last few years

c_t = getZ(chl)
c_u = cellStats(chl, stat = mean, na.rm = TRUE)

sat = data.frame(c_t, c_u)
colnames(sat) = c("time", "chlu")
sat$month = months(sat$time)
sat$month = factor(sat$month, levels = month.name)

sat = data.frame(chl = sat$chlu, month = sat$month, time = sat$time)
sat$ID = "satelite"
sat$ID = as.factor(sat$ID)

# ------------------------------------------------------------------------------
# processing HPLC data
# I should still probably average HPLC data from the same day across depth

# convert the odd datetype to Postix
time   = as.POSIXct(hplc$date, format="%m%d%y")
hplc <- data.frame(apply(hplc[,c("press", "chl")], 
                         MARGIN=2, 
                         FUN=as.numeric))

hplc      <- cbind(hplc, time) # create a dataframe
hplc$chl <- hplc$chl/1000 # convert ng to mg
hplc$time <- as.Date(hplc$time) # that Postix is now a date class
hplc        <- hplc[which(hplc$press < 27), ] # We only want the upper 25m

# Remove old data dates
ind   <- which(hplc$time >   as.Date(min(sat$time)))
#as.Date("1997-06-01")) #as.Date(min(sat$time)))
hplc  <- hplc[ind,] 

# botlle samples used for paper reference
n = length(unique(hplc$time))

# find average of all depth samples 5m and 25m
hplc = aggregate(chl ~ time, data = hplc, FUN = mean, na.rm = TRUE)

# months for idk
hplc$month <- months(hplc$time)
# ------------------------------------------------------------------------------
# be carefull here the months are not ordered correctly
hplc$month = factor(hplc$month, levels = month.name)
hplc$ID = as.factor("hplc")

hplc = hplc[, c("chl", "month", "ID", "time")]
sat = sat[, c("chl", "month", "ID", "time")]

# ------------------------------------------------------------------------------
# Fixin date times
# both need to be true and I'm gutch
max(hplc$time) < max(sat$time)
min(hplc$time) > min(sat$time)

# this can't be right
box_df = rbind(sat, hplc)
box_df = box_df[order(box_df$month),]

# indexing closest value very manually because something is up
j = 0
idx = vector()
for(i in 1:nrow(hplc)){
  idx[i] = which(abs(sat$time - hplc$time[i]) == min(abs(sat$time - hplc$time[i])))
}

j = 0
cu = vector()
for(i in idx){
  j = j+1
  cu[j] = mean(sat$chl[(i-1):(i+1)], na.rm = TRUE)
}

# there should not be any duplicated values
duplicated(idx)

calib = data.frame(time = hplc$time, sat_chl = cu, hplc_chl = hplc$chl)

# ------------------------------------------------------------------------------

# GLM ussing Gamma Log transform x 

x <- calib$hplc_chl
y <- calib$sat_chl
data <- data.frame(y, x)  #dataset

glmGamma <- glm(y ~ x, data = data, family = Gamma(link = "identity"))
glm_0 <- glm(y ~ 1, data = data, family = Gamma(link = "identity"))

model1 <- summary(glmGamma)
anova(glmGamma, glm_0, test = "F")

m  <- round(model1$coefficients[2], 3)
b  <- round(model1$coefficients[1], 3)
r2 <- round(1 - model1$deviance/model1$null.deviance, 3)

# ------------------------------------------------------------------------------

dt = gsub("-", "", as.character(Sys.Date()))
pdf(paste("figures\\hplc_chlsat_", dt, ".pdf", sep = ""),   # The directory you want to save the file in
    width = 5, # The width of the plot in inches
    height = 6,
    pointsize = 10) # The height of the plot in inches

# setEPS()
# postscript(paste("figures\\hplc_chlsat_", dt, ".eps", sep = ""))

# Boxplot
layout( matrix(c(1,
                 2),
               ncol = 1,
               nrow = 2, 
               byrow = TRUE),
        heights = c(1.85, 1))

# i do not know
n = length(y) - sum(is.na(y))
plotit(x = x, y = y, abline = c(m, b, r2),
       tex_x = 0.13,
       tex_y = 0.2,
       pch = 4,
       xlab = expression(HPLC ~ Bottle ~ Samples ~ (mg ~ L^{-1})),
       ylab = expression(GlobColor ~ Daily ~ GSM ~ CHL1 ~ L3 ~ (mg ~ m^{-3})),
       sub = "",
       ylim = c(0.02, 0.25),
       xlim = c(0.02, 0.25),
       mar = c(4, 4, 2, 2))
#title("HPLC to CHLsat Comparison and CHL Seasonality", line = 1, adj = 0, cex.main = 1.5)
text(0.13, 0.17, labels = paste("n =", as.character(n), sep = " "))

# plotting the box by itself
boxit(x2 = box_df$month,
      x1 = box_df$ID,
      col = colvect(c("ivory3", "grey12"), alpha = c(0.5,0.9)),
      add = FALSE,
      y = box_df$chl,
      ylim = c(0.025, 0.22),
      legend = c("HPLC", "Satelite"),
      data = box_df,
      labels = as.character(1:12),
      at = seq(from = 1, to = 24, by = 2) + 0.5,
      sub = "",
      xlab = expression(Month ~ of ~ year),
      ylab = expression(Mean ~ CHL ~ (mg ~ m^{-3})),
      mar = c(4, 4, 0, 2))

dev.off()

# ------------------------------------------------------------------------------
dt = gsub("-", "", as.character(Sys.Date()))
pdf(paste("figures\\hplc_chlsat_", dt, ".pdf", sep = ""),   # The directory you want to save the file in
    width = 6, # The width of the plot in inches
    height = 4,
    pointsize = 10) # The height of the plot in inches

# plotting the box by itself
boxit(x2 = box_df$month,
      x1 = box_df$ID,
      col = colvect(c("ivory4", "grey22"), alpha = c(0.5,0.9)),
      add = FALSE,
      y = box_df$chl,
      ylim = c(0.025, 0.22),
      legend = c("HPLC", "Satelite"),
      data = box_df,
      labels = as.character(1:12),
      at = seq(from = 1, to = 24, by = 2) + 0.5,
      sub = "b)",
      xlab = "Month of year",
      ylab = "Mean CHL [mg/m^3]",
      line = -1.1,
      adj = 0.025,
      mar = c(5, 5, 3, 3))



           