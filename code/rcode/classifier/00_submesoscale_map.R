
source("code//functions.R")
source("code//libraries.R")
library("rerddap")
library("akima")
library("dplyr")
library("ggplot2")
library("mapdata")
library("ncdf4")
library("plot3D")
library("IndexNumR")
library("imputeTS")
library("sp")
library("raster")
library("plotrix") 
rasterOptions(maxmemory = 120e+10, memfrac = 0.9)

# load pre-constructed data set
xyz = read.csv("data\\collated\\xyz_k_20230131.csv")

# ------------------------------------------------------------------------------
# Replace k-means better by day produce table
veron = function(x, y){
  y = scale01(y)
  x = scale01(x)
  
  us = median(x, na.rm = TRUE)
  uf = median(y, na.rm = TRUE)
  ms = mad(x, na.rm = TRUE)
  mf = mad(y, na.rm = TRUE)
  
  cents = rbind(c(us,      uf),
                c(us - ms*3, uf),
                c(us + ms*3, uf),
                c(us, uf + mf))
  
  cents = data.frame(cents)
  colnames(cents) = c("x", "y")
  
  l = list()
  dist = vector()
  for(i in 1:nrow(cents)){
    dist = sqrt((x - cents[i,1]) ^ 2 + (y - cents[i,2]) ^ 2)
    l[[i]] = dist
  }
  dist = do.call(cbind, l)
  k = apply(dist, MARGIN = 1, FUN = which.min)
  k = as.numeric(k)
  k
}

k = veron(x = xyz$slaa, y = xyz$fsle)

xyz$k = k

gc()

# ------------------------------------------------------------------------------
# plotting eddy regions. If I want to be extra I could overlay eddies 
# Make a color scale for each sub-mesoscale  region

df = xyz

temp = split(df, df$k)
mix = temp[[1]]
cyc = temp[[2]]
ant = temp[[3]]
sub = temp[[4]]
rm(temp)

# creating a color intensity based on value of sla and fsle
sub$value = scale(sub$fsle, from = 0.0, to = 1)
sub$colors = colvect(rep("black", length(sub$value)), alpha = sub$value)
sub$lon = jitter(sub$lon, factor = 3)
sub$lat = jitter(sub$lat, factor = 3)

cyc$value = scale(cyc$slaa, from = 0.65, to = 0)
cyc$colors = colvect(rep("navy", length(cyc$value)), alpha = cyc$value)
cyc$lon = jitter(cyc$lon, factor = 3)
cyc$lat = jitter(cyc$lat, factor = 3)

ant$value = scale(ant$slaa, from = 0, to = 0.65)
ant$colors = colvect(rep("firebrick", length(ant$value)), alpha = ant$value)
ant$lon = jitter(ant$lon, factor = 2)
ant$lat = jitter(ant$lat, factor = 2)

df = rbind(cyc, ant, sub)
rm(cyc, ant, sub)

# ------------------------------------------------------------------------------
# Plot function for loop animation ,not looping here.
my_map = function(day, wdmap, col_map, 
                  main = "",
                  e = c(-170, -130, 18, 36)){
  day = split(day, day$k)
  cyc = day[[1]]
  ant = day[[2]]
  sub = day[[3]]
  rm(day)
  
  #lay = matrix(c(1, 2),
  #             nrow = 1,
  #             byrow = TRUE)
  
  #nf = graphics::layout(lay,
  #                      widths = c(14, 1)) # this eh
  
  par(mar = c(4,4,2,1))
  
  plot(0,
       ylim = e[3:4],
       xlim = e[1:2],
       xlab = "",
       ylab = "",
       axes = FALSE)
  # plot fronts
  points(x = ant$lon, y = ant$lat, col = ant$colors, pch = 19, cex = 1)
  # plot fronts
  points(x = cyc$lon, y = cyc$lat, col = cyc$colors, pch = 19, cex = 1)
  # plot fronts
  points(x = sub$lon, y = sub$lat, col = sub$colors, pch = 20, cex = 0.5)
  
  axis(side = 2, las  = 2, lwd  = 1, mgp  = c(1, 0.75, 0), cex.axis = 1)
  
  axis(side = 1,
       at = c(-170, -165, -160, -155, -150, -145, -140, -135, -130, -125),
       las = 1, 
       lwd = 1, 
       mgp = c(2, 1, 0), 
       cex.axis = 1,
       col = "black",
  )
  
  title(main = main, cex.main = 1, line = 0.25, adj = 0)
  
  title(ylab = "Latitude", cex.lab = 1, line = 2.5)
  
  title(xlab = "Longitude", cex.lab = 1, line = 2.5)
  
  plot(wdmap,
       xlim = e[1:2],
       ylim = e[3:4],
       asp = 1,
       bg = "black",
       border = "black",
       col = "black",
       #wrap = c(-180, 180),
       add = TRUE)
  
  box(which = "plot", lty = "solid", lwd = 1.5, col = "black")
}

# ------------------------------------------------------------------------------
# plotting just one day
day = unique(df$time)
day = subset(df, time == day[1])
wdmap = getMap(resolution = "high")

dt = gsub("-", "", as.character(Sys.Date()))
pdf(paste("figures\\eddies_fronts_", dt, ".pdf", sep = ""),   # The directory you want to save the file in
    width = 6, # The width of the plot in inches
    height = 4,
    pointsize = 10) # The height of the plot in inches

my_map(day, wdmap, main = "Eddies and fronts: 2018-06-01",
       e = c(-165, -140, 23, 32.5))

dev.off()


















