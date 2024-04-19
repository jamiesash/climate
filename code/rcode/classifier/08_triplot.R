# USING PREDOWNOADED DATA

# load xyz table ---------------------------------------------------------------
source("code\\functions.R")
source("code\\libraries.R")
library(Thermimage)
rasterOptions(maxmemory = 123e+10)

xyz = read.csv("data\\collated\\xyz_2010_2021_20230218.csv")

gc()

# ------------------------------------------------------------------------------
# functions to prepare data to be plotted 
scale <- function(x, to, from){   
  (x - min(x))/(max(x)-min(x)) * (to - from) + from
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
      # this could be changed to nearest euclidean distance
      box_x <- which(abs(x - xy_grid[[1]][iy, ix]) <= xreach)
      box_y <- which(abs(y - xy_grid[[2]][iy, ix]) <= yreach)
      # I think the grids are dif sizes and should be subet differently
      #index vector of both cox_sla and box_fsle as one
      #box <- sort(match(box_y, box_x))
      
      # I think this is the correct way to do this
      box <- box_x[box_x %in% box_y]
      
      # z_grid[iy, ix] <- mean(z[box], na.rm = TRUE)
      z_grid[iy, ix] <- median(z[box], na.rm = TRUE)
    }
  }
  
  list(z_grid, y_vec, x_vec)
}

# ------------------------------------------------------------------------------
# making a matrix for the contour plot real data

temp = subset(xyz, !is.na(chla))

idx = sample(1:nrow(temp), 100000)
lil = temp[idx,]

lil$fsle = lil$fsle * -1

lil$fsle =  scale(lil$fsle, from = 0, to = 1)
lil$slaa =  scale(lil$slaa, from = 0, to = 1)

# lil$chla = lil$chla - median(lil$chla)

tmat <- tsect(x = lil[,'slaa'], 
              y = lil[,'fsle'], 
              z = lil[,'chla'], 
              xreach = 0.01, 
              yreach = 0.01,
              xlen = 100, 
              ylen = 100)

yvect <- tmat[[2]]
xvect <- tmat[[3]]
# I need to do this better using base R
zmat <- apply(tmat[[1]], 2, rev)

# zmat <- flip.matrix(t(zmat))
zmat = t(zmat)

# ------------------------------------------------------------------------------
# plotting the data as a contour plot
col <- colorRampPalette(c("purple3",  "blue", "cyan", "white", "yellow", "orangered", "red3"))(21)

zlim = c(median(lil$chla)-mad(lil$chla)*2, median(lil$chla) + mad(lil$chla)*2)

dt = gsub("-", "", as.character(Sys.Date()))
pdf(paste("figures\\triplot_", dt, ".pdf", sep = ""),   # The directory you want to save the file in
    width = 6, # The width of the plot in inches
    height = 4,
    pointsize = 10) # The height of the plot in inches

drawPalette(zlim = zlim,
            col  = col, 
            plot = TRUE,
            pos  = 4) 

par(mar = c(5, 5, 3, 4))
image(xvect, yvect, zmat,
      #ylim = c(max(yvect)-0.01, 0.3),
      #xlim = c(-0.08, 1.08),
      zlim = zlim,
      main = "",
      xlab = "SLA",
      ylab = "FSLE",
      col = col)
box(which = "plot", lty = "solid", lwd = 1.5, col = "black")

dev.off()

# ------------------------------------------------------------------------------
# as a veroni diagram
temp = subset(xyz, !is.na(chla))

# increase after test
idx = sample(1:nrow(temp), 1000)
lil = temp[idx,]
lil$fsle = lil$fsle * -1

# if k does not exist ill need to find it
color = c("white","grey32", "firebrick", "navy")
color = factor(lil$k, labels = color)
color = as.character(color)

dt = gsub("-", "", as.character(Sys.Date()))
pdf(paste("figures\\veroni_", dt, ".pdf", sep = ""),   # The directory you want to save the file in
    width = 6, # The width of the plot in inches
    height = 4,
    pointsize = 10) # The height of the plot in inches

plot(0, 
     0.2, 
     ylim = c(0, 1), 
     xlim = c(-0.25, 0.25), 
     col = "white", 
     pch = 18,
     cex = 0.4,
     xlab = "",
     ylab = "",
     axes = FALSE)
grid(nx = NULL, 
     ny = NULL,
     lty = 2,      # Grid line type
     col = "grey69", # Grid line color
     lwd = 1)
points(lil$sla, 
       lil$fsle,
       col = color,
       pch = 18,
       cex = 0.4)
title(xlab = "SLA [m]", line = 2.5, cex.lab = 1)
title(ylab = "FSLE [day]", line = 2.5, cex.lab = 1)
axis(side = 2, 
     las = 3, 
     lwd = 1, 
     mgp = c(1, 0.75, 0), 
     cex.axis = 1)
axis(side = 1, 
     las = 1, 
     lwd = 1,
     mgp = c(2, 1, 0),    
     cex.axis = 1)
# Grid line width
box(which = "plot", lty = "solid", lwd = 1.5, col = "black")
legend(x = 0.15, y = 0.9, 
       legend = c("Mix", "Front","Cyclon", "Anti-cyclon"),
       pch = 22,
       pt.bg = c("white","grey32", "navy", "firebrick"))

dev.off()








