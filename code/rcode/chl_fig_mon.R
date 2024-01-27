# ------------------------------------------------------------------------------
### Liraries and functions
library(cmocean)
library(anytime)
library(terra)

print("2. Librares loaded.")

# ------------------------------------------------------------------------------
### Loading the dataset
states = vect("../../data/states/s_05mr24.shp")
cli = rast("../../data/cli_mon_sum_20240118.nc")
cli = cli[[2]]

print("3. Data loaded")
# ------------------------------------------------------------------------------
### Plotting
# I may need to extract values here. I don't want to lower the suprimum.

# for dayly I use c(-0.00005, 0.00005)
bounds = c(-0.001, 0.001)
cli = clamp(cli, lower=bounds[1], upper=bounds[2], values=TRUE)
print("5. Clamped.")

# -------------------------------------------------------------------------------
colmap = cmocean("delta")(100)
dt = gsub("-", "", as.character(Sys.Date()))
e = ext(cli)

# Plotting the function. 
pdf(paste("../../figures/", "cli_mon_", dt, ".pdf", sep = ""),  
    width = 5.5, # inches
    height = 4,
    pointsize = 8) # inches

image(cli,
      col = colmap,
      ylim = c(16, 40),
      #zlim = c(-0.00005, 0.00005),
      xlim = c(-175, -130),
      asp = 0.83,
      axes = FALSE)
plot(states,
     col = "black",
     #asp = 1,
     add = TRUE)
axis(side = 2, 
     las = 2, 
     lwd = 1, 
     at = c(18, 22, 26, 30, 34, 38),
     mgp = c(1, 0.75, 0), 
     cex.axis = 1)
axis(side = 1, 
     las = 1, 
     lwd = 1,
     mgp = c(2, 1, 0),    
     at = c(-170, -160, -150, -140, -130),
     cex.axis = 1)
box(lwd = 1)
dev.off()


pdf(paste("../../figures/", "cli_mon_", dt, ".pdf", sep = ""),  
    width = 5.5, # inches
    height = 4,
    pointsize = 8) # inches
plot(cli, 
     col = colmap,
     range=c(-0.001, 0.001),
     mar=c(4, 4, 4, 6), 
     #pal=list(shrink=0.9, cex=.8), 
     pax=list(cex.axis=1.25))
plot(states,
     col = "black",
     add = TRUE)
dev.off()


# ------------------------------------------------------------------------------

