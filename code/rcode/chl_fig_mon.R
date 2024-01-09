# ------------------------------------------------------------------------------
### Liraries and functions
library(cmocean)
library(rworldmap)
library(rworldxtra)
library(anytime)
library(terra)

print("2. Librares loaded.")

# ------------------------------------------------------------------------------
### Loading the dataset
cli = rast("/home/jamesash/blooms/data/chl_mon_xxxxxxxx.nc")

print("3. Data loaded")

# ------------------------------------------------------------------------------
### Plotting
# I may need to extract values here. I don't want to lower the suprimum.
 
o = sd(cli, na.rm = TRUE)
l = min(cli, na.rm = TRUE) 
u = max(cli, na.rm = TRUE) 
infi = l # + o*2
supi = u - o*2

cli = clamp(cli, lower=infi, upper=supi, values=TRUE)
print("5. Clamped.")

# -------------------------------------------------------------------------------
wdmap <- getMap(resolution = "high")
colmap = cmocean("delta")(100)
dt = gsub("-", "", as.character(Sys.Date()))
e = ext(cli)

# Plotting the function. 
pdf(paste("/home/jamesash/blooms/figures/", "cli_mon_", dt, ".pdf", sep = ""),  
    width = 5.5, # inches
    height = 4,
    pointsize = 10) # inches

plot(cli[[2]], 
	# ylim = c(16, 40),
	# xlim = c(-175, -130),
	col = colmap, 
	mar = c(3.1, 3.1, 2.1, 7.1),
	plg = list(size = c(1, 1.25)),
	range = c(-0.001, 0.001),
	ylab = "Latitude",
	xlab = "Longitude")
	#breaks = 100)
plot(wdmap,
     ylim = e[3:4],
     xlim = e[1:2],
     asp = 1,
     bg = "black",
     border = "black",
     col = "black",
     add = TRUE,
     lwd = 1)
dev.off()

# ------------------------------------------------------------------------------

