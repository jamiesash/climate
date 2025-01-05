
# ------------------------------------------------------------------------------
### Liraries 
library(cmocean)
# library(rworldmap)
# library(rworldxtra)
library(terra)
# library(anytime)

print("`. Librares loaded.")

# ------------------------------------------------------------------------------
### Loading the dataset
cli = rast("../../../data/chl/cli_day_sum_20240118.nc")

print("2. Data loaded")

# ------------------------------------------------------------------------------
### Plotting prep.
# I may need to extract values here. I don't want to lower the suprimum.
cli = cli[[2]]

o = sd(values(cli), na.rm = TRUE)
l = min(values(cli), na.rm = TRUE) 
u = max(values(cli), na.rm = TRUE) 
infi = l 
supi = u - o

cli = clamp(cli, lower=infi, upper=supi, values=TRUE)

print("3. Clamped.")

# -------------------------------------------------------------------------------
### Plotting
wdmap  = getMap(resolution = "high")
colmap = cmocean("delta")(100)
dt = gsub("-", "", as.character(Sys.Date()))
e  = ext(cli)

# Plotting the function. 
pdf(paste("/home/jamesash/climate/figures/", "cli_day_", dt, ".pdf", sep = ""),  
    width  = 5.5, # inches
    height = 4,
    pointsize = 10) # inches

plot(cli, 
	# ylim = c(16, 40),
	# xlim = c(-175, -130),
	col  = colmap, 
	mar  = c(3.1, 3.1, 2.1, 7.1),
	plg  = list(size = c(1, 1.25)),
	range = c(-0.001, 0.001),
	ylab = "Latitude",
	xlab = "Longitude")
	#breaks = 100)
plot(wdmap,
     ylim = e[3:4],
     xlim = e[1:2],
     asp  = 1,
     bg   = "black",
     border = "black",
     col  = "black",
     add  = TRUE,
     lwd  = 1)
dev.off()

# ------------------------------------------------------------------------------

