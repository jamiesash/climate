
# stlfilter ---------------------------------------------------------------------
# stlfilter() Performs an STL filter on raster object by looping through each 
# grid cell. Filters each girded cell across time and either outputs the 
# seasonal or the residual component as a raster stack. Has a confusing error 
# catch in the for loop. Raster object must have time z value. Output will have
# a Z value as well. Inputs one raster stack, outputs one raster stack. 
# Arguments...
#    ras:   The raster stack object to be filtered. Must have a Z parameter set
#    np:    The length of the period to be seasonal filtered. Default is 91 for  
#           4 day temporal resolution as 91*4 = 364 days in a year. 
#    subst: Short for sub.start argument of prewriten stlplus function. "which 
#           element of sub.labels does the series begin with" - Defoult is 1, 
#           meaning our cycle starts on the first day of a year (for yearly 
#           trends), and our data must start on the 1st day of the year as well. 
#    ot:    Short for the outer argument of the stlplus() function. integer; 
#           "the number of outer robustness iterations. Default is 0, but 
#           recommended if outliers are present." Higher outer value generates 
#           residuals that look more like the input data I think. 
#    wind:  Short for the s.window window argument of the stlplus function. 
#           "either the character string "periodic" or the span (in lags) of the 
#           loess window for seasonal extraction, which should be odd. This has 
#           no default." Length of the loess running window
#    data:  Logical TRUE/FASLE to denote if the resudals or the seasonal data
#           should be saved and output. Should use a string probably
#    rem:   Should the remainder be saved? If false the seasonal cycle is saved.
stlfilter = function(ras, np = 91, subst = 1, ot = 4, wind = 7, rem = TRUE) {
  ## set the memory to handle larger matrices
  rasterOptions(maxmemory = 14e+10)
  ## initialize an empty vector of the same size.
  s     <- dim(ras) ## dimensions
  sdate <- getZ(ras) ## fake datetime vector for stl input
  y     <- array(rep(NA, s[1]*s[2]*s[3]), c(s[1], s[2],s[3])) ## the empty array
  na_id <- is.na(ras)
  ras   <- approxNA(ras)
  for (i in 1:s[1]) {
    if (i == round((s[1]/5)))   print("So much more to go")
    if (i == round(s[1]/2))     print("You're half way there")
    if (i == round((s[1]/5)*4)) print("Last leg")
    for (j in 1:s[2]) {
      tryCatch({
        stlDaily <- stlplus(x         = ras[i,j,], # One time series of the raster
                            t         = sdate, # datetime vector
                            n.p       = np, # give the period of seasonality
                            s.window  = wind, # length of the running window
                            sub.start = subst,
                            outer     = ot) 
        if (rem == FALSE) y[i,j,]  <- stlDaily$data$seasonal  #keep the seasonal
        if (rem == TRUE)  y[i,j,]  <- stlDaily$data$remainder #keep the residual
      }, error=function(e){
        #skips to next iteration if there's an error  
      }) 
    } 
    if (i == max(s[1])) print("Finish!")
    }
  y <- brick(y)
  y[na_id]  <- NA #X and Y are not same dim()
  y <- setZ(y, sdate, name="time")
  extent(y) <- extent(ras)
  y  
} 

# bloomfreq --------------------------------------------------------------------
# Input is a Rasterstack object with a set Z value as a Datetime variable. 
# Output is a vector of 12 values one for each month. Mean value of chl in each 
# month, not frequency. Should be done as a data frame with a separate row of 
# months for each value. This does take the mean of means (grand mean) so should 
# be careful as it may not be accurate. Finds months equal to the month number 
# in a loop then averages them. i in the loop is a string which I love.
# Arguments:
#    x: Raster stack with a Datetime value as the Z value. 
#       Should have at least one value for each month, otherwise an error given
bloomfreq = function(x, func = "mean", valname = "chl") {
  themonths <- c("January",
                 "February", 
                 "March", 
                 "April", 
                 "May",
                 "June", 
                 "July",
                 "August", 
                 "September", 
                 "October", 
                 "November",
                 "December")
  sdate <- getZ(x)
  m     <- months(sdate) 
  clim  <- rep(0, 12)
  j     <- 0
  # looping through each the months and finding the average in the data set
  for (i in themonths) {
    j         <- j + 1
    idx_m     <- which(m==i)
    Z         <- subset(x, idx_m)
    clim[j]   <- mean(cellStats(Z, stat = func, na.rm = TRUE), na.rm = TRUE)
  }
  clim <- data.frame(clim ,c("January", "February", "March", "April", "May",
                             "June", "July", "August", "September", "October", 
                             "November", "December"))
  colnames(clim) <- c(valname, "month")
  clim
  }

# bloomclim --------------------------------------------------------------------
# Input is a Rasterstack object with a set Z value as a Datetime variable. 
# Output is a rasterstack with 12 raster layers  one for each month. To be 
# ploted using levelplot or mylevelplot.
# Arguments:
#    x: Raster stack with a Datetime value as the Z value. 
#       Should have at least one value for each month, otherwise an error given
# bloomclim = function(x) {
#   themonths <- c("January","February", "March", "April", "May","June",  "July",
#                 "August", "September", "October", "November","December")
#  sdate <- getZ(x) # pulling out the Z value to be worked on
#  sdate <- anydate(sdate)
#  m     <- months(sdate) # getting the months from the Z value
#  j     <- 0 # start j at 0 to count loops. i is a string
#  clim  <- stack() #initialize an empty raster stack
#  # looping through each month and finding a spacial average using calc()
#  for (i in themonths) {
#    j      <- j + 1
#    idx_m  <- which(m == i)
#    z      <- subset(x, idx_m)
#    #Z      <- calc(Z, mean, na.rm = TRUE)
#    z      <- calc(z, median, na.rm = TRUE)
#    clim   <- addLayer(clim, stack(z))
#  }
#  # Should add a time stamp and extent to the raster
#  clim
#  }

# subsum -----------------------------------------------------------------------
# Subsets the late summer months from a rasterstack and outputs a rasterstack. 
# Input is a rasterstack with Z values as a Datetime variable. Output is a 
# rasterstack comprised of only the late summer months of the input data set
# one value for each month must be included in the input rasterstack for
# function to work.  
# Arguments...
#    x: Rasterstack with Z values as a Datetime variable. 
subsum    = function(x, mnths = 7:10) {
  sdate     <- getZ(x)
  themonths <- c("January","February",
                 "March", "April",
                 "May","June", 
                 "July","August",
                 "September", "October",
                 "November","December")
  m     <- months(sdate) 
  idx_t <- which(is.element(m, themonths[mnths]))
  x     <- subset(x, idx_t)
  sdate <- sdate[idx_t] 
  x     <- setZ(x, sdate, name = "time")
  x
}

# mycrop -----------------------------------------------------------------------
# Crops a raster stack or layer across a given extent, and timespan. Input is a 
# rasterstack or raster layer with Z value as a Datetime. Output is a 
# rasterstack with an updated extent. For now the region is given as a string,
# but should probably be given as an extent. ALOHA and 30N extent is included in
# the function, but should probably be outside of the function. I honestly don't
# know if this function works properly. It is messy. Looks like it was ported 
# over from one of my matlab functions. 
# Arguments...
#    X:      rasterstack object with Z value as a Daatetime
#    lons:   longitude values as separate variable. This is some matlab shit
#    lats:   latitude values as separate variables, this is some matlab shit
#    sdate:  datetime variable as a separate object. This is some matlab shit
#    start:  time-stamp to begin the cut
#    stop:   times-tamp to end the cut
#    region: either "ALOHA" or anything not "ALOHA" recognized as FALSE
mycrop    = function(X, lons, lats, sdate, start, stop, region) {
  if (string(region) == "ALOHA") {
    i = 1 
  } else {
    i = 2
  }
  
  lat1 = c(c(22,   27),  
           c(27.5, 35))
  lon1 = c(c(-159, -152), 
           c(-159, -132))
  lat  = lat1[i, ]
  lon  = lon1[i, ]
  
  test <- crop(X, extent(lat, lon))
  
  lon = lon + 360 #lons not in correct format
  index_lon = which(lons > lon(1) & lons < lon(2))
  index_lat = which(lats > lat(1) & lats < lat(2))
  
  #index by time of bloom
  start = double(datenum(start)) 
  stop  = double(datenum(stop)) 
  time2 = double(sdate)
  ind_t = find(time2 >= start  & time2 <= stop)
  
  field = X[index_lat, index_lon, ind_t]; #index climotology
  
  p = field;
}

# clean ------------------------------------------------------------------------
# This is a heavy function that can take up a lot of RAM if it's allowed. 
# Removes repeating layers. Keep first iteration of each unique time-stamp. 
# Fills in missing dates. Used to keep all raster's uniform so that indexing and 
# filters work as they should. Still converts a raster to an array which is a 
# bulky computation. 
# Arguments...
#    ras: Rasterstack object with a time Z value as a Datetime variable 
clean     = function(ras) {
  # save the extent
  e     <- extent(ras)
  
  # Remove repeating layers. Keep first iteration of each unique time-stamp.
  ts    <- getZ(ras)
  ts    <- as.Date(substr(ts, 1, 10)) # remove timestamp but keep day year month
  ind_t <- which(!duplicated(ts))
  ras   <- subset(ras, ind_t)
  ts    <- ts[ind_t]
  
  # filling in missing dates
  odate <- seq(from = min(ts), to = max(ts), by=1)
  found <-  is.element(odate, ts)
  miss  <- !is.element(odate, ts)
  fdate <- odate[found]
  mdate <- odate[miss]
  
  if (length(mdate) > 1) {
    s     <- as.numeric(dim(ras))
    y     <- array(rep(NA, s[1]*s[2]*length(mdate)), 
                   c(s[1], s[2], length(mdate)))
    y     <- brick(y)
    extent(y) <- extent(ras)
    ras   <- stack(y, ras)
    temp  <- c(fdate, mdate)
    
    idx   <- sort.int(temp, index.return = TRUE)$ix
    dates <- sort.int(temp, index.return = TRUE)$x
    ras   <- ras[[idx]]
    #inputting the dates in Z field and printing raster
    ras   <- setZ(ras, dates)
    } else {
      ras <- setZ(ras, ts)
    }
  
  extent(ras) <- e
  ras
  }

# vect -------------------------------------------------------------------------
# Takes a rasterstack and converts it to a long format data frame using 
# rastertoPoints and melt functions. This is to prep for a cluster analysis. 
# Arguments...
#    x: rasterstack with the Z values as Datetime variables
vect      = function(x) {
  sdate <- getZ(x)
  x     <- rasterToPoints(x)
  x     <- data.frame(x)
  colnames(x) <- c("lat", "lon", as.character(sdate))
  # may need to use data.table::melt()
  x     <- reshape2::melt(x, id.vars = c("lat", "lon"))
}
# map --------------------------------------------------------------------------
# Maps a single raster layer using imagerp function and worldmap. Plots land as 
# black images. Should be able to plot anywhere in the ocean using the same lat
# lon box. 
# Arguments...
#     val:   Raster layer that has a true lat lon extent, not 0-1
#     caxis: color axis cut offs of the plot. 
#     title: title of the main plot
#     cuts:  How many cuts to use in the continous scale?
#     name:  Plot title
#     col:   CMOCEAN color pallet - algae, tempo, topo, haline, delta
#     clab:  Label for the colormap
map   = function(val, 
                  zlim = c(min(raster::values(val), na.rm=TRUE), 
                           max(raster::values(val), na.rm=TRUE)), 
                  cuts = 25, 
                  col  = 'haline',
                  main = "The Title",
                  mai = c(0.01, 0.01, 0.01, 0.1), 
                  line = 1,
                  adj  = 0,
                  xlab = "Longitude",
                  ylab = "Latitude",
                  zlab = "CHL",
                  xaxes = TRUE,
                  yaxes = TRUE,
                  cex.axis = 1.5,
                  cex.lab  = 1.25,
                  cex.main = 1.5,
                  cex.zlab = 1,
                  plot = TRUE,
                  colorbar = TRUE,
                  add  = FALSE,
                  mar = c(5, 5, 5, 2)) {

  
  wdmap <- getMap(resolution = "high")
  s     <- dim(val)
  e     <- extent(val)
  lons  <- seq(from = e[1], to = e[2], length = s[1])
  lats  <- seq(from = e[3], to = e[4], length = s[2])
  val   <- raster::flip(val, direction = "x")
  val   <- raster::clamp(val, zlim[1], zlim[2])
  
  par(mar = mar)
  if(length(col) == 1) {col <- cmocean(col)(cuts)} 
  if(colorbar == TRUE) {
    drawPalette(zlim = zlim, 
                #mai = mai,
                col  = col, 
                plot = plot,
                pos  = 4,
                zlab = zlab,
                #drawContours = TRUE,
                cex  = cex.zlab) 
  }
  image(lons, lats, as.matrix(val),
        col  = col, # algae, tempo, topo, haline, delta, 
        xlim = e[1:2],
        ylim = e[3:4],
        zlim = zlim,
        main = "",
        xlab = "",
        ylab = "",
        asp = NA,
        axes = FALSE,
        cex.lab  = cex.lab,
        cex.axis = cex.axis,
        add = add)
  title(main = list(main), line = line,  adj = adj, cex.main = cex.main)
  title(xlab = xlab, line = 2.6, cex.lab = 1.75)
  title(ylab = ylab, line = 2.6, cex.lab = 1.75)
  if(yaxes == TRUE) axis(side = 2, 
                         las = 2, 
                         lwd = 2, 
                         at = c(22, 26, 30, 34),
                         mgp = c(1, 0.75, 0), 
                         cex.axis = cex.axis)
  if(xaxes == TRUE) axis(side = 1, 
                         las = 1, 
                         lwd = 2,
                         mgp = c(2, 1, 0),    
                         at = floor(seq(from = e[1]+2, to = e[2]-2, length = 4)),
                         cex.axis = cex.axis)
  box(which = "plot", lty = "solid", lwd = 3, col = "grey25")
  plot(wdmap, 
       xlim = e[1:2], 
       ylim = e[3:4], 
       asp = 1, 
       bg = "white", 
       border = "black", 
       col = "ivory", 
       #wrap = c(-180, 180), 
       add = TRUE)
}
# Line equations  --------------------------------------------------------------
# Equations that take a model summary output and convert the equation to text. 
# Output is a character string of the linear equation. 
# For displaying hte slop and intercept on graphs. 
# Alot of different ones are given because they are not smart functions. 
# Functions...
#    lm0_eqn: 0 intercept linear regression
#    lm_eqn: linear regression
#    nls_eqn: non least squares linear regression.
# Arguments...
#    x:   vector of x values
#    y:   vector of y values
#    dat: data frame x and y come from
#    srt: start values for nls only. m = slope, b = intercept

lm0_eqn <- function(x, y){
  m <- lm(y ~ 0 + x);
  eq <- substitute(italic(y) == a %.% italic(x)*","~~italic(r)^2~"="~r2, 
                   list(a  = format(unname(coef(m)[1]), digits = 2),
                        r2 = format(summary(m)$r.squared, digits = 2)))
  as.character(as.expression(eq));
}

lm_eqn <- function(slope, intercept, r){
  eq <- substitute(italic(y) == b + m %.% italic(x)*","~~italic(r)^2~"="~r2, 
                   list(b  = intercept,
                        m  = slope,
                        r2 = r))
  as.character(as.expression(eq));
}

nls_eqn <- function(x, y, dat, srt = list(m = 1, b = 0.05)){
  m <- nls(y ~ x * m + b, data = dat, start = srt)
  eq <- substitute(italic(y) == a + b %.% italic(x)*","~~italic(r)^2~"="~r2, 
                   list(a  = format(unname(coef(m)[1]), digits = 2),
                        b  = format(unname(coef(m)[2]), digits = 2),
                        r2 = format(sqrt(summary(m)$sigma), digits = 3)))
  as.character(as.expression(eq));
}

# triplot ----------------------------------------------------------------------
# Generates a triangle plot using ggplot scatter plot function. Not actually 
# how it's done in the Guo paper. 
triplot = function(tbl, main, small = FALSE) {  
  # Removing values below FSLE detection limit
  #xyz <- subset(xyz, !xyz$fsle == 0) 
  tbl <- subset(tbl, !is.na(xyz$chl)) 
  
  # Downsizing data frame to test plotting
  if (is.numeric(small)) {
    ind <- sample(1:nrow(tbl), small)
    tbl <- tbl[ind,]
  }
  
  # Creating breaks and color scale manually
  clow  <- round(median(tbl$chl, na.rm = TRUE) - sd(tbl$chl) * 2, 2)
  chigh <- round(median(tbl$chl) + sd(tbl$chl) * 2, 2)
  tbl$chls <- oob_squish(tbl$chl, range = c(clow, chigh), only.finite = TRUE)
  steps <- round(seq(from = clow, to = chigh, length = 20), 2)
  # make mid values more transparent
  trans <- abs(scales::rescale(tbl$chls, to = c(-1, 1))) 
  
  ggplot(data = tbl, aes(sla, fsle*-1, col = chls)) +
    geom_point(stroke = 0, alpha = trans, size = 0.1) +
    geom_vline(xintercept = mean(tbl$sla, na.rm=TRUE), size = 1, alpha = 0.5) +
    geom_vline(xintercept = mean(tbl$sla, na.rm=TRUE) + sd(tbl$sla, na.rm=TRUE),
               alpha = 0.5,
               linetype = "dashed",
               size = 1) +
    geom_vline(xintercept = mean(tbl$sla, na.rm=TRUE) - sd(tbl$sla, na.rm=TRUE), 
               linetype = "dashed",
               alpha = 0.5,
               size = 1) +
    geom_hline(yintercept = mean(tbl$fsle*-1, na.rm=TRUE), 
               size = 1, 
               alpha = 0.5) +
    scale_colour_steps2(low = "darkslategrey", 
                        mid = "white", 
                        high = muted("coral"),  
                        breaks = steps,
                        midpoint = mean(tbl$chl, na.rm = TRUE)) +
    scale_y_continuous(expand = c(0,0)) + 
    guides(colour = guide_coloursteps(show.limits = TRUE)) +
    xlab("SLA [m]") +
    ylab("FSLE [day^-1]") +
    labs(title = main, colour = "CHL [mg/m^3]") +
    #geom_text(x=-0.05, y= 0.27,  label= "(-) Sub") + 
    #geom_text(x= 0.2,  y= 0.27,  label= "(+) Sub") + 
    #geom_text(x=-0.05, y= 0.03,  label= "(-) Meso") + 
    #geom_text(x= 0.2,  y= 0.03,  label= "(+) Meso") + 
    theme_bw() +  
    theme(legend.key.size = unit(2, "cm"), legend.key.width = unit(1,"cm"))
}

# thresher ---------------------------------------------------------------------
# Input is a table with columns sla, and fsle. Thresholds are calculated. 
# Output is a row with percent of points above that threshold
thresher <- function(tbl, mix = TRUE, anom = TRUE) {
  if (mix == TRUE) {
    # remove clouds from chl
    tbl   <- subset(tbl, !is.na(tbl$chl))
    
    eddy  <-   sd(tbl$sla,  na.rm = TRUE)
    usla  <- mean(tbl$sla,  na.rm = TRUE)
    ufsle <- mean(tbl$fsle, na.rm = TRUE)
    uchl  <- median(tbl$chl, na.rm = TRUE)
    
    if (anom == TRUE) tbl   <- subset(tbl, chl > uchl)
    
    mesop <- tbl[tbl$sla > (usla + eddy) & tbl$fsle < ufsle,]
    meson <- tbl[tbl$sla < (usla - eddy) & tbl$fsle < ufsle,]
    subm  <- tbl[tbl$sla < (usla + eddy) & tbl$sla < (usla - eddy) & tbl$fsle > ufsle,]
    mixn  <- tbl[tbl$sla > (usla - eddy) & tbl$sla < (usla + eddy) & tbl$fsle < ufsle,]
    mixp  <- tbl[tbl$sla > (usla - eddy) & tbl$sla < (usla + eddy) & tbl$fsle > ufsle,]
    
    mixp  <-  (nrow(mixp) / nrow(tbl)) * 100 
    mixn  <-  (nrow(mixn) / nrow(tbl)) * 100 
    subm  <-  (nrow(subm) / nrow(tbl)) * 100 
    mesop <- (nrow(mesop) / nrow(tbl)) * 100
    meson <- (nrow(meson) / nrow(tbl)) * 100
    
    df_mix    <- data.frame(mesop, meson, subm, mixn+mixp)
    colnames(df_mix) <- c("mesop", "meson", "subm", "mix")
    df_mix
  }
  
  if (mix == FALSE) {
    
    # remove clouds from chl
    tbl   <- subset(tbl, !is.na(tbl$chl))
    
    usla  <- mean(tbl$sla,  na.rm = TRUE)
    ufsle <- mean(tbl$fsle, na.rm = TRUE)
    uchl  <- median(tbl$chl, na.rm = TRUE)
    
    if (anom == TRUE) tbl <- subset(tbl, chl > uchl)
    
    mesop <- tbl[tbl$sla  > usla  & tbl$fsle < ufsle,]
    meson <- tbl[tbl$sla  < usla  & tbl$fsle < ufsle,]
    subp  <- tbl[tbl$sla  > usla  & tbl$fsle > ufsle,]
    subn  <- tbl[tbl$sla  < usla  & tbl$fsle > ufsle,]
    
    subp  <-  (nrow(subp) / nrow(tbl)) * 100 
    subn  <-  (nrow(subn) / nrow(tbl)) * 100
    mesop <- (nrow(mesop) / nrow(tbl)) * 100
    meson <- (nrow(meson) / nrow(tbl)) * 100
    
    df_no <- data.frame(mesop, meson, subp, subn)
    colnames(df_no) <- c("mesop", "meson", "subp", "subn")
    df_no
  }
}

# tableit ----------------------------------------------------------------------
# Uses the previously defined vect function to generate a long format table from
# raster stack objects. Not sure if this should be using teh vector function now. 
# x: a raster or list of rasters as c(x, y, z) where x, y, and z are rasters 

tableit <- function(x) {
  
  vectorize  <- function(x) {
    sdate <- getZ(x)
    x     <- rasterToPoints(x, na.rm = FALSE)
    x     <- data.frame(x)
    
    colnames(x) <- c("lat", "lon", as.character(sdate))
    x     <- reshape2::melt(x, id.vars = c("lat", "lon"))
    
    colnames(x) <- c("lats", "lons", "time", "val")
    x
  }
  
  x     <- lapply(x, vectorize)
  vals  <- lapply(x, subset, select = "val")
  quord <- lapply(x, subset, select = c("lats", "lons", "time"))
  
  xyz   <- data.frame(quord[1], vals)
  #colnames(xyz) <- c("lats", "lons", "time", "sla", "fsle", "chl")
  
  # removing chl clouds
  #xyz   <- xyz[!is.na(xyz$chl), ]
  xyz
}

# crop3d -----------------------------------------------------------------------
# crop 3d and timesnip seem redundant. I also have a my crop function. 
crop3d  <- function(x, 
                    sdate = "2003-07-01",
                    edate = "2003-10-01", 
                    ext = extent(x)) {
  # subset time
  #if(!is.null(sdate)){
    t <- getZ(x)
    ind_t <- which(t > sdate & t < edate)
    indx_bool <- (t > sdate & t < edate)
    x <- subset(x, ind_t)
    t <- subset(t, indx_bool)
  #}
  
  # subset space
  #if(!is.null(ext)){
    ext <- extent(ext)
    x   <- raster::crop(x, ext)
    x   <- setZ(x, z = t, name= "time")
    x
   # }
}

# ------------------------------------------------------------------------------
timesnip = function(x, 
                     sdate = "2003-07-01",
                     edate = "2003-10-01") {
  # subset time
  t = getZ(x)
  ind_t = which(t >= sdate & t <= edate)
  indx_bool = (t >= sdate & t <= edate)
  x = subset(x, ind_t)
  t = subset(t, indx_bool)
  x = setZ(x, z = t, name = "time")
  x
  }

# findtime ---------------------------------------------------------------------
# Input is a rasterbrick
# out put is a percent each of time spent in a bloom state. 
# Angel wants this value!!!

findtime <- function(x, thresh = 0.15) {
  y <- x
  raster::values(y)[is.na(raster::values(y))]   <- 0 
  if (is.numeric(thresh)) {
    raster::values(y)[raster::values(y) > thresh] <- 1 
    raster::values(y)[raster::values(y) < thresh] <- 0 
  } else {
    raster::values(y)[raster::values(y) > raster::values(thresh)] <- 1 
    raster::values(y)[raster::values(y) < raster::values(thresh)] <- 0 
  }
  blooms <- calc(y, sum, na.rm = TRUE)
  
  y <- x
  raster::values(y)[!is.na(raster::values(y))] <- 0 
  raster::values(y)[is.na(raster::values(y))]  <- 1 
  clouds   <- calc(y, sum) #na.rm should deal with clouds
  days     <- calc(y, length) #na.rm should deal with clouds
  sunydays <- days - clouds
  
  (blooms/sunydays)*100
}

# boot -------------------------------------------------------------------------
boot <- function(tbl, iter = 10, sigma = FALSE, mu = FALSE){
  # remove clouds from chl
  df   <- subset(tbl, !is.na(tbl$chl))
  sz   <- nrow(subset(df, chl > 0)) # size of sub-samples = total anomalies
  
  # creating an empty data frame for each sub-sample iteration
  df_no <- data.frame(rep(NA, iter), rep(NA, iter), rep(NA, iter), rep(NA, iter))
  colnames(df_no) <- c("mesp", "mesn", "sub", "mix")
  for (i in 1:iter) {
    ind   <- sample(1:nrow(df), sz, replace = FALSE)
    anom   <- df[ind, ] # treate these like the anomolies? or find anomolies again  
    #if (anom == TRUE) tbl <- subset(tbl, chl > 0) 
    sub  <- subset(anom, k == 1)
    mesp <- subset(anom, k == 4)
    mesn <- subset(anom, k == 2)
    mix  <- subset(anom, k == 3)
    sub  <- (nrow(sub)  / nrow(anom)) * 100 
    mesp <- (nrow(mesp) / nrow(anom)) * 100 
    mesn <- (nrow(mesn) / nrow(anom)) * 100 
    mix  <- (nrow(mix)  / nrow(anom)) * 100 
    df_no[i, ] <- data.frame(mesp, mesn, sub, mix)[1, ] 
  } 
  # sd and mean of each column. 
  if (sigma == TRUE) test <- apply(df_no, 2, FUN = sd,   na.rm = TRUE)
  if (mu    == TRUE) test <- apply(df_no, 2, FUN = mean, na.rm = TRUE) 
  test
}

# kperc ------------------------------------------------------------------------ 
kperc <- function(df, anom = TRUE){
  tbl   <- subset(df, !is.na(df$chl))
  if (anom == TRUE) tbl <- subset(tbl, chl > 0) 
  sub  <- subset(tbl, k == 1)
  mesp <- subset(tbl, k == 4)
  mesn <- subset(tbl, k == 2)
  mix  <- subset(tbl, k == 3)
  sub  <- (nrow(sub)  / nrow(tbl)) * 100 
  mesp <- (nrow(mesp) / nrow(tbl)) * 100 
  mesn <- (nrow(mesn) / nrow(tbl)) * 100 
  mix  <- (nrow(mix)  / nrow(tbl)) * 100 
  perc <- data.frame(mesp, mesn, sub, mix)
  colnames(perc) <- c("mesp", "mesn", "sub", "mix")
  perc
}

# timeseries -------------------------------------------------------------------
# Time series of year of interest: bloom threshold
# For plotting a timeseries. 
# input is x and y vectors. 
# output is a lineplot.
timeseries <- function(x, 
                       y, 
                       xlim,
                       ylim = c(range(y, na.rm = TRUE)[1] 
                                - sd(y, na.rm = TRUE),
                                range(y, na.rm = TRUE)[2] 
                                + sd(y, na.rm = TRUE)),
                       ylab = "CHLsat Anomaly [mg/m^3]",
                       xlab = "Date",
                       thresh = median(y, na.rm = TRUE) + mad(y, na.rm=TRUE)*2,
                       main   = "Timeseries of CHLsat anomaly",
                       cex.main = 1.5,
                       line = 0.7,
                       adj = 0,
                       xaxis = TRUE,
                       yaxis = TRUE,
                       pretty = TRUE,
                       mar = c(3,3,3,3)) {
  # smooth timeseries
  if (pretty == TRUE) {
    y <- na_ma(y, k = 10, weighting = "exponential")
    y <- smooth(y)
  }
  # bloom threshold
  #thresh <- median(y) + mad(y)
  # all values above threshold are threshold
  y2 <- y
  y2[y2 < thresh] <- NA
  y2[y2 > thresh] <- thresh
  
  ylarge <- max(y, na.rm = TRUE)
  ysmall <- min(y, na.rm = TRUE)
  xlarge <- subset(x, y == ylarge)
  xsmall <- subset(x, y == ysmall)
  
  par(mar=mar)
  plot(x, y, 
       type = 'n', 
       ylim = ylim,
       xlim = xlim,
       xlab = "",
       ylab = "",
       axes = FALSE,
       main = "")
  box(which = "plot", lty = "solid", lwd = 2, col = "grey25")
  title(main = list(main, cex = cex.main),
        line= line,
        adj = adj)
  title(ylab = ylab, 
        cex.lab = 1.5,
        line = 3.5)
  title(xlab = xlab, 
        cex.lab = 1.5,
        line = 2.5)
  if(yaxis) axis(side = 2,
                 las = 2, 
                 lwd = 2, 
                 mgp = c(1, 0.75, 0), 
                 cex.axis = 1.5)
  if(xaxis) axis(side = 1, 
                 #at = as.numeric(x)[seq(from = 10, to = length(x), by = 5)],
                 #labels = labels[seq(from = 10, to = length(x), by = 5)],
                 las = 1, 
                 lwd = 2, 
                 mgp = c(2, 1, 0), 
                 cex.axis = 1.5)
  #text(xlarge - 8, ylarge, as.character(round(ylarge, 2)))
  #text(xsmall - 8, ysmall, as.character(round(ysmall, 2)))
  #points(c(xlarge, xsmall), c(ylarge, ysmall), pch = 21, bg = "grey49")
  #title(main = list(title), line = line, adj = adj)
  clip(x1 = min(x),
       x2 = max(x), 
       y1 = thresh, 
       y2 = max(y))
  polygon(c(min(x), x, max(x)), c(min(y), y, min(y)), col = "darkseagreen4")
  clip(par("usr")[1], par("usr")[2], par("usr")[3], par("usr")[4]) # reset clipping region
  lines(x, y,  type='l')
  lines(x, y2, lwd = 1)
  abline(h = thresh, col = "coral4")
}
# hotdata ----------------------------------------------------------------------
# converts the hot data text file to a data frame. 
hotdata = function(path = "..//data//infiles//hotdogs//", 
                   file = "HPLC_HOT.txt",
                   cutdate = "2002-01-01",
                   varname = "chl",
                   id = c("press", "hplc"),
                   depth = NULL){
  
  dtfile <- paste(path, file, sep = "")
  # Load the HPLC Data text file from HOT_DOGS webserver
  dt <- data.frame(fread(dtfile, 
                         header = TRUE),
                   stringsAsFactors = FALSE)
  
  dt <- data.frame(fread(dtfile, 
                         colClasses = c(rep("character", ncol(dt))), 
                         header = TRUE),
                   stringsAsFactors = FALSE)
  date <- as.POSIXct(dt$date, format = "%m%d%y")
  
  dt   <- data.frame(apply(dt[, id], 
                           MARGIN = 2, 
                           FUN = as.numeric))
  dt      <- cbind(dt, date) # create a dataframe
  colnames(dt) <- c(id, "time")
  dt$time <- as.Date(dt$time) # that Postix is now a date class
  if(is.numeric(depth)) dt <- subset(dt, press < depth)# subset upper 5m
  ind     <- which(dt$time > as.Date(cutdate))
  dt      <- dt[ind,] 
  dt
}

# ------------------------------------------------------------------------------
# creates a color vector for ploting that has an alpha value
colvect <- function(x = c("white", "lightblue", "coral3", "purple"), 
                    alpha = 0.5) {
  cols <- vector()
  temp <- col2rgb(x)
  for (i in 1:ncol(temp)) {
    if (length(alpha) > 1) {
      cols[i] <- rgb(temp[1,i], temp[2,i], temp[3,i], 
                     alpha = alpha[i] * 255, 
                     max = 255)
    } else {
      cols[i] <- rgb(temp[1,i], temp[2,i], temp[3,i], 
                     alpha = alpha * 255, 
                     max = 255)
    }
  }
  cols
}

# ------------------------------------------------------------------------------
shadedline = function(x, y1, y2, y0 = NULL,
                      ylim = c(min(y1, na.rm = TRUE), max(y2, na.rm = TRUE)),
                      xlim = c(min(x, na.rm = TRUE), max(x, na.rm = TRUE)),
                      ylab = "",
                      xlab = "",
                      line = -1.5,
                      adj = 0.05,
                      cex.main = 1.25,
                      main = "CHL Climatology",
                      col = "darkseagreen",
                      labels = format(x, "%b"),
                      xaxis = TRUE,
                      yaxis = TRUE,
                      mar = c(3, 3, 3, 3),
                      title = "") {
  
  mycolor <- colvect(x = col, alpha = 0.5)
  gridcol <- colvect(x = "gray", alpha = 0.6)
  m <- as.Date(x)
  
  par(mar=mar)
  plot(x, y1, 
       type = "l",
       ylim = ylim, 
       xlim = xlim,
       ylab = "",
       xlab = "",
       axes = FALSE)
  lines(x, y2, type = "l", col = 2)
  if(yaxis) axis(side = 2,
                 las = 2, 
                 lwd = 2, 
                 mgp = c(1, 0.75, 0), 
                 cex.axis = 1.15)
  if(xaxis) axis(side = 1, 
                 at = as.numeric(m)[seq(from = 1, to = length(m), by = 2)],
                 labels = format(m, "%b")[seq(from = 1, to = length(m), by = 2)],
                 las = 1, 
                 lwd = 2, 
                 mgp = c(2, 1, 0), 
                 cex.axis = 1.15)
  grid(nx = NULL, 
       ny = NULL,
       lty = 1,      # Grid line type
       col = gridcol, # Grid line color
       lwd = 1)      # Grid line width
  box(which = "plot", lty = "solid", lwd = 3, col = "grey22")
  title(main = list(main, cex = cex.main),
        line= line,
        adj = adj)
  title(main = list(title, cex = cex.main),
        line= 0.7,
        adj = 0)
  title(ylab = ylab, 
        cex.lab = 1.25,
        line = 3.5)
  title(xlab = xlab, 
        cex.lab = 1.25,
        line = 2.5)
  # Fill area between lines
  polygon(c(x, rev(x)), c(y2, rev(y1)),
          col = mycolor, lty = 0)
  
  # Redraw the lines
  lines(x, y0, col = "grey22", lty = 1, lwd = 2)
  lines(x, y1, col = col, lwd = 1)
  lines(x, y2, col = col, lwd = 1)
}

# load_edds PARELLEL -----------------------------------------------------------
load_edds = function(
    path = "..\\data\\eddies\\",
    file = "Eddy_trajectory_nrt_3.2exp_cyclonic_20180101_20220118.nc",
    sdate = as.Date("2018-10-02"),
    edate = as.Date("2018-10-02") + 10,
    domain = extent(-158, -130, 23, 35)){
  
  nc_data <- nc_open(paste(path, file, sep = ""))
  # variables <- names(nc_data$var)
  lat = ncvar_get(nc_data, "latitude")
  lon = ncvar_get(nc_data, "longitude")
  time = ncvar_get(nc_data, "time")
  age  = ncvar_get(nc_data, "observation_number")
  track = ncvar_get(nc_data, "track")
  
  clon = ncvar_get(nc_data, "effective_contour_longitude")
  clat = ncvar_get(nc_data, "effective_contour_latitude")
  
  names(nc_data$var)
  
  nc_close(nc_data)
  
  time  = unlist(time)
  track = unlist(track)
  age   = unlist(age)
  lat   = unlist(lat)
  lon   = unlist(lon)
  
  clat = c(clat)
  clon = c(clon)
  
  lon = rep(lon, each = 20)
  lat = rep(lat, each = 20)
  time = rep(time, each = 20)
  age = rep(age, each = 20)
  track = rep(track, each = 20)
  
  time  = as.Date(time, origin = as.Date("1950-01-01"))
  age   = as.numeric(age)
  track = as.numeric(track)
  lon   = lon - 360
  clon  = clon - 360
  
  edds = data.frame(time, lat, lon,  clat, clon, track, age)
  domain = extent(domain)
  edds = subset(edds, lon > domain[1] & lon < domain[2])
  edds = subset(edds, lat > domain[3] & lat < domain[4])
  edds
  }

# boxit ------------------------------------------------------------------------
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
                  line = -1.5,
                  adj = 0.05,
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
          outline = FALSE, 
          axes = FALSE)
  grid(nx = 6, # X-axis divided in two sections
       ny = 3, # Y-axis divided in three sections
       lty = 2, col = colvect(c("gray69"), alpha = 0.5), lwd = 1)
  legend(x = 20, y = 0.2, legend = legend, pch = 22,
         pt.bg = col)
  if(yaxes) axis(side = 2,
                 las = 2, 
                 lwd = 2, 
                 mgp = c(1, 0.75, 0), 
                 cex.axis = 1.5)
  if(xaxes) axis(side = 1, 
                 las = 1, 
                 lwd = 2, 
                 at = at,
                 labels = labels, 
                 mgp = c(2, 1, 0), 
                 cex.axis = 1.5)
  box(which = "plot", lty = "solid", lwd = 3, col = "grey25")
  title(main = list(sub, cex = cex.main),
        line= line,
        adj = adj)
  title(ylab = ylab, 
        cex.lab = cex.lab,
        line = 3.5)
  title(xlab = xlab, 
        cex.lab = cex.lab,
        line = 2.5)
}



# ------------------------------------------------------------------------------
plotit <- function(x, y, 
                   asp = 1,
                   abline = c(m, b, r),
                   tex_x = 0.15,
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
                   line = -1.5,
                   adj = 0.05,
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
       lty = 2, col = colvect(c("gray69"), alpha = 0.6), 
       lwd = 0.75)
  lines(xm, ym, lwd = 1,  col = "grey12")
  points(x, y, pch = pch)
  if(yaxes) axis(side = 2,
                 las = 2, 
                 lwd = 2, 
                 mgp = c(1, 0.75, 0), 
                 cex.axis = 1.5)
  if(xaxes) axis(side = 1, 
                 las = 1, 
                 lwd = 2, 
                 mgp = c(2, 1, 0), 
                 cex.axis = 1.5)
  box(which = "plot", lty = "solid", lwd = 2, col = "grey12")
  title(main = list(sub, cex = cex.main),
        line= line,
        adj = adj)
  title(ylab = ylab, 
        cex.lab = cex.lab,
        line = 3.5)
  title(xlab = xlab, 
        cex.lab = cex.lab,
        line = 2.5)
}

# ------------------------------------------------------------------------------
# This function will requir three data inputs (x, y, z) and one units input 
# (grid spacing). Possibly a padding input.
# Dates should be a numeric type
# x, y, and z should be individual nemeric vectors
# use the range of x, and y to determine the x-yunits
# units = 0.005

loadfloat <- function(fields =  c("temp",
                                  "temp_qc",
                                  "temp_adjusted", 
                                  "chla",
                                  "doxy", 
                                  "psal", 
                                  "psal_adjusted", 
                                  "pres", 
                                  "pres_adjusted", 
                                  "time", 
                                  "float_serial_no",
                                  "cycle_number", 
                                  "latitude", 
                                  "longitude",
                                  "platform_number"), 
                      sdate = "2018-07-01", 
                      edate = "2018-10-01", 
                      lons = c(-160, -145), 
                      lats = c(28, 32),
                      base_url = "https://erddap.ifremer.fr/erddap",
                      id   = "ArgoFloats-synthetic-BGC") {
  
  sdate <- paste("time>=", as.character(sdate), "T00:00:00Z", sep = "")
  edate <- paste("time<=", as.character(edate), "T00:00:00Z", sep = "")
  slon <- paste("longitude>=", lons[1], sep = "")
  elon <- paste("longitude<=", lons[2], sep = "")
  slat <- paste("latitude>=", lats[1], sep = "")
  elat <- paste("latitude<=", lats[2], sep = "")
  
  data_info <- rerddap::info(id, url = base_url)
  
  floats <- tabledap(x = data_info, 
                     fields = fields, 
                     url = base_url,
                     sdate, 
                     edate,
                     slon, 
                     elon,
                     slat, 
                     elat)
  floats      <- data.frame(floats)
  floats 
}

# ------------------------------------------------------------------------------
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
# mldchu -----------------------------------------------------------------------
# Method for determining the mixed layer by Chu and Fan (2010)
# ctd is a data frame with at least a "pressure/pres" column and variable column
# n is the length of the running window
# variable is the dependent variable that determines the mixed layer
# output is the row of the ctd data set that the mixed layer occurs at 
mldchu = function(ctd, n = 5, y = "pres", x = "temperature"){
  pressure = ctd[[y]]
  x = ctd[[x]]
  ndata = length(pressure)
  E1 = rep(NA, ndata)
  E2 = E1
  E2overE1 = E2
  kstart = min(n,3)
  for (k in seq(kstart, ndata-n,1)){
    above = seq.int(1,k)
    below = seq.int(k+1, k+n)
    fit = lm(x~pressure, subset = above)
    E1[k] = sd(predict(fit) - x[above])
    pBelow = data.frame(pressure = pressure[below])
    E2[k] = abs(mean(predict(fit, newdata = pBelow) -x[below]))
    E2overE1[k] = E2[k] / E1[k]
  }
  MLDindex = which.max(E2overE1)
  return(ctd[MLDindex,])
  # MLDindex = MLDindex, #E1 = E1, #E2 = E2
}

# ------------------------------------------------------------------------------
scale01 <- function(x) {
  (x - min(x, na.rm = TRUE)) / (max(x, na.rm = TRUE) - min(x, na.rm = TRUE))}

# anom <- function(x) x - mean(x, na.rm = TRUE)

# ------------------------------------------------------------------------------
# This is all to calculate the mixed layer depth. Not sure if it's needed
# The mixed layer depth (MLD) is detected following de Boyer Montegut et al 
# (2004): a reference value for density is taken near the surface and the water 
# column is considered to be mixed until the depth at which density deviates 
# from this reference by more than 0.03 kg/m^3. Optionally, when no depth 
# satisfies this criterion, a second criterion of 0.01 kg/m^3 can be considered 
# (as is by default). 
# In addition, here, when a range of depths is provided as reference, the 
# reference density is the average of the densities recorded in this depth range
# x:             vector of the variable of interest, usually potential or in situ density.
# depth:         vector of depths at which x is measured.
# ref.depths:    depth(s) of reference, near the surface; when ref.depths is a vector, the value of x is averaged between those depths.
# criteria:      value(s) considered as thresholds for the computation of the depth of the mixed layer. The thresholds are tried successively.
# default.depth: when no threshold is crossed, return this value instead; a usual value is the full depth of the water column (or max(depth)).
# n.smooth:      integer, number of times to smooth the data before applying the mixed layer criteria.
# k:             order of the moving average used for smoothing; the window size is 2k+1. NB: when data is smoothed, it should have been recorded at approximately regular intervals for the moving average to make sense

slide <- function(x, k, fun, n=1, ...) {
  # make sure to get a function as the `fun` argument (i.e. avoid name masking)
  if (!is.function(fun)) {
    fun <- get(as.character(substitute(fun)), mode="function")
  }
  
  if (n>=1) {
    # repeat n times
    for (t in 1:n) {
      # pad the extremities of data to be able to compute over the whole vector
      x <- c(rep(NA, times=k), x, rep(NA, times=k))
      
      # apply the rolling function (and remove padding at the extremities)
      x <- sapply((k+1):(length(x)-k), function(i) {
        fun(x[(i-k):(i+k)], ...)
      })
    }
  }
  
  return(x)
}

get_depth <- function(i, depth) {
  if (length(i) > 0) {
    if (!is.null(depth)) {
      i <- depth[i]
    }
  } else {
    i <- NA
  }
  return(i)
}

smooth <- function(x, k=1, n=1) {
  # compute centered weights
  w <- c(1:k,k+1,k:1)
  w <- w / sum(w)
  # compute the (running) weighted moving average
  slide(x, k=k, stats::weighted.mean, na.rm=TRUE, w=w, n=n)
}

check_input <- function(x, depth=NULL) {
  ok <- TRUE
  # check the input
  if (all(is.na(x))) {
    ok <- FALSE
  }
  if (!is.null(depth)) {
    if (length(depth) != length(x)) {
      ok <- FALSE
      stop("The vector of data (n=", length(x), ") should be as long as the vector of depths (n=", length(depth), ")")
    }
  }
  return(ok)
}

mld <- function(x, depth, 
                ref.depths = 5:10, 
                criteria = c(0.03, 0.01), 
                default.depth = NA, 
                n.smooth = 0, 
                k = 2) {
  # check input
  ok <- check_input(x, depth)
  if (!ok) { return(NA) }
  
  # smooth the profile (if requested)
  x <- smooth(x, k = k, n = n.smooth)
  
  # compute the reference value
  iref <- which(depth >= min(ref.depths) & depth <= max(ref.depths))
  xref <- mean(x[iref], na.rm = TRUE)
  if (is.na(xref)) {
    warning("No data at reference depth(s).")
    m <- NA
  } else {
    for (crit in criteria) {
      i <- which(x > (xref + crit) & depth > max(ref.depths)) - 1
      # NB: we want the element previous to meeting this criterion
      if (length(i) > 0) {
        i <- min(i)
        break
      }
    }
    
    # extract the corresponding depth
    m <- get_depth(i, depth)
    
    # replace by the default value when no criterion is met
    if (is.na(m)) {
      m <- default.depth
    }
  }
  return(m)
}

# anomalize --------------------------------------------------------------------
# anomalize = function(ras, detrend = FALSE, f = 0.6){
#   # find the monthly climotology of the data set 
#  ras_clim = bloomclim(ras)
   
  # subtract each month from corresponding daily data set
#   ogt      = getZ(ras)
#  ogt      = anydate(ogt)
  # I should be able to do this with week of year or day of year onj a larger
  # dataset
#  mon_raw  = month(ogt)
  #mon_clim = month(getZ(ras_clim))
  
  # subtract each day from corresponding daily data set
  #day_raw  = day(getZ(ras))
  #day_clim = day(getZ(ras_clim))
  
#  s     = dim(ras_clim)
#  chla  = stack() #initialize an empty raster stack
#  t     = vector()
  
#  for (mon in 1:s[3]) {
#    ind  = which(mon_raw == mon)
#    z    = ogt[ind]
#    temp = ras[[ind]] - ras_clim[[mon]]
#    chla = addLayer(chla, temp)
#    t    = c(t, z)
#  }
  
#  rm(temp, z)
  
#  chla = setZ(chla, t, name = "time")
#  chla = chla[[order(t)]]
#  t    = getZ(chla)
#  chla = setZ(chla, anydate(getZ(chla)), name = "time")
#  
#  if(detrend == TRUE) {
#    clim = smooth.time.series(chla, f = f, smooth.data = TRUE)
#    chla = chla - clim
#    chla = setZ(chla, z = anydate(t), name = "time")
#    }
#  chla
#}
# anom -------------------------------------------------------------------------
# x is a raster with a time atribute
# output is  araster missing leap days
# will likely take 365/12 times longer than the other function 
# x is a raster with a time atribute
# output is  araster missing leap days
# will likely take 365/12 times longer than the other function 
anom = function(x, detrend = FALSE){
  t = getZ(x)
  y = year(t)
  m = month(t) 
  d = yday(t)
  
  lyears = y %% 4 == 0
  
  idx = which(!(lyears & d == 60))
  
  # removing teh leap days
  x = x[[idx]]
  t = t[idx]
  d = d[idx]
  lyears = lyears[idx]
  
  # creating a day of year to do climitology
  idx = which(lyears & (d > 60))
  d[idx] = d[idx] - 1
  
  # just make a seq of 1:365 as the doy from then on
  anomaly = stack()
  time = vector()
  for (i in 1:365) {
    idx = which(d == i)
    z = subset(x, idx)
    t = getZ(z)
    z = calc(z, median, na.rm = TRUE)
    z = x[[idx]] - z
    
    anomaly = addLayer(anomaly, stack(z))
    time = c(time, t)
  }
  rm(z, t, idx, x)
  
  # re-order amonlay by time
  idx = order(time)
  time = as.Date(time[idx])
  anomaly = anomaly[[idx]]
  
  # set the time atribute
  anomaly = setZ(anomaly, z = time, name = "time") 
  
  if(detrend == TRUE) {
    clim = smooth.time.series(anomaly, f = 0.6, smooth.data = TRUE)
    anomaly = anomaly - clim
    anomaly = setZ(anomaly, z = time, name = "time")
  }
  anomaly
}

# polymask ---------------------------------------------------------------------
polymask <- function(ras, 
                     lon_coord = c(-170, -152, -157, -162, -170, -170), 
                     lat_coord = c(  20,   20,   22,  24.5, 28, 20)){
  t   <- getZ(ras)
  ras <- raster::flip(ras, direction = "x")
  ras <- raster::flip(ras, direction = "y")
  
  Sr1 <- Polygon(cbind(lon_coord, lat_coord))
  
  spp <- SpatialPolygons(list(Polygons(list(Sr1), "s1")))
  spp_ras <- raster::rasterize(spp, ras, getCover = TRUE)
  spp_ras[spp_ras == 1] <- NA
  
  ras <- raster::mask(ras, spp_ras)
  
  # put it back in reverse order
  ras <- raster::flip(ras, direction = "y")
  ras <- raster::flip(ras, direction = "x")
  
  ras <- setZ(ras, z = t, name = "time")
  ras
}

# dap --------------------------------------------------------------------------
# the data loading function we've all been waiting for
dap = function(sdate, edate, e, 
               id = 'jplMURSST41anom1day', 
               url = "https://coastwatch.pfeg.noaa.gov/erddap/",
               var = c("sstAnom", "latitude", "longitude", "time")) {
  # input is a datetime and an extent. It grabs that raster
  dap_that_ass = function(x, data_info, e){
    data = griddap(data_info, 
                   latitude = e[3:4], 
                   longitude = e[1:2], 
                   time = c(as.character(x),as.character(x)), 
                   fields = 'all')
    data = nc_open(data$summary$filename)
    ras  = ncvar_get(data, varid = var[1])
    lats = ncvar_get(data, varid = var[2])
    lons = ncvar_get(data, varid = var[3])
    time = ncvar_get(data, varid = var[4])
    time = as.Date(as.POSIXct(time, origin = "1970-01-01"))
    nc_close(data)
    rm(data)
    
    ras = raster(ras)
    extent(ras) = extent(min(lons), max(lons), min(lats), max(lats))
    ras = setZ(ras, z = time, name = "time")
    crs(ras) = "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"
    ras
  }
  
  time = seq(sdate, edate, by = 1)
  time = as.list(time)
  
  data_info = info(id, 
                   url = url)
  
  ras  = lapply(time, FUN = dap_that_ass, data_info = data_info, e = e)
  time = lapply(ras, getZ)
  time = unlist(time)
  time = as.Date(time)
  ext    = lapply(ras, extent)
  ras  = stack(ras)
  extent(ras) = ext[[1]]
  ras = setZ(ras, z= time, name = "time")
  ras
}
# dap ----------------------------------------------------------------------
# load errdapp data function. 
# set up for ssta.
# only good for daily data now, but probs not necissary for monthy data

dap = function(sdate, edate, e, 
               id = 'jplMURSST41anom1day', 
               url = "https://coastwatch.pfeg.noaa.gov/erddap/",
               var = c("sstAnom", "latitude", "longitude", "time")) {
  # input is a datetime and an extent. It grabs that raster
  dap_that_ass = function(x, data_info, e){
    data = griddap(data_info, 
                   latitude = e[3:4], 
                   longitude = e[1:2], 
                   time = c(as.character(x),as.character(x)), 
                   fields = 'all')
    data = nc_open(data$summary$filename)
    ras  = ncvar_get(data, varid = var[1])
    lats = ncvar_get(data, varid = var[2])
    lons = ncvar_get(data, varid = var[3])
    time = ncvar_get(data, varid = var[4])
    time = as.Date(as.POSIXct(time, origin = "1970-01-01"))
    nc_close(data)
    rm(data)
    
    ras = raster(ras)
    extent(ras) = extent(min(lons), max(lons), min(lats), max(lats))
    ras = setZ(ras, z = time, name = "time")
    crs(ras) = "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"
    ras
  }
  
  time = seq(sdate, edate, by = 1)
  time = as.list(time)
  
  data_info = info(id, 
                   url = url)
  
  ras  = lapply(time, FUN = dap_that_ass, data_info = data_info, e = e)
  time = lapply(ras, getZ)
  time = unlist(time)
  time = as.Date(time)
  ext    = lapply(ras, extent)
  ras  = stack(ras)
  extent(ras) = ext[[1]]
  ras = setZ(ras, z= time, name = "time")
  ras
}

# ------------------------------------------------------------------------------

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
  e = extent(lons, lats)
  
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
bufcoast = function(ras, 
                    region = "Hawaiian Islands", 
                    path = "../data/USMaritimeLimitsAndBoundariesSHP"){
  path.eez.usa = (path)
  fnam.eez.usa = "USMaritimeLimitsNBoundaries.shp"
  eez.usa = readOGR(dsn = path.eez.usa, layer = file_path_sans_ext(fnam.eez.usa))
  idx = eez.usa$REGION == "Hawaiian Islands"
  hawaii = eez.usa[idx,]
  #idx = hawaii$CZ == 1
  idx = hawaii$TS == 1
  hawaii = hawaii[idx,]
  hawaii = st_as_sf(hawaii)
  hawaii = st_polygonize(hawaii)
  hawaii = as(hawaii, "Spatial")
  
  ras = raster::mask(ras, hawaii, inverse = TRUE)
}

# ------------------------------------------------------------------------------
# input is a raster
fronts = function(ras,
                  downsize = 10,
                  I = 2){
  ras  = raster::flip(ras, direction = "x")
  thresh = median(ras, na.rm = TRUE) + mad(ras, na.rm = TRUE)*I
  idx  = which(values(ras) > thresh)
  e    = extent(ras)
  s    = dim(ras)
  lons = seq(from = e[1], to = e[2], length = s[1]) 
  lats = seq(from = e[3], to = e[4], length = s[2]) 
  grid = pracma::meshgrid(lons, lats)
  fronts = data.frame(lons = grid$X[idx], lats = grid$Y[idx], value = ras[idx])
  idx = seq(1, nrow(fronts), downsize)
  fronts = fronts[idx, ]
  fronts
}

# ------------------------------------------------------------------------------
# input is a raster output is a boolian raster of bloom/not bloom
bool = function(x){
  u <- calc(x, fun = median, na.rm = TRUE)
  o <- calc(x, fun = mad, na.rm = TRUE)
  boo <- x > (u + o)
  extent(boo) <- extent(x)
  boo <- setZ(boo, z = getZ(x), name = "time")
  boo
}

# oreant -----------------------------------------------------------------------

oreant = function(ras, flip = NULL, t1 = FALSE, t2 = FALSE){
  e    = extent(ras)
  time = getZ(ras)
  if(t1) ras = raster::t(ras)
  if(!is.null(flip)) ras = raster::flip(ras,  direction = flip)
  if(t2) ras = raster::t(ras)
  extent(ras) = e
  ras = setZ(ras, z = time, name = "time")
  ras
}

# scale ------------------------------------------------------------------------
# plotting the data
scale <- function(x, to, from){   
  (x - min(x, na.rm = TRUE))/(max(x, na.rm = TRUE)-min(x, na.rm = TRUE)) * (to - from) + from
}

# bloomclim --------------------------------------------------------------------
# Input is a Rasterstack object with a set Z value as a Datetime variable. 
# Output is a rasterstack with 12 raster layers  one for each month. To be 
# ploted using levelplot or mylevelplot.
# Arguments:
#    x: Raster stack with a Datetime value as the Z value. 
#       Should have at least one value for each month, otherwise an error given
bloomclim = function(x) {
  themonths <- c("January","February", "March", "April", "May","June",  "July",
                 "August", "September", "October", "November","December")
  sdate <- time(x)
  sdate <- anydate(sdate)
  m     <- months(sdate) # getting the months from the Z value
  j    = 0 # start j at 0 to count loops. i is a string
  # placeholder raster
  clim = x[[1]]

  # looping through each month and finding a spacial average using calc()
  for (i in themonths) {
    j = j + 1
    idx = which(m == i)
    z =  x[[idx]]
    z = median(z, na.rm = TRUE)
    clim = c(clim, z)
  }
  # Should add a time stamp and extent to the raster
  clim[[2:13]]
  }

# anomalize --------------------------------------------------------------------
anomalize = function(ras, detrend = FALSE, f = 0.6){
  themonths <- c("January","February", "March", "April", "May","June",  "July",
                 "August", "September", "October", "November","December")

  # find the monthly climotology of the data set 
  ras_clim = bloomclim(ras)
  
  # subtract each month from corresponding daily data set
  ogt = time(ras)
  ogt = anydate(ogt)
  mon_raw  = months(ogt)
  s = dim(ras_clim)
  
  # placeholder raster
  chla  = ras[[1]] 
  j    = 0 # start j at 0 to count loops. i is a string

  for (mon in themonths) {
    j = j + 1
    ind  = which(mon_raw == mon)
    temp = ras[[ind]] - ras_clim[[j]]
    chla = c(chla, temp)
  }
  
  rm(temp)
  chla = chla[[2:nlyr(chla)]]
  idx = order(time(chla))
  chla = chla[[idx]]
  
  if(detrend == TRUE) {
    clim = smooth.time.series(chla, f = f, smooth.data = TRUE)
    chla = chla - clim
    chla = setZ(chla, z = anydate(t), name = "time")
    }

  chla
}

# match ------------------------------------------------------------------------

# match = function(x, y){
#   ty = getZ(y)
#   ty = as.numeric(ty)
#   tx = getZ(x)
#   tx = as.numeric(tx)
#   idx = which(is.element(tx, ty))
#   tx = subset(tx, is.element(tx, ty))
#   x = raster::subset(x, idx)
#   x = setZ(x, z = as.Date(tx), name = "time")
#   x
# }

# example 
# sla = match(x = sla, y = chl)

timesnip = function(x,
                     sdate = "2003-07-01",
                     edate = "2003-10-01") {
  # subset time
  t = time(x)
  ind_t = which(t >= sdate & t <= edate)
  indx_bool = (t >= sdate & t <= edate)
  x = subset(x, ind_t)
  t = subset(t, indx_bool)
  time(x) = t
  x
  }

