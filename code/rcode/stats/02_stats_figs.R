source("code//functions.R")
source("code//libraries.R")

# ------------------------------------------------------------------------------
# loads the summary table created in stat_table.R file and makes some fiugures 
# from the information. 

sumtbl = read.csv("data/summary_statistics/full_sum_20230218.csv")

# making a pretty table --------------------------------------------------------

tbl = sumtbl
rownames(tbl) = NULL
#idx = !year(tbl$sdate) == 2008
#tbl = tbl[idx,]
kbl(tbl, 
    booktabs = T, 
    latex_options = c("striped"),
    caption = "Summary statistics of each bloom around St. ALOHA and 30N between 2002 to 2019. Magnitude taken as the maximum CHL $mg/m^3$ value the blooms reached, and area is the maximum area the bloom reached.") %>%
  kable_styling(full_width = F, latex_options = "HOLD_position") %>%
  column_spec(1:2, bold = F, color = "black") %>%
  column_spec(1:7, width = "10em")
# row_spec(c(4, 5, 8, 11, 12, 17), bold = F, color = "black", background = "#add8e6")
# add_header_above(c("Table of (sub)mesoscale contributions: St. ALOHA" = 5), bold=T)

# ------------------------------------------------------------------------------
# Big Bloom Magnitudes
tbl = sumtbl
tbl[y == 2010, ]$mag = 0.142

dt = gsub("-", "", as.character(Sys.Date()))
pdf(paste("figures\\bloom_magnitude_", dt, ".pdf", sep = ""),   # The directory you want to save the file in
    width = 6, # The width of the plot in inches
    height = 4,
    pointsize = 10) # The height of the plot in inches

main = ""
# ylab = TeX(r'(CHL $[mg/mg^3]$)')
ylab = expression(Bloom ~ Magnitude ~ CHL ~ (mg ~ m^{-3}))
xlab = "Year"
# this changes the linewidth outsidet he plot

opar = par(lwd = 1)
par(mar = c(5,5,1,2))
barplot(tbl$mag, 
        names.arg = year(tbl$sdate), 
        main = "",
        ylab = "", 
        xlab = "",
        col = "grey30",
        ylim = c(0, max(tbl$mag) + max(tbl$mag)*0.05),
        axes = FALSE,
        cex.names = 1,
        las = 1,
        cex.axis = 1)
box(which = "plot", lty = "solid", lwd = 1.5, col = "black")
grid(nx = NA, # X-axis divided in two sections
     ny = NULL,
     lty = 2, 
     col = colvect(c("gray69"), alpha = 0.5), lwd = 1)
axis(side = 2,
     las  = 2,
     lwd  = 1,
     cex.axis = 1,
     col = "black")
title(main = main,
      cex.lab = 2,
      line= 0.75,
      adj = 0)
title(ylab = ylab, cex.lab = 1, line = 2.5)
title(xlab = xlab, cex.lab = 1, line = 2.5)
dev.off()

# ------------------------------------------------------------------------------
# Big Bloom Areas
tbl = sumtbl
tbl[y == 2022, ]$area = 682983.1

dt = gsub("-", "", as.character(Sys.Date()))
pdf(paste("figures\\bloom_area_", dt, ".pdf", sep = ""),   # The directory you want to save the file in
    width = 6, # The width of the plot in inches
    height = 4,
    pointsize = 10) # The height of the plot in inches

# this changes the linewidth outsidet he plot
par(mar = c(5,5,1,2))
opar = par(lwd = 1) #2
barplot(tbl$area/100000, 
        axes = FALSE,
        names.arg = year(tbl$sdate), 
        main = "",
        ylab = "", 
        col = "grey30",
        xlab = "",
        ylim = c(0, max(tbl$area/100000) + max(tbl$area/100000)*0.05))
box(which = "plot", lty = "solid", lwd = 1.5, col = "black")
grid(nx = NA, # X-axis divided in two sections
     ny = NULL, # Y-axis divided in three sections
     lty = 2, 
     col = colvect(c("gray69"), alpha = 0.5), lwd = 1)
axis(side = 2,
     las  = 2, 
     lwd  = 1, 
     mgp  = c(1, 0.75, 0), 
     cex.axis = 1,
     col = "black")
# title(main = main,
#       cex.lab = 1,
#       line= 0.75,
#       adj = 0)
title(ylab = expression(Bloom ~ Area ~ (km * 10^{5})), 
      cex.lab = 1, 
      line = 2.5)
title(xlab = "Year", 
      cex.lab = 1, 
      line = 2.5)
dev.off()

# ------------------------------------------------------------------------------

tbl = sumtbl
y = year(tbl$sdate)
tbl[y == 2010, ]$mag = 0.142
tbl[y == 2022, ]$area = 682983.1

dt = gsub("-", "", as.character(Sys.Date()))
pdf(paste("figures\\bloom_magnitude_area_", dt, ".pdf", sep = ""),   # The directory you want to save the file in
    width = 8, # The width of the plot in inches
    height = 8,
    pointsize = 10) # The height of the plot in inches

layout( matrix(c(1,
                 2),
               ncol = 1,
               nrow = 2, 
               byrow = TRUE),
        heights = c(1, 1))

opar = par(lwd = 1)
par(mar = c(3,5,3,2))
barplot(tbl$mag, 
        names.arg = year(tbl$sdate), 
        main = "",
        ylab = "", 
        xlab = "",
        ylim = c(0, max(tbl$mag) + max(tbl$mag)*0.05),
        axes = FALSE,
        cex.names = 1,
        col = "grey30",
        las = 1,
        cex.axis = 1
        #density = 1
)
grid(nx = NA, # X-axis divided in two sections
     ny = NULL,
     lty = 2, 
     col = colvect(c("gray69"), alpha = 0.5), lwd = 1)
axis(side = 2,
     las  = 2,
     lwd  = 1,
     cex.axis = 1,
     col = "black")
box(which = "plot", lty = "solid", lwd = 1.5, col = "black")
title(ylab = expression(Bloom ~ Magnitude ~ CHL ~ (mg ~ m^{-3})), 
      cex.lab = 1, 
      line = 2.5)

# this changes the linewidth outsidet he plot
par(mar = c(5, 5, 0.5, 2))
opar = par(lwd = 1) #2
barplot(tbl$area/100000, 
        axes = FALSE,
        names.arg = year(tbl$sdate), 
        main = "",
        ylab = "", 
        xlab = "",
        col =  "grey30",
        ylim = c(0, max(tbl$area/100000) + max(tbl$area/100000)*0.05))
grid(nx = NA, # X-axis divided in two sections
     ny = NULL, # Y-axis divided in three sections
     lty = 2, 
     col = colvect(c("gray69"), alpha = 0.5), lwd = 1)
axis(side = 2,
     las  = 2, 
     lwd  = 1, 
     mgp  = c(1, 0.75, 0), 
     cex.axis = 1,
     col = "black")
title(ylab = expression(Bloom ~ Area ~ (km ~ x ~10^{5})), 
      cex.lab = 1, 
      line = 2.5)
title(xlab = "Year", 
      cex.lab = 1, 
      line = 2.5)
box(which = "plot", lty = "solid", lwd = 1.5, col = "black")

dev.off()


#-------------------------------------------------------------------------------
# Bloom Centers
dt = gsub("-", "", as.character(Sys.Date()))
pdf(paste("figures\\bloom_centers_", dt, ".pdf", sep = ""),   # The directory you want to save the file in
    width = 8, # The width of the plot in inches
    height = 5.5,
    pointsize = 10) # The height of the plot in inches

tbl = sumtbl
e = extent(-168, -125, 15, 42)
wdmap <- getMap(resolution = "high")
temp = tbl
#idx = tbl$lon > -165
#tbl = tbl[idx,]
#idx = !year(tbl$sdate) == 2015
#tbl = tbl[idx,]

# radius in km
a = tbl$area
r = sqrt(a/pi) 
r = r/110.574 
tbl$rad = r 
tbl$y = as.character(year(tbl$sdate))

idx = order(tbl$area, decreasing = TRUE)
txt = tbl[idx,]
txt = txt[1,]

idx = year(tbl$sdate) == 2018
tbl = tbl[!idx,]

x = e[1:2]
at = round(seq(range(x)[1]+2, range(x)[2]-2, length = 5), 1)

jamie_theme(e[1:2], e[3:4],
            main = "", 
            ylab = expression(paste("Latitude (",degree,"N)")), 
            xlab = expression(paste("Longitude (",degree,"W)")),
            dt = FALSE,
            xaxes = TRUE,
            yaxes = TRUE,
            mar = c(5,5,2,3)
)
axis(1, at = at, tck = 1, lty = 2, col = "grey", labels = NA)
points(tbl$lon, tbl$lat, pch = 4)
for(i in 1:nrow(tbl)) draw.circle(tbl$lon[i], 
                                  tbl$lat[i], 
                                  tbl$rad[i],
                                  border = colvect("black", alpha = 0.9),
                                  lwd = 1,
                                  nv  = 500,
                                  lty = 2)
points(txt$lon, txt$lat, pch = 10, col = "firebrick")
draw.circle(txt$lon[1], 
            txt$lat[1], 
            txt$rad[1],
            border = "firebrick",
            lwd = 1,
            nv  = 500,
            lty = 2)
plot(wdmap, 
     xlim = e[1:2], 
     ylim = e[3:4], 
     asp = 1, 
     bg = "black", 
     border = "black", 
     col = "black", 
     #wrap = c(-180, 180), 
     add = TRUE)
#text(txt$lon, txt$lat, label = txt$y, adj = 1.2)
box(which = "plot", lty = "solid", lwd = 1.5, col = "black")
legend(-138, 21, 
       legend = c("bloom centers", "2018 bloom"), 
       col = c(colvect("black", alpha = 0.9), "firebrick"),
       pch = c(4, 10))
dev.off()

# ------------------------------------------------------------------------------
# Bloom Duration
tbl = sumtbl

dt = gsub("-", "", as.character(Sys.Date()))
pdf(paste("figures\\bloom_duration_", dt, ".pdf", sep = ""),   # The directory you want to save the file in
    width = 6, # The width of the plot in inches
    height = 6,
    pointsize = 12) # The height of the plot in inches

start = yday(tbl$sdate)
end   = yday(tbl$edate)
# end[end < 100] = 365
# end[length(end)] = end[length(end)] + 60

durations = data.frame(id = as.character(year(tbl$sdate)), start, end)
durations$middle = apply(durations[, c("start", "end")], 1, FUN = mean, na.rm = TRUE)
durations$year = year(tbl$sdate)

# seperating years with two blooms 
# idx_dup   = duplicated(durations$year)
# twobloom  = durations[idx_dup,]
# durations = durations[!idx_dup,]
# 
# # doing it again for three bloom years
# idx_dup    = duplicated(twobloom$year)
# threebloom = twobloom[idx_dup,]
# twobloom   = twobloom[!idx_dup,]

ggplot(durations) +
  geom_boxplot(
    fill = "grey30",
    colour = "black",
    stat = "identity",
    aes(x = id, 
        lower  = start, 
        middle = start, 
        upper  = end, 
        ymin   = start, 
        ymax   = end
    )) +
  # geom_boxplot(data = twobloom,
  #              fill = "grey100",
  #              colour = colvect("grey22", alpha = 0.9),
  #              size = 0.75,
  #              stat = "identity",
  #              aes(x = id, 
  #                  lower  = start, 
  #                  middle = middle, 
  #                  upper  = end, 
  #                  ymin   = start, 
  #                  ymax   = end
  #              )) +
  # geom_boxplot(data = threebloom,
  #              fill = "grey100",
  #              colour = colvect("grey22", alpha = 0.9),
  #              size = 0.75,
  #              stat = "identity",
  #              aes(x = id, 
  #                  lower  = start, 
  #                  #middle = middle, 
  #                  upper  = end, 
  #                  ymin   = start, 
  #                  ymax   = end
  #              )) +
  xlab("Bloom year") +
  ylab("Day of year") +
  labs(title = "") + 
  theme_bw() + 
  theme(
    panel.border = element_rect(colour = "black", 
                                fill = NA,  
                                size = 1),
    axis.title.x = element_text(size = rel(1)),
    axis.title.y = element_text(size = rel(1)),
    axis.text.x = element_text(size = rel(1)),
    axis.text.y = element_text(size = rel(1)),
    plot.title = element_text(size  = rel(1))
  ) +
  coord_flip()
dev.off()

# ------------------------------------------------------------------------------













