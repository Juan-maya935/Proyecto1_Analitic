
#Install libraries
install.packages("terra")

#Libraries
library(terra)

#Work directory



################################################################################
# Introduction of spatial data
################################################################################

name <- LETTERS[1:10]
longitude <- c(-116.7, -120.4, -116.7, -113.5, -115.5,
               -120.8, -119.5, -113.7, -113.7, -110.7)
latitude <- c(45.3, 42.6, 38.9, 42.1, 35.7, 38.9,
              36.2, 39, 41.6, 36.9)
stations <- cbind(longitude, latitude)

# Simulated rainfall data
set.seed(0)
precip <- round((runif(length(latitude))*10)^3)


psize <- 1 + precip/500
plot(stations, cex=psize, pch=20, col='red', main='Precipitation')
# add names to plot
text(stations, name, pos=4)
# add a legend
breaks <- c(100, 250, 500, 1000)
legend.psize <- 1+breaks/500
legend("topright", legend=breaks, pch=20, pt.cex=legend.psize, col='red', bg='gray')


lon <- c(-116.8, -114.2, -112.9, -111.9, -114.2, -115.4, -117.7)
lat <- c(41.3, 42.9, 42.4, 39.8, 37.6, 38.3, 37.6)
x <- cbind(lon, lat)
plot(stations, main='Precipitation')
polygon(x, col='blue', border='light blue')
lines(stations, lwd=3, col='red')
points(x, cex=2, pch=20)
points(stations, cex=psize, pch=20, col='red', main='Precipitation')


wst <- data.frame(longitude, latitude, name, precip)
  wst


################################################################################
# Vector data
################################################################################

#Points
longitude <- c(-116.7, -120.4, -116.7, -113.5, -115.5, -120.8, -119.5, -113.7, -113.7, -110.7)
latitude <- c(45.3, 42.6, 38.9, 42.1, 35.7, 38.9, 36.2, 39, 41.6, 36.9)
lonlat <- cbind(longitude, latitude)

pts <- vect(lonlat)
class (pts)
pts

crdref <- "+proj=longlat +datum=WGS84"
pts <- vect(lonlat, crs=crdref)
pts
crs(pts)


precipvalue <- runif(nrow(lonlat), min=0, max=100)
df          <- data.frame(ID=1:nrow(lonlat), precip=precipvalue)


ptv <- vect(lonlat, atts=df, crs=crdref)
ptv

#Lines and polygons
lon <- c(-116.8, -114.2, -112.9, -111.9, -114.2, -115.4, -117.7)
lat <- c(41.3, 42.9, 42.4, 39.8, 37.6, 38.3, 37.6)
lonlat <- cbind(id=1, part=1, lon, lat)
lonlat

lns <- vect(lonlat, type="lines", crs=crdref)
lns

pols <- vect(lonlat, type="polygons", crs=crdref)
pols

plot(pols, las=1)
plot(pols, border='blue', col='yellow', lwd=3, add=TRUE)
points(pts, col='red', pch=20, cex=3)


################################################################################
# Raster data
################################################################################

r <- rast(ncol=10, nrow=10, xmin=-150, xmax=-80, ymin=20, ymax=60)
r

values(r) <- runif(ncell(r))
r

values(r) <- 1:ncell(r)
r


plot(r)
# add polygon and points
lon <- c(-116.8, -114.2, -112.9, -111.9, -114.2, -115.4, -117.7)
lat <- c(41.3, 42.9, 42.4, 39.8, 37.6, 38.3, 37.6)
lonlat <- cbind(id=1, part=1, lon, lat)
pts <- vect(lonlat)
pols <- vect(lonlat, type="polygons", crs="+proj=longlat +datum=WGS84")
points(pts, col="red", pch=20, cex=3)
lines(pols, col="blue", lwd=2)


r2 <- r * r
r3  <- sqrt(r)
s <- c(r, r2, r3)
s
plot(s)

################################################################################
#Handle spatial data
################################################################################

library(terra)
filename <- system.file("ex/lux.shp", package="terra")
basename(filename)

s <- vect(filename)
s

#Writing
outfile <- "shp_test.shp"
writeVector(s, outfile, overwrite=TRUE)

#Reading raster data
f <- system.file("ex/logo.tif", package="terra")
basename(f)

r <- rast(f)
r

r2 <- r[[2]]
r2

#Writing raster data
x <- writeRaster(r, "test_output.tif", overwrite=TRUE)
x

plot(x)

###############################################################################
# Coordinate Reference Systems (CRS)
###############################################################################

library(terra)

f <- system.file("ex/lux.shp", package="terra")
p <- vect(f)
p
crs(p)


###############################################################################
#Vector manipulation
###############################################################################

f <- system.file("ex/lux.shp", package="terra")
p <- vect(f)
p


plot(p)
plot(p, "NAME_2")
plot(p, "POP")
plot(p, "AREA")


#Convert to data.frame
d <- as.data.frame(p)
d

#Extrac geometric
g <- geom(p)
g

#SpatVector with "NAME_2"
p$NAME_2
p[, "NAME_2"]

#Add new feature
p$counts <- rpois(nrow(p), 3)
as.data.frame(p)

#Merge
dfr <- data.frame(District=p$NAME_1, Canton=p$NAME_2, Value=round(runif(length(p), 100, 1000)))
dfr
dfr <- dfr[order(dfr$Canton), ]
dfr

pm <- merge(p, dfr, by.x=c('NAME_1', 'NAME_2'), by.y=c('District', 'Canton'))
head(pm)

#Select
g <- p[which(p$NAME_1 == 'Grevenmacher'),]
head(g)


###############################################################################
#Raster manipulation
###############################################################################

x <- rast()
x

plot(x)


x <- rast(ncol=36, nrow=18, xmin=-1000, xmax=1000, ymin=-100, ymax=900)

res(x) #resolution

ncol(x)
ncol(x) <- 18
ncol(x)
res(x)

#Add coordinates
crs(x) <- "+proj=utm +zone=48 +datum=WGS84"
x
plot(x)


r <- rast(ncol=10, nrow=10)
ncell(r)

#add values
values(r) <- 1:ncell(r)

plot(r)


set.seed(123)
values(r) <- runif(ncell(r))
plot(r)

#show values
values(r)[1:10]



#Load file
filename <- system.file("ex/elev.tif", package="terra")
basename(filename)


r <- rast(filename)

plot(r)
values(r)

###############################################################################
# Temporal autocorrelation
###############################################################################

set.seed(123)
d <- sample(100, 10)
d


a <- d[-length(d)]
b <- d[-1]
plot(a, b, xlab='t', ylab='t-1')

cor(a, b)


d <- sort(d)
d
a <- d[-length(d)]
b <- d[-1]
plot(a, b, xlab='t', ylab='t-1')

cor(a, b)

acf(d)


###############################################################################
# Spatial autocorrelation
###############################################################################

library(terra)
p <- vect(system.file("ex/lux.shp", package="terra"))
p <- p[p$NAME_1=="Diekirch", ]
p$value <- c(10, 6, 4, 11, 6)
as.data.frame(p)


par(mai=c(0,0,0,0))
plot(p, col=2:7)
xy <- centroids(p)
points(xy, cex=6, pch=20, col='white')
text(p, 'ID_2', cex=1.5)


#Adjacent polygons
w <- adjacent(p, symmetrical=TRUE)
class(w)


plot(p, col='gray', border='blue', lwd=2)
p1 <- xy[w[,1], ]
p2 <- xy[w[,2], ]
lines(p1, p2, col='red', lwd=2)

