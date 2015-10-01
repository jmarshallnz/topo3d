library(maptools)

contour <- readShapeSpatial("shp/contour")

breaks = quantile(contour$elevation,seq(0,1,length=255))
cols <- gray(seq(0,1,length=255))
quant = findInterval(contour$elevation, breaks)

cols = (contour$elevation - 240)/20 * 2
png(file="blah.png", width=2000, height=4000, antialias = "none")
plot(contour, col = rgb(red=0, green=cols, blue=0, max=255), lwd=3)
dev.off()

library(tripack)

# we're going to triangulate the areas between contours.
# get all the coordinates on the contours.
# then get their associated elevation.
# and unique-ise them (though there's probably not many duplicates)
# then do a triangulation with tri.mesh.

con_geom = geometry(contour)
xye = NULL
for (i in 1:length(con_geom)) {
  line = unlist(coordinates(con_geom[i]), recursive=FALSE)
  elev = contour$elevation[i]
  for (j in 1:length(line)) {
    xye = rbind(xye, cbind(line[[j]], elev))
  }
  cat("done contour", i, "of", length(con_geom), "\n")
}
xyu = unique(xye)

# triangulate...
library(tripack)
trimesh = tri.mesh(xyu[,1], xyu[,2])
triang  = triangles(trimesh)
# generate a set of x,y,z points based on these triangles in sequence...
trianp = cbind(xyu[triang[,1],1:3],xyu[triang[,2],1:3],xyu[triang[,3],1:3])
triamp = t(matrix(t(trianp), nrow=3))
dim(triamp)
head(triamp,1000)

triamp <- read.csv("tri/triangles.csv")
library(rgl)
triangles3d(triamp[,1], triamp[,2], triamp[,3], col="green")

library(proj4)
src.proj = '+proj=tmerc +lat_0=0 +lon_0=173 +k=0.9996 +x_0=1600000 +y_0=10000000 +ellps=WGS84'
dst.proj = '+proj=longlat  +datum=WGS84 +no_defs'
latlong = project(triamp[,1:2], src.proj, inverse=TRUE, ellps.default = NA)
xlim = range(latlong$x)
ylim = range(latlong$y)

library(RgoogleMaps)
map = GetMap.bbox(xlim, ylim, maptype="satellite", verbose=TRUE, SCALE=4)

# grab the subset of the map we need to map to our coords
map_size = map$size
map_ll   = map$BBOX$ll
map_ur   = map$BBOX$ur

x1 = round((xlim[1] - map_ll[2]) / (map_ur[2] - map_ll[2]) * map_size[1])
x2 = round((xlim[2] - map_ll[2]) / (map_ur[2] - map_ll[2]) * map_size[1])

y1 = round((ylim[1] - map_ll[1]) / (map_ur[1] - map_ll[1]) * map_size[2])
y2 = round((ylim[2] - map_ll[1]) / (map_ur[1] - map_ll[1]) * map_size[2])

# hmm, maybe we can alter the size...
map$myTile = map$myTile[x1:x2, y1:y2]
class(map$myTile) = 'nativeRaster'
attr(map$myTile, "dim") = c(y2-y1+1, x2-x1+1)
attr(map$myTile, "channels") = 4

#' TODO: 
#' 1. Get better quality image from google maps
#' 2. ...
#' 3. Profit!

# try splitting the image up and get a composite one...
nsplit = 3
size = 640
xl = seq(xlim[1], xlim[2], length = nsplit+1)
yl = seq(ylim[1], ylim[2], length = nsplit+1)
v <- list()
for (y in 1:1) {
  # download the maps
  maps <- list()
  for (x in 1:nsplit) {
    maps[[x]] = GetMap.bbox(xl[x+0:1], yl[y+0:1], maptype="satellite")$myTile #destfile=paste0("map_", x, "_", y, ".png"))$myTile
  }
  # process the maps
  for (x in 2:nsplit) {
    nr = dim(maps[[x-1]])[1]
    nc = dim(maps[[x-1]])[2]
    size = dim(maps[[x]])[1]
    # TODO: Need to optimise this
    v = numeric(size/2)
    for (i in 0:(size/2)) {
      v[i] = cor(maps[[x-1]][1,nc+1:100-(i+100)],maps[[x]][1,1:100])
    }
    lag = which.max(v)+100
    # now take 1:(nc - lag - 100) from the first
    # and the rest overlaps by lag+100 I guess?
    ccff = ccf(maps[[x-1]][,1], maps[[x]][,1], lag.max = nr, plot=FALSE)
    lag  = ccff$lag[which.max(ccff$acf)]
    new_image = matrix(integer(0), nrow=size+lag, ncol=nr)
    new_image[1:lag,] = maps[[x-1]][1:lag,]
    new_image[size+1:lag,] = maps[[x]][1:lag + size - lag,]
    new_image[(lag+1):size,] = blend(maps[[x-1]][(lag+1):size,], maps[[x]][1:(size-lag),], 1)
    # back to nativeRaster
    class(new_image) = 'nativeRaster'
    attr(new_image, "dim") = c(nc, lag+size)
    attr(new_image, "channels") = 4
    maps[[x]] <- new_image
  }
  v[[y]] <- maps[[nsplit]]
}

m12 = GetMap.bbox(xl[2:3], yl[1:2], maptype="satellite")
m21 = GetMap.bbox(xl[1:2], yl[2:3], maptype="satellite")
m22 = GetMap.bbox(xl[2:3], yl[2:3], maptype="satellite")
PlotOnStaticMap(m11)
PlotOnStaticMap(m12)
PlotOnStaticMap(m21)
PlotOnStaticMap(m22)

# figure out the overlap
ccff = ccf(m11$myTile[,1], m12$myTile[,1], lag.max = 300, plot=FALSE)
lag  = ccff$lag[which.max(ccff$acf)]
cbind(m12$myTile[1:100,1], m11$myTile[1:100+lag,1])

# now merge the two images
# 1+lag maps to 1 on second image
# so we have an image of width
new_image = matrix(integer(0), nrow=lag + m11$size[1], ncol=m11$size[2])
new_image[1:lag,] = m11$myTile[1:lag,]
new_image[m11$size[2]+1:lag,] = m12$myTile[1:lag + m11$size[1] - lag,]
# the overlap is then a smooth transition...

blend = function(A, B, margin = 2) {

  # convert A + B to argb...
  a1 = bitwAnd(bitwShiftR(A, 24),255)
  r1 = bitwAnd(bitwShiftR(A, 16),255)
  g1 = bitwAnd(bitwShiftR(A, 8),255)
  b1 = bitwAnd(A,255)

  a2 = bitwAnd(bitwShiftR(B, 24),255)
  r2 = bitwAnd(bitwShiftR(B, 16),255)
  g2 = bitwAnd(bitwShiftR(B, 8),255)
  b2 = bitwAnd(B,255)

  # compute the blend on those channels individually
  size = dim(A)[margin]
  scale = matrix((1:size-0.5)/size, nrow=nrow(A), ncol=ncol(A), byrow = margin == 2)
  a = bitwAnd(as.integer(a2*scale + a1*(1-scale)),255)
  r = bitwAnd(as.integer(r2*scale + r1*(1-scale)),255)
  g = bitwAnd(as.integer(g2*scale + g1*(1-scale)),255)
  b = bitwAnd(as.integer(b2*scale + b1*(1-scale)),255)
  
  # bitshift everything back again (this would be way easier in c++)
  bitwShiftL(a,24) + bitwShiftL(r,16) + bitwShiftL(g,8) + b

#  c12 = unlist(lapply(strsplit(sprintf("%x", m12$myTile[,1]), split=""), function(x) { paste0("#",x[7],x[8],x[5],x[6],x[3],x[4],x[1],x[2]) }))
#  c11 = unlist(lapply(strsplit(sprintf("%x", m11$myTile[,1]), split=""), function(x) { paste0("#",x[7],x[8],x[5],x[6],x[3],x[4],x[1],x[2]) }))
  
#  as.integer((A+B)/2)
}

test_nothing = matrix(0, m11$size[2]-lag, m11$size[1])
#new_image[(lag+1):m11$size[2],] = blend(m11$myTile[(lag+1):m11$size[2],], test_nothing, 2)

test = m11
test$myTile = new_image
PlotOnStaticMap(test)

c12 = sprintf("#%x", bitwAnd(bitwShift(m12$myTile[,1],2), 255))

image(cols=cols)
image(matrix(1, 640,640), col=cols)

x1 = (xlim[1] - map_ll[2]) / 
x2 = (xlim[2] - map_ll[2]) / (map_ur[2] - map_ll[2])
y1 = (ylim[1] - map_ll[1]) / (map_ur[1] - map_ll[1])
y2 = (ylim[2] - map_ll[1]) / (map_ur[1] - map_ll[1])

texcoords.x = (latlong$x - map_ll[2]) / (map_ur[2] - map_ll[2])
texcoords.y = (latlong$y - map_ll[1]) / (map_ur[1] - map_ll[1])

triamp = read.csv("tri/triangles.csv")
library(rgl)
triangles3d(triamp[,1], triamp[,2], triamp[,3], lwd=0.01, col="white", texture="map.png", textype="rgb",
            texcoords=triamp[,4:5], shininess = 20)
