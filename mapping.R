library(maptools)

contour <- readShapeSpatial("contour")

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
trimesh = tri.mesh(xyu[,1], xyu[,2])
triang  = triangles(trimesh)
# generate a set of x,y,z points based on these triangles in sequence...
trianp = cbind(xyu[triang[,1],1:3],xyu[triang[,2],1:3],xyu[triang[,3],1:3])
triamp = t(matrix(t(trianp), nrow=3))
dim(triamp)
head(triamp,1000)

library(rgl)
triangles3d(triamp[,1], triamp[,2], triamp[,3], col="green")
