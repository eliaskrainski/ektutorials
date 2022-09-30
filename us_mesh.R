### load packages
library(maps)
library(INLA)

if(!any(ls()=='plotting'))
    plotting <- TRUE

### get a US map
map <- map('state', plot=FALSE, fill=TRUE)
mapID <- sapply(strsplit(map$names, ":"), "[", 1L)

### convert to SpatialPolygons 
map.sp <- maptools::map2SpatialPolygons(map, ID=mapID)
proj4string(map.sp) <- "+proj=longlat +datum=WGS84"

### to mollweide projection
map.moll <- spTransform(map.sp, '+proj=moll +units=km') 

### define a boundary from US map
bound <- inla.sp2segment(map.moll)

### build a mesh on the US boundary
mesh <- inla.mesh.2d(
    boundary=bound,
    max.edge=c(25, 200),
    offset=c(25, 300),
    cutoff=25,
    n=25,
    min.angle=25)
mesh$n

if(plotting) {
    
    plot(mesh, asp=1)

}

### project into a grid
(bb <- bbox(map.moll))
(r <- apply(bb, 1, diff))
(grid.dim <- round(r/4))

### create a projector 
system.time(grid.proj <- inla.mesh.projector(
                mesh, xlim=bb[1,], ylim=bb[2,], dims=grid.dim))

### identify pixels outside the map
system.time(id.grid.out <- which(is.na(over(
    SpatialPoints(grid.proj$lattice$loc, map.moll@proj4string),
    map.moll))))

