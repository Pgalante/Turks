library(dismo)
library(spocc)
library(rgeos)
library(dplyr)
library(sp)
library(spThin)
library(scrubr)
library(ENMeval)
## Get Species name
species<-"Meleagris gallopavo"
## Get subspecies names
subspecies<-c("silvestris","osceola","intermedia","merriami","mexicana","gallopavo")

## Source function that queries GBIF
TurkGetter<-function(subspecies, species){
  lapply(subspecies, function(x){
    return(spocc::occ(paste(species, x, sep=' '), from = 'gbif', has_coords=T))
  })
}
## Run function
Turks<-TurkGetter(subspecies, species)
## Cleaning up some data
for (i in 1:length(Turks)){
  Turks[[i]]<-unique(na.omit(data.frame(occ2df(Turks[[i]]))))[,1:3]
}
names(Turks)<-subspecies

## Scrub data localities
scrubTurks <- lapply(Turks, function(x) dframe(x %>% coord_unlikely()))

## Thin
thinTurks <- lapply(scrubTurks, function(x) thin(loc.data = data.frame(x), lat.col = "latitude",long.col = "longitude", spec.col = "name",
                                                 thin.par = 10, reps = 1, locs.thinned.list.return = T, write.files = F))

## Get environmental data
#env<-getData(name = 'worldclim', download = T, path = "/home/pjg/worldClim/wc2-5", res = 2.5, var = 'bio', lon = )
env<-stack(list.files(path = "/home/pjg/worldClim/wc2-5", pattern = "\\.tif$", full.names = T))
## Crop backgrounds
# Source Minimum Convex Hull function
mcp <- function (xy) { 
  xy <- as.data.frame(coordinates(xy))
  coords.t <- chull(xy[, 1], xy[, 2])
  xy.bord <- xy[coords.t, ]
  xy.bord <- rbind(xy.bord[nrow(xy.bord), ], xy.bord)
  return(SpatialPolygons(list(Polygons(list(Polygon(as.matrix(xy.bord))), 1))))
}
# Create hulls for each subspecies
hulls <- lapply(thinTurks, function(x) mcp(x))
# Mask env by each of the hulls. List of rasterStacks
turkEnv <- lapply(hulls, function(x) mask(env, x))
# Crop environmental data from global to local
turkCrop <- list()
for (i in 1:length(turkEnv)){
  turkCrop[[i]] <- crop(turkEnv[[i]], extent(hulls[[i]]))
}

## Model
turkRes<-list()
for (i in 1:length(turkCrop)){
  turkRes[[i]]<-ENMevaluate(occ=data.frame(thinTurks[[i]][[1]]), env = turkCrop[[i]], method='block', parallel=T, numCores=4,
              fc=c("L", "LQ", "H"), RMvalues=seq(0.5,4,0.5), rasterPreds=F)@results
}

sum(is.na(extract(x = turkCrop[[1]], y = thinTurks[[1]][[1]])))
plot(turkCrop[[1]][[1]])
points(thinTurks[[1]][[1]])


<- lapply(thinTurks, function(x) ENMevaluate(occ=data.frame(x), env = turkCrop, method='block', parallel=T, numCores=4,
                                            fc=c("L", "LQ", "H"), RMvalues=seq(0.5,4,0.5), rasterPreds=F))

test<-function(x,y){
  Turkresults<-ENMevaluate(occ=data.frame(x), env = y, method='block', parallel=T, numCores=4,
                        fc=c("L", "LQ", "H"), RMvalues=seq(0.5,4,0.5), rasterPreds=F)
  return(Turkresults@results)
}
test(thinTurks[[1]], turkEnv[[1]])

Turkresults<-ENMevaluate(occ=data.frame(thinTurks[[1]]), env = turkEnv[[1]], method='block', parallel=T, numCores=4,
                         fc=c("L", "LQ", "H"), RMvalues=seq(0.5,4,0.5), rasterPreds=F, bg.coords = NULL)
