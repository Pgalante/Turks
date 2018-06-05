library(spocc)
library(ENMeval)
library(dplyr)
library(rgeos)
library(sp)
library(spThin)
library(scrubr)


library(dismo)

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
# Scrub data localities
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
  return(gBuffer(SpatialPolygons(list(Polygons(list(Polygon(as.matrix(xy.bord))), 1))), width = 1))
}
# Create hulls for each subspecies
hulls <- lapply(thinTurks, function(x) mcp(x))
# Mask env by each of the hulls. List of rasterStacks
turkEnv <- lapply(hulls, function(x) mask(env, x))
# Crop environmental data from global to local
turkCrop<-lapply(1:length(turkEnv), function(x) crop(turkEnv[[x]], hulls[[x]]))

## Model
turkRes<-list()
for (i in 1:length(turkCrop)){
  turkRes[[i]]<-ENMevaluate(occ=data.frame(thinTurks[[i]][[1]]), env = turkCrop[[i]], method='block', parallel=T, numCores=4,
              fc=c("L", "LQ", "H"), RMvalues=seq(0.5,4,0.5), rasterPreds=F)@results
}

resTable<-turkRes[[1]]

## Select optimal model and create logistic predictions
optimize<- function(mod.Table){
  no.zero.param <- mod.Table[mod.Table$nparam != 0,]
  bestOR = no.zero.param[order(no.zero.param[,'Mean.ORmin']),]
  bestAUC = bestOR[order(bestOR[,'Mean.AUC'], decreasing=TRUE),]
  opt.mod = bestAUC[1,]
  return(opt.mod)
}

mod.Table<-turkRes[[6]]

m<-list()
pred<-list()
opt.mod<-list()
for (i in 1:length(turkRes)){
  opt.mod[[i]]<-optimize(turkRes[[i]])
  b.m<-opt.mod[[i]]$rm
  beta.mulr<- paste('betamultiplier=',b.m,sep='')
  false.args<-c('noautofeature','noproduct','nothreshold','noquadratic','nohinge','nolinear')
  feat<-strsplit(as.vector(opt.mod[[i]][,2]), ',')[[1]]
  if (feat == 'L'){feats = 'linear'} else if (feat == 'LQ'){feats = c('quadratic', 'linear')} else if (feat == 'H'){feats = 'hinge'} else if (feat == 'P'){feats = 'product'} else if (feat == 'T'){feats = 'threshold'} else if (feat == 'LQH'){feats = c('linear', 'quadratic', 'hinge')} else if (feat == 'LQHP'){feats = c('linear', 'quadratic', 'hinge', 'product')} else if (feat == 'LQHPT'){feats = c('linear', 'quadratic', 'hinge', 'product', 'threshold')}
  for (j in 1:length(feats)){false.args[which(sub('no','',false.args)==feats[j])] = feats[j]}
  m[[i]] <-maxent(turkCrop[[i]], thinTurks[[i]][[1]], args=c(false.args, beta.mulr, 'noremoveduplicates', 'noaddsamplestobackground'))
  pred[[i]]  <- predict(object= m[[i]], x=turkCrop[[i]], na.rm=TRUE, format='GTiff',overwrite=TRUE, progress='text',args='logistic')
}

## Threshold predictions to 10% OR
modThreshs <- lapply(opt.mod, function(x) return(x$Mean.OR10))
for (i in 1:length(modThreshs)){
  pred[[i]][pred[[i]] < modThreshs[[i]]] <- NA
}

## Compare niches

