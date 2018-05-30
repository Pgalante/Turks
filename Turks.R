library(dismo)
library(spocc)
library(rgeos)
library(dplyr)
library(sp)
library(spThin)
library(scrubr)
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
thinTurks <- lapply(scrubTurks, function(x) thin(loc.data = data.frame(scrubTurks[[1]]), lat.col = "latitude",long.col = "longitude", spec.col = "name",
                                                 thin.par = 10, reps = 1, locs.thinned.list.return = T, write.files = F))

## Model