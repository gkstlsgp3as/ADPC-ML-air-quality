library(ncdf4) # package for netcdf manipulation
#library(ncdf) # package for netcdf manipulation
library(raster) # package for raster manipulation
#library(rgdal) # package for geospatial analysis
library(ggplot2) # package for plotting
library(corrplot)
library(randomForest)
library(mgcv)
library(gstat)
library(reshape)
library(tidyverse)

setwd('C:/Users/com/Documents/ADPC/ML/satellite')

files<-list.files(path="./",pattern=NULL)
files
for (file in files[1:8]){
  s5p <- list.files(path=file)

  ## load rasters
  o3_data <- raster(paste0('./',file,'/',s5p[11]))
  o3_data@crs <- CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs")
  
  so2_data <- raster(paste0('./',file,'/',s5p[12]))
  so2_data@crs <- CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs")
  
  no2_data <- raster(paste0('./',file,'/',s5p[10]))
  no2_data@crs <- CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs")
  
  co_data <- raster(paste0('./',file,'/',s5p[8]))
  co_data@crs <- CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs")
  
  hcho_data <- raster(paste0('./',file,'/',s5p[9]))
  hcho_data@crs <- CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs")
  
  ai_data <- raster(paste0('./',file,'/',s5p[5]))
  ai_data@crs <- CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs")
  
  dem_data <- raster(paste0('./',file,'/',s5p[1]))
  dem_data@crs <- CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs")
  
  ## stack rasters 
  rst = setNames(stack(o3_data, so2_data, no2_data, hcho_data, ai_data, co_data, dem_data),#, lon_data, lat_data),
                 c('o3', 'so2', 'no2', 'hcho', 'ai', 'co', 'dem'))#, 'Lon','Lat'))
  
  ## load ground data
  df_kr <- read.csv(paste0('./',file,'/',s5p[4]), header=T)
  df_kr <- df_kr[,c(3,8,9)]
  head(df_kr)
  names(df_kr) <- c('PM25','Lon','Lat')
  coordinates(df_kr) <- ~Lon+Lat
  
  ## extract values from raster files
  rasValue=raster::extract(rst, df_kr)
  combinePointValue=cbind(df_kr,rasValue)
  
  write.table(combinePointValue,file=paste0(paste0('./',file,'/combined_',s5p[4])), append=FALSE, sep= ",", row.names = FALSE, col.names=TRUE)
}

## extract lat and lon raster files
# a <- o3_data
# dem_crs <- projectRaster(dem_data, crs=crs(o3_data))
# s<-resample(dem_data,o3_data, method='ngb')
# writeRaster(s, './satellite/S5P_mean_mar/dem_new.tif')

## load final data
setwd('C:/Users/com/Documents/ADPC/ML/satellite')
files<-list.files(path="./",pattern=NULL)
files
file <- files[1]
file 
s5p <- list.files(path=file)
s5p
df_kr <- read.csv(paste0('./',file,'/',s5p[4]), header=T)
df_kr <- df_kr[,c(-11)]
head(df_kr)

lat_data <- raster(paste0('./',file,'/',s5p[2]))
lat_data@crs <- CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs")

lon_data <- raster(paste0('./',file,'/',s5p[3]))
lon_data@crs <- CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs")

o3_data <- raster(paste0('./',file,'/',s5p[12]))
o3_data@crs <- CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs")

so2_data <- raster(paste0('./',file,'/',s5p[13]))
so2_data@crs <- CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs")

no2_data <- raster(paste0('./',file,'/',s5p[11]))
no2_data@crs <- CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs")

co_data <- raster(paste0('./',file,'/',s5p[9]))
co_data@crs <- CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs")

hcho_data <- raster(paste0('./',file,'/',s5p[10]))
hcho_data@crs <- CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs")

ai_data <- raster(paste0('./',file,'/',s5p[6]))
ai_data@crs <- CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs")

dem_data <- raster(paste0('./',file,'/',s5p[1]))
dem_data@crs <- CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs")

rst = setNames(stack(o3_data, so2_data, no2_data, hcho_data, ai_data, co_data, dem_data),#, lon_data, lat_data),
               c('o3', 'so2', 'no2', 'hcho', 'ai', 'co', 'dem'))#, 'Lon','Lat'))

