library(raster) # package for raster manipulation
library(ggplot2) # package for plotting
library(corrplot)
library(randomForest)
library(mgcv)
library(gstat)
library(reshape)
library(sp)
library(maptools)
library(rgdal)
library(sf)
library(dplyr)
library(caret)
library(dismo)

setwd('C:/Users/com/Documents/ADPC/ML')

head(df_kr)

## ----------------------------------------------------------------------------- ##
## Prepare Dataset  ----
## ----------------------------------------------------------------------------- ##
df_kr <- na.omit(df_kr)
row.names(df_kr) <- NULL

proj4Str <- "+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs"
statPoints <- SpatialPointsDataFrame(coords      = df_kr[,c('Lon','Lat')], 
                                     data        = df_kr,
                                     proj4string = CRS(proj4Str))


## ----------------------------------------------------------------------------- ##
## Calculate Correlation among Variables and P-values  ----
## ----------------------------------------------------------------------------- ##
corMat <- cor(df_kr[complete.cases(df_kr[,c(1:10)]),c(1:10)], method="pearson")
plot.new()
corrplot.mixed(corMat, number.cex=0.7, tl.cex = 0.7, tl.col = "black", 
               outline=FALSE, mar=c(0,0,2,20), upper="square", bg=NA)

corname <- list(colnames(df_kr),colnames(df_kr))
cortest <- matrix(data = NA, nrow = 10, ncol =10, byrow = F, dimnames = corname)
for (i in 1:ncol(df_kr[complete.cases(df_kr[,c(1:10)]),c(1:10)])){
  for (j in 1:ncol(df_kr[complete.cases(df_kr[,c(1:10)]),c(1:10)])){
    cortest[i,j] <- round(cor.test(as.numeric(df_kr[,i]),as.numeric(df_kr[,j]), method="pearson", use="complete.obs")$p.value,4)
    #print(i,j,cor.test(as.numeric(na.omit(df_kr[,i])),as.numeric(na.omit(df_kr[,j])), method="pearson", use="complete.obs")$p.value)
  }
}
cortest

## ----------------------------------------------------------------------------- ##
## Define Grid  ----
## ----------------------------------------------------------------------------- ##
bb <- c(
  "xmin" = 97,
  "ymin" = 5.5,
  "xmax" = 106,
  "ymax" = 20.5
)

grid <- expand.grid(
  X = seq(from = bb["xmin"], to = bb["xmax"], by = 0.01),
  Y = seq(from = bb["ymin"], to = bb["ymax"], by = 0.01) # 20 m resolution
)

x<-grid[,1]
y<-grid[,2]
grid<-data.frame(x,y)
coordinates(grid) <- ~ x+y
gridded(grid) <- TRUE

## ----------------------------------------------------------------------------- ##
## Split Data into Test and Train Set  ----
## ----------------------------------------------------------------------------- ##
kfoldSplit <- function(x, k=10, train=TRUE){
  x <- sample(x, size = length(x), replace = FALSE)
  out <- suppressWarnings(split(x, factor(1:k)))
  if(train) out <- lapply(out, FUN = function(x, len) (1:len)[-x], len=length(unlist(out)))
  return(out)
}

df_kr.pm25 <- subset(df_kr, !is.na(PM25), c(PM25, o3, so2,no2,hcho,ai,co, dem, Lon, Lat))
row.names(df_kr.pm25) <- NULL
statPoints <- subset(df_kr, !is.na(PM25), c(PM25, o3, so2,no2,hcho,ai,co, dem, Lon, Lat))

statPoints <- SpatialPointsDataFrame(coords      = statPoints[,c("Lon","Lat")], 
                                     data        = statPoints,
                                     proj4string = CRS(proj4Str))

resid.RF <- function(x) return(x$y - x$predicted)

set.seed(5132015)

k <- 10

kfolds <- kfoldSplit(1:nrow(df_kr.pm25), k = 10, train = TRUE)

evalData <- matrix(NA, nrow=k, ncol=12, 
                   dimnames = list(1:k, c("IDW_RMSE", "IDW_SI","IDW_R2",
                                          "OK_RMSE", "OK_SI", "OK_R2",
                                          "RF_RMSE", "RF_SI","RF_R2",
                                          "RK_RMSE","RK_SI", "RK_R2")))

avg <- mean(df_kr.pm25$PM25)
avg

## ----------------------------------------------------------------------------- ##
## Evaluation with 10 Test and Train Sets  ----
## ----------------------------------------------------------------------------- ##

for(i in 1:k){
  cat("K-fold...",i,"of",k,"....\n")
  
  # TRAIN indices as integer
  idx <- kfolds[[i]]
  
  # TRAIN indices as a boolean vector
  idxBool <- (1:nrow(df_kr.pm25)) %in% idx
  
  # Observed test data for the target variable
  obs.test <- df_kr.pm25[!idxBool, "PM25"]

  ## ----------------------------------------------------------------------------- ##
  ## Inverse Distance Weighted Method  ----
  ## ----------------------------------------------------------------------------- ##
  
  IDW <-idw(formula= PM25 ~ 1, locations=statPoints[idxBool, ], newdata=statPoints[!idxBool, ], idp=2)
  IDW.output=as.data.frame(IDW)
  names(IDW.output)[1:3]<-c("long","lat","PM25.pred")
  head(IDW.output)
  
  idw.pred.test <- IDW.output$PM25.pred
  evalData[i,"IDW_RMSE"] <- sqrt(mean((idw.pred.test - obs.test)^2))
  evalData[i,"IDW_SI"] <- sqrt(mean((idw.pred.test - obs.test)^2))/avg*100
  evalData[i,"IDW_R2"] <- cor(idw.pred.test,obs.test)^2 # R squared
    
  ## ----------------------------------------------------------------------------- ##
  ## Ordinary Kriging ----
  ## ----------------------------------------------------------------------------- ##
  
  # Make variogram
  formMod <- PM25 ~1
  mod <- vgm(nugget=0,psill=150,range=1000,model="Sph")
  plot(variog <- variogram(formMod, statPoints[idxBool, ],cutoff=1100,width=1100/11))
  #plot(variogram(PM25 ~1, data=df_kr.pm25,cutoff=10, width=10/11))
  
  # Variogram fitting by Ordinary Least Sqaure
  variogFitOLS<-fit.variogram(variog, model = mod,  fit.method = 6)
  plot(variog, variogFitOLS, main="OLS Model")
  variogFitOLS
  
  # kriging predictions
  OK <- krige(formula = formMod ,
              locations = statPoints[idxBool, ], 
              model = variogFitOLS,
              newdata = statPoints[!idxBool, ],
              debug.level = 0)
  
  ok.pred.test <- OK@data$var1.pred
  evalData[i,"OK_RMSE"] <- sqrt(mean((ok.pred.test - obs.test)^2))
  evalData[i,"OK_SI"] <- sqrt(mean((ok.pred.test - obs.test)^2))/avg*100
  evalData[i,"OK_R2"] <- cor(ok.pred.test,obs.test)^2 # R squared
  
  ## ----------------------------------------------------------------------------- ##
  ## Random Forest ----
  ## ----------------------------------------------------------------------------- ##
  
  RF <- randomForest(y = df_kr.pm25[idx, "PM25"], 
                     x = df_kr.pm25[idx, c(-1)],
                     ntree = 500,
                     mtry = 2)
  
  rf.pred.test <- predict(RF, newdata = df_kr.pm25[-idx,], type="response")
  evalData[i,"RF_RMSE"] <- sqrt(mean((rf.pred.test - obs.test)^2))
  evalData[i,"RF_SI"] <- sqrt(mean((rf.pred.test - obs.test)^2))/avg*100
  evalData[i,"RF_R2"] <- cor(rf.pred.test,obs.test)^2 # R squared
  
  # Ordinary Kriging of Random Forest residuals
  #
  statPointsTMP <- statPoints[idxBool, ]
  statPointsTMP@data <- cbind(statPointsTMP@data, residRF = resid.RF(RF))
  
  formMod <- residRF ~ 1
  mod <- vgm(nugget=0,psill=15,range=200,model="Sph")
  plot(variog <- variogram(formMod, statPointsTMP,cutoff=1100,width=1100/10))
  
  # Variogram fitting by Ordinary Least Sqaure
  variogFitOLS<-fit.variogram(variog, model = mod,  fit.method = 6)
  plot(variog, variogFitOLS, main="OLS Model")
  variogFitOLS
  
  # kriging predictions
  RF.OK <- krige(formula = formMod ,
                 locations = statPointsTMP, 
                 model = variogFitOLS,
                 newdata = statPoints[!idxBool, ],
                 debug.level = 0)
  
  rf.ok.pred.test <- rf.pred.test + RF.OK@data$var1.pred
  evalData[i,"RK_RMSE"] <- sqrt(mean((rf.ok.pred.test - obs.test)^2))
  evalData[i,"RK_SI"] <- sqrt(mean((rf.ok.pred.test - obs.test)^2))/avg*100
  evalData[i,"RK_R2"] <- cor(rf.ok.pred.test,obs.test)^2 # R squared
  
}
round(apply(evalData,2,FUN = function(x,...) c(mean(x,...),sd(x,...))),3)

## ----------------------------------------------------------------------------- ##
## Producing PM2.5 Map ----
## ----------------------------------------------------------------------------- ##
## 1) idw
coordinates(df_kr.pm25)<-~Lon+Lat 
coordinates(df_kr)<-~Lon+Lat 
IDW <-idw(formula= PM25 ~ 1, locations=df_kr.pm25, newdata=grid, idp=2)
IDW.output=as.data.frame(IDW)
names(IDW.output)[1:3]<-c("long","lat","PM25.pred")

# ggplot(data=IDW.output,aes(x=long,y=lat))+
#   geom_tile(data=IDW.output,aes(fill=PM25.pred))+
#   scale_fill_gradient(low="#0000FF", high="#FF0000")+coord_equal()

int_bk_gs<-data.frame('x'=IDW.output$long, 'y'=IDW.output$lat, 'z'=IDW.output$PM25.pred)
rst<-rasterFromXYZ(int_bk_gs, res=c(NA,NA), crs='+proj=longlat +datum=WGS84', digits=5)
rf <-writeRaster(rst,'./interp/harp/all_Feb_idw_samp.tif')#,overwrite=TRUE)

## 2) ordinary kriging
# Make variogram
formMod <- PM25 ~1
mod <- vgm(psill=20,range=100,model="Sph")
plot(variog <- variogram(formMod, statPoints,cutoff=1100,width=1100/11))

# Variogram fitting by Ordinary Least Sqaure
variogFitOLS<-fit.variogram(variog, model = mod,  fit.method = 6)
plot(variog, variogFitOLS, main="OLS Model")
variogFitOLS

OK <- gstat::krige(formula = PM25 ~ 1 ,
            locations = df_kr, 
            model = variogFitOLS,
            newdata = grid,
            debug.level = 0)

ok.pred.test <- OK@data$var1.pred
OK.output=as.data.frame(OK)
summary(OK.output)
int_bk<-data.frame('x'=OK.output$x, 'y'=OK.output$y, 'z'=OK.output$var1.pred)
rst<-rasterFromXYZ(int_bk, res=c(NA,NA), crs='+proj=longlat +datum=WGS84', digits=5)
rf <-writeRaster(rst,'./interp/harp/all_Jul_OK_samp2.tif',overwrite=TRUE)

## 3) regression kriging and random forest 
df_kr.pm25 <- as.data.frame(df_kr.pm25)
RF <- randomForest(y = df_kr.pm25[, "PM25"], 
                   x = df_kr.pm25[, c(-1)],
                   ntree = 500,
                   mtry = 2)
importance(RF)
varImpPlot(RF)

rst = setNames(stack(o3_data, so2_data, no2_data, hcho_data, ai_data, co_data, dem_data, lon_data, lat_data),
               c('o3', 'so2', 'no2', 'hcho', 'ai', 'co', 'dem', 'Lon','Lat'))

rstPredRF <- predict(rst, RF, type="response")
rstPixDF <- as(rst[[1]], "SpatialPixelsDataFrame")

# Create a temporary SpatialPointsDF object to store rf residuals
statPointsTMP <- statPoints
crs(statPointsTMP) <- crs(rstPixDF)

statPointsTMP@data <- cbind(statPointsTMP@data, residRF = resid.RF(RF))

# Define the kriging parameters and fit the variogram using OLS
formMod <- residRF ~ 1
mod <- vgm(nugget=40,psill=80,range=800,model="Sph")
plot(variog <- variogram(formMod, statPointsTMP,cutoff=1100,width=1100/10))
variogFitOLS<-fit.variogram(variog, model = mod,  fit.method = 6)
variogFitOLS #check variogram parameters

# Plot the variogram
plot(variog, variogFitOLS, main="Semi-variogram of RF residuals")

residKrigMap <- krige(formula = formMod ,
                      locations = statPointsTMP, 
                      model = variogFitOLS,
                      newdata = rstPixDF)

residKrigRstLayer <- as(residKrigMap, "RasterLayer")

gamKrigMap <- rstPredRF #+ residKrigRstLayer

# plot.new()
# plot(gamKrigMap, main="Annual average air temperature\n(RF regression-kriging)",
#      xlab="Longitude", ylab="Latitude", cex.main=0.8, cex.axis=0.7, cex=0.8,
#      xlim=c(100.048,101.042), ylim=c(13.455,14.292))

rf <-writeRaster(gamKrigMap,'./interp/harp/all_Jul_RF_samp2.tif',overwrite=TRUE)

