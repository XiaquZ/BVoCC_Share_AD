
##Load packages
library(parallelly)
library(doParallel)
library(dplyr)
library(data.table)
library(terra)
library(tictoc)


## Get current and future climate raster

pre <- rast("I:/InputData/Temp/PresentMicro/NorthEU/ForestMAT_NorthEU_2_50kmBufferRing_2000-2020_EPSG3035_100m.tif"); names(pre) <- 'pre' ## Mean annual temperature (1995) 
fut <- rast("I:/InputData/Temp/FutureMicro/NorthEU/ForestMAT_NorthEU_2_BVoCC50km_2071-2100_SSP370_EPSG3035_100m.tif"); names(fut) <- 'fut' ## Mean annual temperature (2085)

show(pre)
show(fut)
#####Subset a raster smaller than a value.####
#The climTol is +-0.25 degree.
##For Forward Velocity, if pre + 0.25 is smaller than the minimum future temperature, then these type of 
#pre cells will never found their climate analogues within EU, since there is no overlap in their ranges. 

##For backward velocity, if fut-0.25 is still larger than the maximum present temperature, 
#then the fut cells won't be able to find the analog present location.
minmax_pre <- data.frame(minmax(pre))
fut01 <- fut - 0.25 #The minimum climate analog value for the focus location/cell.
fut_avai <- ifel(fut01 > minmax_pre[2,1], NA, fut)
plot(fut)
plot(fut_avai)
plot(pre)
##Create a stack for only pre and fut data, so 2 layers to only run the VoCC 
##algorithm within the
##range of pre map.
# the.stack <-  c(pre,fut) #DON't RUN this in this case.
# plot(the.stack) #DON't RUN this in this case.

####Create initial climate data frame.####

##future microclimate data
climpre <- pre %>%
  values() %>%
  data.frame(cid = 1:ncell(pre)) %>%
  na.omit()
##present microclimate data 
climfut <- fut_avai %>%
  values() %>%
  data.frame(cid = 1:ncell(fut_avai)) %>%
  na.omit()

climpre[,c("x", "y")] <- xyFromCell(pre, climpre$cid) 
climfut[,c("x", "y")] <- xyFromCell(fut_avai, climfut$cid) 
# Add a new column combining x and y in both data frames
climpre$xy_combined <- paste(climpre$x, climpre$y, sep = "_")
climfut$xy_combined <- paste(climfut$x, climfut$y, sep = "_")
clim <- merge(climpre,climfut, by = "xy_combined", all.x = TRUE, all.y = TRUE)
# Remove the xy_combined column (if you don't need it anymore)
clim <- clim[, c("pre", "fut", "cid.x","cid.y", "x.x", "y.x", "x.y", "y.y")]
colnames(clim) <- c("pre", "fut", "cid.pre","cid.fut", "x_pre", "y_pre", "x_fut", "y_fut")

# clim <- the.stack %>%
#   values() %>%
#   data.frame(cid = 1:ncell(the.stack)) %>%
#   na.omit()

#Get the coordination for each pixel.

##Round the pre1 and fut30 to keep 1 decimal.
clim$pre <- round(clim$pre,1) 
clim$fut <- round(clim$fut, 1)
head(clim)
clim <- as.data.table(clim)

#### Make a foreach parallel processing around the code you want to run for each temperature value####

##Some input indices.
n=1 #n=1 #The number of climate variable, in this case, only one climate variable (mean annual temp)
tdiff = 75 #tdiff = 90 #Time interval between present and future.
method = "Single" #Single threshold method, i.e. "climTol" is a constant value.
climTol = 0.25 ##climTol = +- 0.25 degrees C ; total bin diameter = 0.5 degrees C 
geoTol = 150000 #geoTol = 180000 #Searching radius for the analog pixel. (2 km yr-1)
#cost.penalty <- 2 ## For LCP method: Two penalty units per degree C dissimilarity from temperature of interest

##################################################################################
#' Distance-based velocity based on geographically closest climate analogue
#' Initial dVoCC code modified by Aggeliki Doxa for Xiaqu Zhou analysis
##################################################################################

dat <- data.table(clim)
filtered_dat <- dat[!is.na(dat$fut),]

## set things up for parallel processing.
ncores = detectCores()-2
cuts <- cut(1:nrow(filtered_dat), ncores, labels = FALSE) #Skip the rows with NA in  future data when running the foreach.
#cuts <- cut(1:nrow(dat), ncores, labels = FALSE)
cl <- makeClusterPSOCK(ncores, autoStop = TRUE)
registerDoParallel(cl)
print(ncores)

# profvis({
tic("B_VoCC")  
result <- foreach(x = 1:ncores,.packages = c('terra', 'data.table'),
                  .combine = rbind, 
                  .multicombine = TRUE) %dopar% { 
                    
                    a <- x
                    Dat <- filtered_dat[cuts == a,]
                    
                    resu <- data.table(focal = Dat$cid.fut, target = as.integer(NA), geoDis_m = as.double(NA), vel_m_yr = as.double(NA))

                    i <- 0
                    while(i <= nrow(Dat)){
                      #for (i in 1:nrow(Dat)) {
                      i <- i+1
                      
                      #### for each focal cell subset neighbor cells # Added by AD
                      pres_XY <- as.numeric(Dat[i,(ncol(Dat)-1):ncol(Dat)])  # Added by AD
                      dif_XY <- data.table(sweep(dat[,(ncol(dat)-3):(ncol(dat)-2)], 2, pres_XY, "-"))  # Added by AD
                      dif_XY$Eucl_dis <- sqrt((dif_XY$x_pre)^2 + (dif_XY$y_pre)^2) #Added by X.Z
                      # Keep only cells with X,Y differences < than k threshold # Added by AD
                      #hist(dif_XY$x) #check the distribution of differences in X.Added by XZ
                      #hist(dif_XY$y) #check the distribution of differences in Y.Added by XZ
                      upperXY = colnames(dif_XY[,3])  # Added by AD. edited by X.Z
                      l_XY <- lapply(upperXY, function(x) call("<=", call("abs", as.name(x)), 50000))  # Added by AD. Remember to change this when changing searching radius!!!
                      #Comment from Xiaqu: If I understood correctly, the "0.1" here is the threshold of the search radius, and the 0.1 is around 11 km.
                      
                      ii_XY = Reduce(function(c1, c2) substitute(.c1 & .c2, list(.c1=c1, .c2=c2)), l_XY)  # Added by AD
                      
                      anacid_XY <- dat$cid.pre[dif_XY[eval(ii_XY), which=TRUE]]  # cids neighbor cells from future map # Added by AD
                      
                      #### for each focal cell subset target cell analogues (within ClimTol) WITHIN THE SPECIFIED DISTANCE # Added by AD
                      futs <- as.numeric(Dat[i, 2])
                      pres <- dat[dat$cid.pre %in% anacid_XY, 1] # Added by AD;Changed by XZ
                      dif <- data.table(sweep(pres, 2, futs, "-"))
                      
                      # Identify future analogue cells
                      upper = colnames(dif)
                      l <- lapply(upper, function(x) call("<=", call("abs", as.name(x)), climTol[grep(x, colnames(dif))]))
                      ii = Reduce(function(c1, c2) substitute(.c1 & .c2, list(.c1=c1, .c2=c2)), l)
                      anacid <- anacid_XY[dif[eval(ii), which=TRUE]]  # cids analogue cells/changed by XZ
                      
                      
                      # LOCATE CLOSEST ANALOGUE #Created by Xiaqu.
                      if(length(anacid)>0){
                        # check which of those are within distance and get the analogue at minimum distance
                        x <- dat$x_pre[dat$cid.pre %in% anacid]
                        match <- as.data.frame(cbind(Dat$x_fut[i],Dat$y_fut[i], dat$x_pre[dat$cid.pre %in% anacid], dat$y_pre[dat$cid.pre %in% anacid]))
                        colnames(match)[1] <- "from.x" #Change the column name of data frame
                        colnames(match)[2] <- "from.y"
                        colnames(match)[3] <- "to.x" #Change the column name of data frame
                        colnames(match)[4] <- "to.y"
                        match$d <- sqrt((match$from.x-match$to.x)^2 + (match$from.y-match$to.y)^2) # in x/y units
                        
                        
                        an <- anacid[match$d < geoTol]       # cids analogue cells within search radius
                        dis <- match$d[match$d < geoTol]     # distance to candidate analogues
                        if (length(an) > 0){
                          resu[i, target := an[which.min(dis)]]   # cid of geographically closest climate analogue
                          #resu[i, climDis := mean(as.numeric(dif[which(anacid == resu[i, target]),]))]  # mean clim difference for the closest analogue
                          resu[i, geoDis_m := min(dis)]
                          #resu[i, ang := geosphere::bearing(Dat[i, c("x","y")], dat[cid == resu[i, target], c("x","y")])]
                          resu[i, vel_m_yr := (resu$geoDis_m[i])/tdiff]
                        }
                      }
                      
                      #}
                      
                      # setTxtProgressBar(pb, i) #added by AD
                    }
                    return(resu)
                  }
#close(pb) #added by AD
stopCluster(cl)

#return(result)
# }
toc()
gc()

install.packages("writexl")
library(writexl)
write_xlsx(result, "E:/Output/VoCC/Results_xlsx_csv/BVoCC_NorthEU_2.xlsx", col_names = T, format_headers = T)
write.csv(result, "E:/Output/VoCC/Results_xlsx_csv/BVoCC_NorthEU_2.csv", row.names=FALSE)
#Generate the Euclidean fvocc raster map
Eucli_VoCC <- rast(fut_avai);names(Eucli_VoCC) <- 'Euclidean BVoCC'
Eucli_VoCC[result$focal] <-  result$vel_m_yr
Eucli_VoCC <- round(Eucli_VoCC, digits = 1 )

writeRaster(Eucli_VoCC, filename = paste0("/lustre1/scratch/348/vsc34871/output/VoCC/", number, "_BVoCC_NorthEU_Euclidean_SR10km_100m_ssp370.tif"), overwrite=TRUE)
BVoCC <- rast(paste0("/lustre1/scratch/348/vsc34871/output/VoCC/", number, "_BVoCC_NorthEU_Euclidean_SR10km_100m_ssp370.tif"))
print(BVoCC)
