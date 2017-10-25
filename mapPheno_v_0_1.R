# sessionInfo()
# library(phenex)
# library(rgdal)
# library(raster)
# library(stringr)
# library(plyr)
# library(dplyr)
# 
# Arguments of the mapPheno function:
# File_List = List raster files
# PhenoFactor = type of dataset to analyze (e.g., "VI", "Snow")
# phase = What phenophase are you monitoring
# threshold = GWI threshold of analysis
# year = Year of analysis
# NDVI = Band number of NDVI band in raster files
# VIQ = Band number of VI Quality layer in raster files
# DOY = Band number of Composite Day of Year layer in raster files
# PR = Band Number of PR layer in raster files
# verbose = verbose output (True/False)

mapPheno <- function(File_List = NA, PhenoFactor = NA, 
                     phase = NA, threshold = NA, year = NA,
                     NDVI = NA, VIQ = NA, DOY = NA, PR = NA, 
                     verbose = F){
  
  annualcrops <- File_List
  pos <- environment(mapPheno)
  
  if(verbose){print(paste0("Creating data tiles for each timepoint... ", Sys.time()))}
  
  # Create object for each datatype at each timepoint
  if(!is.na(NDVI) && !is.na(VIQ) && !is.na(DOY) && !is.na(PR)){
    for(i in 1:length(annualcrops)){
      TileDOY <- substr(annualcrops[i], start = 15, stop = 17)
      assign(x = paste("NDVITile", TileDOY, sep=""), 
             value = raster(annualcrops[i], band=NDVI),
             pos = pos)
      assign(x = paste("VIQTile", TileDOY, sep=""), 
             value = raster(annualcrops[i], band=VIQ),
             pos = pos)
      assign(x = paste("DOYTile", TileDOY, sep=""), 
             value = raster(annualcrops[i], band=DOY),
             pos = pos)
      assign(x = paste("PRTile", TileDOY, sep=""), 
             value = raster(annualcrops[i], band=PR),
             pos = pos)
    }
  }
  
  # Create lists of each datatype
  NDVI.List <- ls(pattern="NDVITile", pos = pos)[1:length(ls(pattern="NDVITile", pos = pos))]
  VIQ.List <- ls(pattern="VIQTile", pos = pos)[1:length(ls(pattern="VIQTile", pos = pos))]
  DOY.List <- ls(pattern="DOYTile", pos = pos)[1:length(ls(pattern="DOYTile", pos = pos))]
  PR.List <- ls(pattern="PRTile", pos = pos)[1:length(ls(pattern="PRTile", pos = pos))]
  
  ## This function transforms a list of strings returned by searching for objects
  ## into a list of objects (required for proper stacking)
  listobjecter <- function(x){
    FinalList <- lapply(X = x, FUN = get)
    return(FinalList)
  }
  
  print(paste0("Detected array; dimensions = "))
  print(dim(listobjecter(NDVI.List)[[1]]))
  
  if(verbose){print(paste0("Sorting data into 3-dimensional arrays... ", Sys.time()))}
  
  if(PhenoFactor == "VI"){
    
    if(verbose){print(paste0("VI option selected... ", Sys.time()))}
    
    if(is.na(NDVI) || is.na(VIQ) || is.na(DOY) || is.na(PR)){
      stop("Crucial Dataset missing (have you included NDVI, VIQ, DOY, and PR layers?)")
    }
    else{
      if(!is.na(NDVI)){
        if(verbose){print(paste0("Creating NDVI array... ", Sys.time()))}
        NDVI.Objects.List <- listobjecter(NDVI.List)
        NDVI.Stack <- stack(NDVI.Objects.List)
        NDVI.Array <- as.array(NDVI.Stack)
        
        ## If NDVI is less than 0, it needs to be brought up to 0
        NDVIcleaner <- function(x) {x[x<0] <- 0; x}
        NDVIcleaned <- apply(X = NDVI.Array, MARGIN = c(1,2,3), FUN = NDVIcleaner)
      }
      if(!is.na(VIQ)){
        if(verbose){print(paste0("Creating VI Quality array... ", Sys.time()))}
        VIQ.Objects.List <- listobjecter(VIQ.List); if(verbose){print("1")}
        VIQ.Stack <- stack(VIQ.Objects.List); if(verbose){print("2")}
        VIQ.Array <- as.array(VIQ.Stack); if(verbose){print("3")}
        
        first_k_bits <- function(int, k=16) {
          ## https://lpdaac.usgs.gov/products/modis_products_table/mod13q1 TABLE 2:
          ## MOD13Q1 VI Quality: "Bit 0 is the least significant (read bit words right to left)"
          integer_vector <- as.integer(intToBits(int))[1:k]
          return(paste(as.character(integer_vector), collapse=""))
        }
        
        VIQbinary <- apply(X = VIQ.Array, MARGIN = c(1,2,3), FUN = first_k_bits); print("4")
      }
      if(!is.na(DOY)){
        if(verbose){print(paste0("Creating Day of Year array... ", Sys.time()))}
        DOY.Objects.List <- listobjecter(DOY.List); if(verbose){print("1")}
        DOY.Stack <- stack(DOY.Objects.List); if(verbose){print("2")}
        DOY.Array <- as.array(DOY.Stack); if(verbose){print("3")}
      }
      if(!is.na(PR)){
        if(verbose){print(paste0("Creating Pixel Reliability array... ", Sys.time()))}
        PR.Objects.List <- listobjecter(PR.List); if(verbose){print("1")}
        PR.Stack <- stack(PR.Objects.List); if(verbose){print("2")}
        PR.Array <- as.array(PR.Stack); if(verbose){print("3")}
      }
      
      # Clean NDVI data based on VIQ and PR data
      PRcleaner <- function(NDVISet, PRSet){
        Data <- NDVISet
        Meta <- PRSet
        for (i in 1:length(Data)){
          if (!is.na(Meta[i])){
            if (Meta[i] > 1){
              Data[i] <- NA
            }
            if(Meta[i] == -1){
              Data[i] <- NA
            }
          }
        }
        return(Data)
      }
      VIQcleaner <- function(NDVISet, VIQSet){
        Data <- NDVISet
        Meta <- VIQSet
        MODLAND.QA <- as.numeric(substring(Meta, 1, 2))
        VI.Usefulness <- as.numeric(substring(Meta, 3, 6))
        Aerosol.Quantity <- as.numeric(substring(Meta, 7, 8))
        Adjacent.Clouds <- as.numeric(substring(Meta, 9))
        BRDF.Correction <- as.numeric(substring(Meta, 10))
        Mixed.Clouds <- as.numeric(substring(Meta, 11))
        Land.Water.Flag <- as.numeric(substring(Meta, 12, 14))
        Possible.Snow.Ice <- as.numeric(substring(Meta, 15))
        Possible.Shadow <- as.numeric(substring(Meta, 16))
        for (i in 1:length(Data)){
          if (MODLAND.QA[i] > 10){
            Data[i] <- NA
          }
          if (VI.Usefulness[i] > 1100){
            Data[i] <- NA
          }
        }
        return(Data)
      }
      
      if(verbose){print(paste0("Cleaning NDVI data... ", Sys.time()))}
      
      FinalNDVI.PR <- PRcleaner(NDVISet = NDVIcleaned, PRSet = PR.Array)
      FinalNDVI.VIQ <- VIQcleaner(NDVISet = FinalNDVI.PR, VIQSet = VIQbinary)
      FinalNDVI <- FinalNDVI.VIQ
      
      ## Prep your NDVI values for Phenex processing
      ## Begin by creating a 3-D array in which
      ## NDVI values will are sorted by date of composite creation
      if(verbose){print(paste0("Sorting NDVI data into complete annual space-time cube... ", Sys.time()))}
      year_array <- array(data = NA, dim = c(nrow(FinalNDVI), ncol(FinalNDVI), 366))
      for (p in 1:dim(FinalNDVI)[3]){
        zindex <- p
        for (i in 1:nrow(FinalNDVI)){
          rowindex <- i
          for (j in 1:ncol(FinalNDVI)){
            colindex <- j
            if (!is.na(DOY.Array[rowindex, colindex, zindex])){
              year_array[rowindex, colindex, DOY.Array[rowindex, colindex, zindex]] <- FinalNDVI[rowindex,colindex, zindex]
            }
          }
        }
      }
      
      if(verbose){
        print(paste0("Converting annual space-time cube into phenoscape... ", Sys.time()))
      }
      
      ## We want to apply the phenex functions through
      ## time in the FinalNDVI array;
      ## Meaning use aaply in dimensions c(1,2)
      NDVIFunction <- function(a){
        pixelID <- a
        ndvi <- modelNDVI((pixelID/10000), correction = "none",
                          method = "DLogistic", year=as.integer(year))[[1]]
        springtime <- phenoPhase(ndvi, phase = phase, method = "local", threshold = threshold)[[1]]
        return(springtime)
      }
      
      
      phenodump <- aaply(.data = year_array,
                         .fun = NDVIFunction,
                         .margins = c(1,2),
                         .progress = "text")
      
      
      if(verbose){print(paste0("Converting phenoscape to raster... ", Sys.time()))}
      
      ## Turn final array back into a raster based on original raster properties
      phenoscape <- raster(x = phenodump,
                           template = raster(File_List[[1]]))
      
      
    }
  }
  
  if(verbose) {print(Sys.time())}
  return(phenoscape)
}
