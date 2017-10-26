#' Convert a series of raster files to a single phenology raster.
#'
#' @param File_List List of raster files
#' @param PhenoFactor Character string; type of dataset to analyze (e.g., "VI", "Snow")
#' @param phase Character string; name of phenophase to be measured
#' @param threshold Float threshold GWI value to be projected
#' @param year Integer Year (YYYY)
#' @param NDVI Integer Band number of NDVI band in raster files
#' @param VIQ Integer Band number of VI Quality layer in raster files
#' @param DOY Integer Band number of Composite Day of Year layer in raster files
#' @param PR Integer Band Number of PR layer in raster files
#' @param SnowExtent Integer Band number of Maximum_Snow_Extent in raster files
#' @param verbose T/F verbose
#' @return Raster object with extent=extent(raster(File_List)[1]) and CRS = crs(raster(File_List)[1]).
#' @examples
#' list.files()
#' File_List <- list.files(pattern=c("TinyCrop_"))
#'
#' PhenoFactor = "VI"
#' phase = "greenup"
#' threshold = 0.05
#' year = 2016
#' NDVI = 1
#' VIQ = 3
#' DOY = 4
#' PR = 5
#' verbose = TRUE
#'
#' Sample.Greenup <- mapPheno(File_List = File_List, PhenoFactor = PhenoFactor,
#'                            phase = phase, threshold = threshold, year = year,
#'                            NDVI = NDVI, VIQ = VIQ, DOY = DOY, PR = PR,
#'                            verbose = verbose)
#'
#' @import dplyr
#' @import rgdal
#' @import stringr
#' @importFrom raster raster
#' @importFrom raster stack
#' @importFrom plyr aaply
#' @importFrom phenex modelNDVI
#' @importFrom phenex phenoPhase
#' @export

mapPheno <- function(File_List = NA, PhenoFactor = NA,
                     phase = NA, threshold = NA, year = NA,
                     NDVI = NA, VIQ = NA, DOY = NA, PR = NA,
                     SnowExtent = NA,
                     verbose = F){

  annualcrops <- File_List
  pos <- new.env()

  ## This function transforms a list of strings returned by searching for objects
  ## into a list of objects (required for proper stacking)
  if(verbose){print(paste0("Sorting data into 3-dimensional arrays... ", Sys.time()))}

  if(PhenoFactor == "VI"){

    if(verbose){print(paste0("VI option selected... ", Sys.time()))}

    if(is.na(NDVI) || is.na(VIQ) || is.na(DOY) || is.na(PR)){
      stop("Crucial Dataset missing (have you included NDVI, VIQ, DOY, and PR layers?)")
    }
    if(!is.na(NDVI)){
      if(verbose){print(paste0("Creating NDVI array... ", Sys.time()))}
      NDVI.Stack <- raster::stack(raster(annualcrops, bands = NDVI))

      if(verbose){
        print(class(NDVI.Stack))
        print("NDVI.Stack dimensions ")
        print(dim(NDVI.Stack))
        print("2")
        }
      NDVI.Array <- as.array(NDVI.Stack)
      if(verbose){
        print(dim(NDVI.Array))
        print("3")
        }


      ## If NDVI is less than 0, it needs to be brought up to 0
      NDVIcleaner <- function(x) {x[x<0] <- 0; x};if(verbose){print("cleaner created")}
      NDVIcleaned <- apply(X = NDVI.Array, MARGIN = c(1,2,3), FUN = NDVIcleaner)
    };if(verbose){print("NDVI cleaned")}
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

      VIQbinary <- apply(X = VIQ.Array, MARGIN = c(1,2,3), FUN = first_k_bits); if(verbose){print("4")}
    }
    if(!is.na(DOY)){
      if(verbose){print(paste0("Creating Day of Year array... ", Sys.time()))}
      DOY.Objects.List <- listobjecter(DOY.List); if(verbose){print("1")}
      DOY.Stack <- stack(DOY.Objects.List); if(verbose){print("2")}
      DOY.Array <- as.array(DOY.Stack); if(verbose){print("3")}
    }
    if(!is.na(PR)){
      if(verbose){print(paste0("Creating Pixel Reliability array... ", Sys.time()))}
      PR.Objects.List <- listobjecter(PR.List)
      PR.Stack <- stack(PR.Objects.List)
      PR.Array <- as.array(PR.Stack)
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
        if (MODLAND.QA[i] > 5){
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
      #print(year_array[120:130,230:240,42])
    }

    ## Create a matrix to be filled with phenology data
    phenodump <- array(data=NA, dim = c(nrow(FinalNDVI), ncol(FinalNDVI)))

    ## We want to apply the phenex functions through
    ## time in the FinalNDVI array;
    ## Meaning use aaply in dimensions c(1,2)
    NDVIFunction <- function(a){
      pixelID <- a
      ndvi <- modelNDVI((pixelID/10000), correction = "none",
                        method = "DLogistic", year=as.integer(year))[[1]]
      springtime <- phenoPhase(ndvi, phase = phase, method = "local", threshold = threshold)
      return(springtime)
    }


    if(verbose){print(paste0("Converting phenoscape to raster... ", Sys.time()))}
    phenodump <- aaply(.data = year_array,
                       .fun = NDVIFunction,
                       .margins = c(1,2),
                       .progress = "text")


    ## Turn final array back into a raster based on original raster properties
    phenoscape <- raster(x = phenodump,
                         template = NDVI.Objects.List[[1]])



  }

  if(PhenoFactor == "Snow"){
    if(is.null(SnowExtent)){
      stop("Crucial Dataset missing (have you included a SnowExtent layer?)")
    }
    if(!is.null(Extent.List)){
      Extent.Objects.List <- listobjecter(Extent.List)
      Extent.Stack <- stack(Extent.Objects.List)
      Extent.Array <- as.array(Extent.Stack)

      Complete.Extent.Array <- array(dim = c(dim(Extent.Array)[1], dim(Extent.Array)[2], 366))
      #Sort Extent.Array into daily slices (rather than 8-day slices)
      for(p in (1:dim(Extent.Array)[3])){
        zindex <- p
        zfinal <- p*8-7
        for(q in (1:nrow(Extent.Array))){
          rowindex <- q
          for(r in (1:ncol(Extent.Array))){
            colindex <- r
            Complete.Extent.Array[rowindex,colindex,zfinal] <- Extent.Array[rowindex,colindex,zindex]
          }
        }
      }

      # Create a phenodump array to deposit snowmelt phenology data
      phenodump.extent <- array(data=NA, dim = c(nrow(Complete.Extent.Array), ncol(Complete.Extent.Array)))

      # Fill in phenodump with the first date at which x == 25, that is,
      # the first time slice when pixels are detected as completely snow-free
      phenodump.extent <- apply(X = Complete.Extent.Array,
                                MARGIN = c(1,2),
                                FUN = function(x) min(which(x == 25)))

      phenoscape <- raster(x = phenodump.extent,
                           template = Extent.Objects.List[[1]])
    }

  }
  if(verbose) {print(Sys.time())}
  return(phenoscape)
}

#
# setwd("/Users/christianjohn/Documents/Post_Lab/Phenoscape/test")
# list.files()
# File_List <- list.files(pattern=c("TinyCrop_"))
# PhenoFactor = "VI"
# phase = "greenup"
# threshold = 0.5
# year = 2016
# NDVI = 1
# VIQ = 3
# DOY = 4
# PR = 5
# verbose = TRUE
#
# library(phenex)
# library(rgdal)
# library(raster)
# library(stringr)
# library(plyr)
# library(dplyr)
#
#
# Sample.Greenup <- mapPheno(File_List = File_List, PhenoFactor = PhenoFactor,
#                            phase = phase, threshold = threshold, year = year,
#                            NDVI = NDVI, VIQ = VIQ, DOY = DOY, PR = PR,
#                            verbose = verbose)
#
