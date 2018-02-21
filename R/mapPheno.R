#' Convert a series of raster files to a single phenology raster.
#'
#' @param File_List List of raster files
#' @param PhenoFactor Character string; type of dataset to analyze (e.g., "VI", "Snow")
#' @param phase Character string; name of phenophase to be measured (e.g., "greenup", "snowmelt", "senescence" or other arguments passed to phenex::phenophase())
#' @param threshold Float threshold GWI value to be projected. Use only for VI option.
#' @param year Integer Year (YYYY)
#' @param NDVI Integer Band number of NDVI band in raster files
#' @param VIQ Integer Band number of VI Quality layer in raster files
#' @param DOY Integer Band number of Composite Day of Year layer in raster files
#' @param PR Integer Band Number of PR layer in raster files
#' @param SnowExtent Integer Band number of Maximum_Snow_Extent in raster files
#' @param parallel TRUE or FALSE (Default = FALSE) if TRUE, use parallel backend through plyr::aaply
#' @param n.cores Integer number of cores to be used for parallel processing (only use if parallel = TRUE)
#' @param verbose TRUE or FALSE (Default = FALSE)
#' @return Raster object with extent=extent(raster(File_List)[1]) and CRS = crs(raster(File_List)[1]).  Digital numbers are expressed as Day of Year.
#' @examples
#' \dontshow{
#' fpath <- system.file("extdata", package="phenomap")
#' File_List <- paste(fpath, list.files(path = fpath, pattern=c("ExCrop_")), sep="/")
#' File_List <- File_List[1:2]
#'
#' PhenoFactor = "VI"
#' phase = "greenup"
#' threshold = 0.5
#' year = 2016
#' NDVI = 1
#' VIQ = 3
#' DOY = 4
#' PR = 5
#' parallel = FALSE
#' n.cores = NA
#' verbose = FALSE
#'
#' Sample.Greenup <- mapPheno(File_List = File_List, PhenoFactor = PhenoFactor,
#'                            phase = phase, threshold = threshold, year = year,
#'                            NDVI = NDVI, VIQ = VIQ, DOY = DOY, PR = PR,
#'                            SnowExtent=SnowExtent,
#'                            parallel = parallel, n.cores = n.cores,
#'                            verbose = verbose)
#' }
#' \dontrun{
#'
#' fpath <- system.file("extdata", package="phenomap")
#' File_List <- paste(fpath, list.files(path = fpath, pattern=c("TinyCrop_")), sep="/")
#' File_List
#'
#' PhenoFactor = "VI"
#' phase = "greenup"
#' threshold = 0.5
#' year = 2016
#' NDVI = 1
#' VIQ = 3
#' DOY = 4
#' PR = 5
#' parallel = FALSE
#' n.cores = NA
#' verbose = TRUE
#'
#' Sample.Greenup <- mapPheno(File_List = File_List, PhenoFactor = PhenoFactor,
#'                            phase = phase, threshold = threshold, year = year,
#'                            NDVI = NDVI, VIQ = VIQ, DOY = DOY, PR = PR,
#'                            SnowExtent = SnowExtent,
#'                            parallel = parallel, n.cores = n.cores,
#'                            verbose = verbose)
#'
#' }
#' @import dplyr
#' @import rgdal
#' @import stringr
#' @importFrom raster raster
#' @importFrom raster stack
#' @importFrom plyr aaply
#' @importFrom phenex modelNDVI
#' @importFrom phenex phenoPhase
#' @importFrom doParallel registerDoParallel
#' @export
# version 1.1.4 - glacial filter added
mapPheno<- function(File_List = NA, PhenoFactor = NA,
                    phase = NA, threshold = NA, year = NA,
                    NDVI = NA, VIQ = NA, DOY = NA, PR = NA,
                    SnowExtent = NA,
                    parallel = FALSE, n.cores=NA,#v1.1.1 update; for aaply func; reqs doParallel
                    verbose = FALSE){

  if(parallel == TRUE){
    if(verbose){print("parallel option selected...")}
    if(verbose){print(paste0(n.cores, " cores selected"))}
    registerDoParallel(cores=n.cores)
    if(verbose){print("Parallel backend registered")}
    paropts <- list(.packages = "phenex")
  }

  annualcrops <- File_List
  pos <- new.env()

  ## This function transforms a list of strings returned by searching for objects
  ## into a list of objects (required for proper stacking)
  if(verbose){print(paste0("Sorting data into 3-dimensional arrays... ", Sys.time()))}

  if(PhenoFactor == "VI"){

    if(verbose){print(paste0("VI option selected... ", Sys.time()))}

    if(is.na(NDVI) || is.na(VIQ) || is.na(DOY) || is.na(PR)){
      stop("Crucial Dataset missing (have you included NDVI, VIQ, DOY, and PR layers?)")
    } # Check presence of crucial bands

    # Create 3D NDVI Array and preliminarily screen by moving negative values up to 0 (0 <= NDVI <= 1)
    if(!is.na(NDVI)){
      if(verbose){
        print(paste0("Creating NDVI array... ", Sys.time()))
        print("Array dimensions = ")
        print(c(dim(raster::as.matrix(raster(annualcrops[1])))[1],
                dim(raster::as.matrix(raster(annualcrops[1])))[2],
                length(annualcrops)))
      }
      NDVI.Array <- array(data=NA,
                          dim=c(dim(raster::as.matrix(raster(annualcrops[1]), dimnames = NULL))[1],
                                dim(raster::as.matrix(raster(annualcrops[1]), dimnames = NULL))[2],
                                length(annualcrops)),
                          dimnames=NULL)
      if(verbose){print("Created blank array")}
      for(i in 1:length(annualcrops)){
        NDVI.Array[,,i] <- raster::as.matrix(raster(annualcrops[i], band=NDVI))
      }; if(verbose){print("Filled array")}

      ## If NDVI is less than 0, it needs to be brought up to 0
      NDVIfilter <- function(x) {
        if(!is.na(x)){
          x[x<0] <- 0; x
        } else{x<-NA; return(x)}
      }; if(verbose){print("NDVI negative value filter created")}

      NDVI.filtered <- apply(X = NDVI.Array, MARGIN = c(1,2,3), FUN = NDVIfilter)
      if(verbose){print("NDVI cleared of negative values")}
    }

    if(verbose){raster::plot(raster(NDVI.filtered[,,5], template=raster(File_List[5])), main="NDVI timeslice 5")}

    # Create VIQ Binary array
    if(!is.na(VIQ)){
      if(verbose){print(paste0("Creating VI Quality array... ", Sys.time()))}
      VIQ.Array <- array(dim=c(dim(raster::as.matrix(raster(annualcrops[1])))[1],
                               dim(raster::as.matrix(raster(annualcrops[1])))[2],
                               length(annualcrops)),
                         data=NA); if(verbose){print("Created blank array")}
      for(i in 1:length(annualcrops)){
        VIQ.Array[,,i] <- raster::as.matrix(raster(annualcrops[i], band=VIQ))
      }
      if(verbose){
        print("Filled VIQ Array")
        print("VIQ.Array dimensions")
        print(dim(VIQ.Array))
        raster::plot(raster(VIQ.Array[,,5], template=raster(File_List[5])), main = "VIQ timeslice 5")
      }
      ## first_k_bits derived from:
      ## https://gis.stackexchange.com/questions/144441/how-can-i-parse-modis-mod13q1-quality-layers-in-r/144487
      first_k_bits <- function(int, k=16) {
        if(is.na(int)) {return(NA)}
        else{
          ## https://lpdaac.usgs.gov/products/modis_products_table/mod13q1 TABLE 2:
          ## MOD13Q1 VI Quality: "Bit 0 is the least significant (read bit words right to left)"
          integer_vector <- as.integer(intToBits(int))[1:k]
          return(paste(as.character(integer_vector), collapse=""))
        }
      } ; if(verbose){print("first_k_bits function created")}

      VIQbinary <- apply(X = VIQ.Array, MARGIN = c(1,2,3), FUN = first_k_bits)
      if(verbose){
        print("VIQ array translated to reverse binary; read left to right")
        print(VIQbinary[(nrow(VIQbinary)/2),((ncol(VIQbinary)/2)),])
      }
    }

    # Create DOY Array
    if(!is.na(DOY)){
      if(verbose){print(paste0("Creating Day of Year array... ", Sys.time()))}
      DOY.Array <- array(dim=c(dim(raster::as.matrix(raster(annualcrops[1])))[1],
                               dim(raster::as.matrix(raster(annualcrops[1])))[2],
                               length(annualcrops)),
                         data=NA); if(verbose){print("Empty DOY array created")}
      for(i in 1:length(annualcrops)){
        DOY.Array[,,i] <- raster::as.matrix(raster(annualcrops[i], band=DOY))
      }
      if(verbose){
        print("DOY array filled")
        raster::plot(raster(DOY.Array[,,5], template=raster(File_List[5])), main="DOY timeslice 5")
      }

    }

    # Create PR Array
    if(!is.na(PR)){
      if(verbose){print(paste0("Creating Pixel Reliability array... ", Sys.time()))}
      PR.Array <- array(dim=c(dim(raster::as.matrix(raster(annualcrops[1])))[1],
                              dim(raster::as.matrix(raster(annualcrops[1])))[2],
                              length(annualcrops)),
                        data=NA); if(verbose){print("Empty PR array created")}
      for(i in 1:length(annualcrops)){
        PR.Array[,,i] <- raster::as.matrix(raster(annualcrops[i], band=PR))
      }; if(verbose){print("PR array filled")}

    }

    # Clean NDVI data based on VIQ metadata
    VIQcleaner <- function(NDVISet, VIQSet){
      Data <- NDVISet
      Meta <- VIQSet
      MODLAND.QA <- as.numeric(substring(Meta, 1, 2))
      VI.Usefulness <- as.numeric(substring(Meta, 3, 6))
      Aerosol.Quantity <- as.numeric(substring(Meta, 7, 8))
      Adjacent.Clouds <- as.numeric(substring(Meta, 9, 9))
      BRDF.Correction <- as.numeric(substring(Meta, 10, 10))
      Mixed.Clouds <- as.numeric(substring(Meta, 11, 11))
      Land.Water.Flag <- as.numeric(substring(Meta, 12, 14))
      Possible.Snow.Ice <- as.numeric(substring(Meta, 15, 15))
      Possible.Shadow <- as.numeric(substring(Meta, 16, 16))

      for (i in 1:length(Data)){
        # if(verbose){print(Data[i])}
        if(is.na(Data[i])) {Data[i] <- NA}
        else{
          if (MODLAND.QA[i] > 10){
            Data[i] <- NA
          }
          if (VI.Usefulness[i] > 1100){
            Data[i] <- NA
          }
          if(Land.Water.Flag[i] == 0){ # shallow ocean
            Data[i] <- NA
          }
          if(Land.Water.Flag[i] == 101){ # deep inland water
            Data[i] <- NA
          }
          if(Land.Water.Flag[i] == 110){ # moderate or continental ocean
            Data[i] <- NA
          }
          if(Land.Water.Flag[i] == 111){ # deep ocean
            Data[i] <- NA
          }
          else{Data[i] <- Data[i]}
        }
      }
      if(verbose){print(Data[100:200])}
      return(Data)
    }

    if(verbose){print(paste0("Cleaning NDVI data based on VIQ layer... ", Sys.time()))}

    NDVI.VIQ <- VIQcleaner(NDVISet = NDVI.filtered, VIQSet = VIQbinary); if(verbose){print("NDVI data cleaned")}

    if(verbose){raster::plot(raster(NDVI.VIQ[,,5], template=raster(File_List[5])), main="Final NDVI timeslice 5")}

    if(verbose){print(paste0("Filtering out glacial features... ", Sys.time()))}

    GlacierFilter <- function(NDVISet, VIQSet){
      Data <- NDVISet
      Meta <- VIQSet

      for(i in 1:nrow(Meta)){
        for(j in 1:ncol(Meta)){
          Possible.Snow.Ice <- as.numeric(substring(Meta[i,j,], 15, 15))
          if(sum(Possible.Snow.Ice == 1) > 18){
            Data[i,j,] <- NA
          }
          else{
            Data[i,j,] <- Data[i,j,]
          }
        }
      }

      return(Data)
    }; if(verbose){print("GlacierFilter created")}

    NDVI.ice <- GlacierFilter(NDVISet = NDVI.VIQ, VIQSet = VIQbinary)

    if(verbose){
      print("Raster cleared of glaciers")
      raster::plot(raster(NDVI.ice[,,8],template=raster(File_List[5])),main="NDVI Frame 8,\nglaciers removed")
      raster::plot(raster(NDVI.ice[,,13],template=raster(File_List[5])),main="NDVI Frame 13,\nglaciers removed")
      raster::plot(raster(NDVI.ice[,,18],template=raster(File_List[5])),main="NDVI Frame 18,\nglaciers removed")
    }

    FinalNDVI <- NDVI.ice

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
            if(DOY.Array[rowindex, colindex, zindex] > 0){
              year_array[rowindex, colindex, DOY.Array[rowindex, colindex, zindex]] <- FinalNDVI[rowindex,colindex, zindex]
            }
          }
        }
      }
    }

    if(verbose){
      print(paste0("Converting annual space-time cube into phenoscape... ", Sys.time()))
      print("NDVI.Array pixel example")
      print(NDVI.Array[(nrow(VIQbinary)/2),((ncol(VIQbinary)/2)),])
      print("VIQ Processed NDVI pixel example")
      print(FinalNDVI[(nrow(VIQbinary)/2),((ncol(VIQbinary)/2)),])
      print("Composite Day of Year pixel example")
      print(DOY.Array[(nrow(VIQbinary)/2),((ncol(VIQbinary)/2)),])
      print("Annual space-time cube pixel example")
      print(year_array[(nrow(VIQbinary)/2),((ncol(VIQbinary)/2)),])
      raster::plot(raster(year_array[,,148],template=raster(File_List[5])),main="NDVI Julian Day 148")
      raster::plot(raster(year_array[,,149],template=raster(File_List[5])),main="NDVI Julian Day 149")
      raster::plot(raster(year_array[,,150],template=raster(File_List[5])),main="NDVI Julian Day 150")
      raster::plot(raster(year_array[,,151],template=raster(File_List[5])),main="NDVI Julian Day 151")
      raster::plot(raster(year_array[,,152],template=raster(File_List[5])),main="NDVI Julian Day 152")
      raster::plot(raster(year_array[,,153],template=raster(File_List[5])),main="NDVI Julian Day 153")
      raster::plot(raster(year_array[,,154],template=raster(File_List[5])),main="NDVI Julian Day 154")
      raster::plot(raster(year_array[,,155],template=raster(File_List[5])),main="NDVI Julian Day 155")
      raster::plot(raster(year_array[,,156],template=raster(File_List[5])),main="NDVI Julian Day 156")
      raster::plot(raster(year_array[,,157],template=raster(File_List[5])),main="NDVI Julian Day 157")
      raster::plot(raster(year_array[,,158],template=raster(File_List[5])),main="NDVI Julian Day 158")
      raster::plot(raster(year_array[,,159],template=raster(File_List[5])),main="NDVI Julian Day 159")
      raster::plot(raster(year_array[,,160],template=raster(File_List[5])),main="NDVI Julian Day 160")
      raster::plot(raster(year_array[,,161],template=raster(File_List[5])),main="NDVI Julian Day 161")
    }

    ## Create a matrix to be filled with phenology data
    phenodump <- array(data=NA, dim = c(nrow(FinalNDVI), ncol(FinalNDVI)))

    ## We want to apply the phenex functions through
    ## time in the FinalNDVI array;
    ## Meaning use aaply in dimensions c(1,2)
    NDVIFunction <- function(a){
      pixelID <- a
      if(sum(!is.na(pixelID) > 10)){ # ignore poorly represented pixels
        ndvi <- modelNDVI((pixelID/10000), correction = "none",
                          method = "DLogistic", year.int=as.integer(year))[[1]]
        springtime <- phenoPhase(ndvi, phase = phase, method = "local", threshold = threshold)[[1]]
      }
      else{springtime <- NA}
      return(springtime)
    }

    if(verbose){
      print("Phenex dependency initiated")
      print("Dumping threshold date")
    }

    phenodump <- aaply(.data = year_array,
                       .fun = NDVIFunction,
                       .margins = c(1,2),
                       #.progress = "text",
                       .parallel = parallel, .paropts = paropts)

    if(verbose){
      print("Phenoscape produced...")
      print(class(phenodump))
      print(dim(phenodump))
      print(paste0("Converting phenoscape to raster... ", Sys.time()))
    }

    ## Turn final array back into a raster based on original raster properties
    phenoscape <- raster(x = phenodump,
                         template = raster(annualcrops[[1]]))



  }

  if(PhenoFactor == "Snow"){
    if(verbose){print("Snowmelt option selected...")}
    if(is.null(SnowExtent)){
      stop("Crucial Dataset missing (have you included a SnowExtent layer?)")
    }
    Extent.Stack <- raster::stack(annualcrops, bands = SnowExtent)
    Extent.Array <- raster::as.array(Extent.Stack)

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
    phenodump.extent <- apply(.data = Complete.Extent.Array,
                              .margins = c(1,2),
                              .fun = function(x) min(which(x == 25)))

    phenoscape <- raster(x = phenodump.extent,
                         template = raster(annualcrops[[1]]))


  }
  if(verbose) {print(Sys.time())}
  return(phenoscape)
}

