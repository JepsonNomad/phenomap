#' Convert a series of raster files to a single phenology raster.
#'
#' @param File_List List of phenology raster files (i.e. those produced in `mapPheno`)
#' @param Year_List Integer Year (YYYY)
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
#' verbose = TRUE
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
#' @export
# mapTrend added to phenomap v1.1.6

library(raster)
library(plyr)
library(doParallel)

mapTrend <- function(File_List,
                     Year_List,
                     parallel = FALSE,
                     n.cores = NULL,
                     verbose = FALSE){
  if(parallel == TRUE){
    if(verbose){print("parallel option selected...")}
    if(verbose){print(paste0(n.cores, " cores selected"))}
    registerDoParallel(cores=n.cores)
    if(verbose){print("Parallel backend registered")}
  }

  myfiles <- File_List

  if(verbose){print(myfiles)}
  if(verbose){print(Sys.time())}

  phenostack <- stack(lapply(X = myfiles, MARGIN = 1, FUN = raster))
  if(verbose){
    print("Stack created")
  }

  phenostack.array <- raster::as.array(phenostack)
  if(verbose){
    print(str(phenostack.array))
  }

  # create function to determine pixel-specific linear trend through time series
  # requires minimum 5 non-NA cases
  trend.extraktor <- function(x){
    pixel.timeseries <- data.frame(x,Year_List)

    if(sum(!is.na(pixel.timeseries$x)) > 5){
      pixel.lm <- lm(pixel.timeseries$x ~ pixel.timeseries$Year_List)
      ts.coeff <- summary(pixel.lm)$coefficients[2,1]

      extraktmat <- ts.coeff

      return(extraktmat)
    }
    else{return(NA)}
  }
  sig.extraktor <- function(x){
    pixel.timeseries <- data.frame(x,Year_List)

    if(sum(!is.na(pixel.timeseries$x)) > 5){
      pixel.lm <- lm(pixel.timeseries$x ~ pixel.timeseries$Year_List)
      ts.sig <- summary(pixel.lm)$coefficients[2,4]

      extraktmat <- ts.sig

      return(extraktmat)
    }
    else{return(NA)}
  }
  se.extraktor <- function(x){
    pixel.timeseries <- data.frame(x,Year_List)

    if(sum(!is.na(pixel.timeseries$x)) > 5){
      pixel.lm <- lm(pixel.timeseries$x ~ pixel.timeseries$Year_List)
      ts.se <- summary(pixel.lm)$coefficients[2,2]

      extraktmat <- ts.se

      return(extraktmat)
    }
    else{return(NA)}
  }
  R2.extraktor <- function(x){
    pixel.timeseries <- data.frame(x,Year_List)

    if(sum(!is.na(pixel.timeseries$x)) > 5){
      pixel.lm <- lm(pixel.timeseries$x ~ pixel.timeseries$Year_List)
      ts.R2 <- summary(pixel.lm)$r.squared

      extraktmat <- ts.R2

      return(extraktmat)
    }
    else{return(NA)}
  }

  if(verbose){
    print("Trend extraktor created")
    print("Dumping regression metadata")
    print(Sys.time())
  }

  # extract the linear trend; dump to a trend array
  trend.array <- aaply(.data = phenostack.array,
                       .margins = c(1,2),
                       .fun = trend.extraktor,
                       .parallel = parallel)
  if(verbose){
    print("Trend data dumped")
    print(Sys.time())
  }

  sig.array <- aaply(.data = phenostack.array,
                     .margins = c(1,2),
                     .fun = sig.extraktor,
                     .parallel = parallel)
  if(verbose){
    print("Significance data dumped")
    print(Sys.time())
  }

  se.array <- aaply(.data = phenostack.array,
                    .margins = c(1,2),
                    .fun = se.extraktor,
                    .parallel = parallel)
  if(verbose){
    print("Standard error data dumped")
    print(Sys.time())
  }

  R2.array <- aaply(.data = phenostack.array,
                    .margins = c(1,2),
                    .fun = R2.extraktor,
                    .parallel = parallel)
  if(verbose){
    print("R^2 data dumped")
    print(Sys.time())
  }


  if(verbose){
    print("Regression metadata dumped")
    print("Converting trend.array to trend.raster")
    print(Sys.time())
  }

  trend.list <- list(trend.array,sig.array,se.array,R2.array)
  raster.list <- lapply(trend.list, raster, template = raster(File_List[[1]]))
  trend.raster <- stack(raster.list)

  names(trend.raster) <- c("Coefficient","P","Std.Error","R2")

  if(verbose){
    print("trend.raster created; function completed")
    print(Sys.time())
  }

  return(trend.raster)
}

