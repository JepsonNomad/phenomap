#' Convert a series of phenology terra::raster files to a single long-term trend terra::raster.
#'
#' @param File_List List of phenology terra::raster files (i.e. those produced in `mapPheno`)
#' @param Year_List Vector of Integer Year (YYYY) with length > 5
#' @param parallel TRUE or FALSE (Default = FALSE) if TRUE, use parallel backend through plyr::aaply
#' @param n.cores Integer number of cores to be used for parallel processing (only use if parallel = TRUE)
#' @param verbose TRUE or FALSE (Default = FALSE)
#' @return terra::raster object with extent=ext(rast(File_List)[1]) and CRS = crs(rast(File_List)[1]).  Layer 1 is the slope estimate of the linear model relating green-up timing (Day of Year) to time (Year).  Layer 2 is the p-value of the slope estimate.  Layer 3 is the standard error of the slope estimate.  Layer 4 is the r-squared value for the linear model.
#' @examples
#' \dontshow{
#' fpath <- system.file("extdata", package="phenomap")
#' File_List.Trend <- paste(fpath, list.files(path = fpath, pattern=c("Sample_Greenup_")), sep="/")[5:6]
#'
#' Year_List <- 2015:2016 # Tell it what years you're using
#' n.cores <- NA
#'
#' phenotrend <- mapTrend(File_List = File_List.Trend,
#'                              Year_List = Year_List,
#'                              parallel = FALSE,
#'                              n.cores = n.cores,
#'                              verbose=FALSE)
#'
#' }
#' \dontrun{
#'
#' fpath <- system.file("extdata", package="phenomap")
#' File_List.Trend <- paste(fpath, list.files(path = fpath, pattern=c("Sample_Greenup_")), sep="/")
#'
#' Year_List <- 2011:2016 # Tell it what years you're using
#' n.cores <- 4 # Set up parallel computing
#'
#' phenotrend <- mapTrend(File_List = File_List.Trend,
#'                              Year_List = Year_List,
#'                              parallel = TRUE,
#'                              n.cores = n.cores,
#'                              verbose=TRUE)
#'
#' }
#' @import stringr
#' @import terra
#' @importFrom plyr aaply
#' @importFrom doParallel registerDoParallel
#' @importFrom phenex phenoPhase
#' @export
# mapTrend added to phenomap v1.1.6

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

  phenostack <- terra::rast(lapply(X = myfiles,
                                   FUN = terra::rast))
  if(verbose){
    print("Stack created")
  }

  phenostack.array <- terra::as.array(phenostack)
  if(verbose){
    print(utils::str(phenostack.array))
  }

  # create function to determine pixel-specific linear trend through time series
  # requires minimum 5 non-NA cases
  trend.extraktor <- function(x){
    pixel.timeseries <- data.frame(x,Year_List)

    if(sum(!is.na(pixel.timeseries$x)) > 5){
      pixel.lm <- stats::lm(pixel.timeseries$x ~ pixel.timeseries$Year_List)
      ts.coeff <- summary(pixel.lm)$coefficients[2,1]

      extraktmat <- ts.coeff

      return(extraktmat)
    }
    else{return(NA)}
  }
  sig.extraktor <- function(x){
    pixel.timeseries <- data.frame(x,Year_List)

    if(sum(!is.na(pixel.timeseries$x)) > 5){
      pixel.lm <- stats::lm(pixel.timeseries$x ~ pixel.timeseries$Year_List)
      ts.sig <- summary(pixel.lm)$coefficients[2,4]

      extraktmat <- ts.sig

      return(extraktmat)
    }
    else{return(NA)}
  }
  se.extraktor <- function(x){
    pixel.timeseries <- data.frame(x,Year_List)

    if(sum(!is.na(pixel.timeseries$x)) > 5){
      pixel.lm <- stats::lm(pixel.timeseries$x ~ pixel.timeseries$Year_List)
      ts.se <- summary(pixel.lm)$coefficients[2,2]

      extraktmat <- ts.se

      return(extraktmat)
    }
    else{return(NA)}
  }
  R2.extraktor <- function(x){
    pixel.timeseries <- data.frame(x,Year_List)

    if(sum(!is.na(pixel.timeseries$x)) > 5){
      pixel.lm <- stats::lm(pixel.timeseries$x ~ pixel.timeseries$Year_List)
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
  trend.array <- plyr::aaply(.data = phenostack.array,
                             .margins = c(1,2),
                             .fun = trend.extraktor,
                             .parallel = parallel)
  if(verbose){
    print("Trend data dumped")
    print(Sys.time())
  }

  sig.array <- plyr::aaply(.data = phenostack.array,
                           .margins = c(1,2),
                           .fun = sig.extraktor,
                           .parallel = parallel)
  if(verbose){
    print("Significance data dumped")
    print(Sys.time())
  }

  se.array <- plyr::aaply(.data = phenostack.array,
                          .margins = c(1,2),
                          .fun = se.extraktor,
                          .parallel = parallel)
  if(verbose){
    print("Standard error data dumped")
    print(Sys.time())
  }

  R2.array <- plyr::aaply(.data = phenostack.array,
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
  raster.list <- lapply(trend.list, terra::rast)
  trend.raster <- terra::rast(raster.list)

  names(trend.raster) <- c("Coefficient","P","Std.Error","R2")

  if(verbose){
    print("trend.raster created; function completed")
    print(Sys.time())
  }

  return(trend.raster)
}

