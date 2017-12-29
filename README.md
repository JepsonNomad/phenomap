# phenomap

This takes in a series of multi-layer raster files and returns a phenology projection raster.  `phenomap` is capable of analyzing satellite-derived NDVI and snowmelt time series on a regional or global scale.  Currently 'phenomap' is optimized for data that derives from the MOD 13 and MOD 10 data products, but it is expected that this project will eventually expand to also support GIMMS AVHRR and perhaps Landsat datasets.

## Motivation

Though landscape phenology is increasingly becoming a focal point of investigations of migration timing, hitherto no R packages existed that were able to reconstruct satellite-derived phenological metrics in space.  Thus, it was my motivation to develop this package *en passant* with John 2016 (M.S. thesis in ecology) in order to enable a broader group of researchers to project landscape timing measures in space.  

## Prerequisites

`phenomap` requires at least R version 3.4.2, and has been tested on Mac OSX El Capitan v10.11.6 and Windows 7.

This package uses methods from the `dplyr`, `phenex`, `plyr`, `raster`, `stringr`, and `rgdal` packages.  I recommend loading the `raster` package so that the product of `mapPheno()` can then be inspected and further processed as needed.

`phenomap` does not support .hdf files, and therefore it is recommended that users download and mosaic data.  I find that [pyModis](https://github.com/lucadelu/pyModis) is a reliable tool for such pursuits.  

## Warning

`phenomap` requires a considerable amount of computing power and even for small datasets it may take a while to process the data.  I therefore recommend cropping datasets to the minimum possible extent and projecting to the coarsest possible resolution.  This dev version of `phenomap` aims to support parallel processing, but it is expected to be especially buggy in early versions, but improve in future versions.

## Installation

Use devtools to install this package...
```
devtools::install_github("JepsonNomad/phenomap")
```
... until it is available on CRAN, at which point install `phenomap` using
```
install.packages(phenomap)
```

## Code example
```
library(phenomap)
File_List.VI <- list.files(pattern="TeenyCrop") # these files can be found in the /inst/extdata folder
PhenoFactor = "VI"
threshold = 0.5
year = 2016
NDVI = 1
VIQ = 3
DOY = 4
PR = 5
verbose = TRUE
```
To project the timing of 50% GWI, use mapPheno with `phase = "greenup"`, and `threshold = 0.5`.
```
Sample.Greenup <- mapPheno(File_List = File_List.VI, PhenoFactor = PhenoFactor,
                           phase = "greenup", threshold = threshold, year = year,
                           NDVI = NDVI, VIQ = VIQ, DOY = DOY, PR = PR,
                           verbose = verbose)
plot(Sample.Greenup)
```

To project the timing of senescence, use `phase = "senescence"`:
```
Sample.Senescence <- mapPheno(File_List = File_List.VI, PhenoFactor = PhenoFactor,
                           phase = "senescence", threshold = threshold, year = year,
                           NDVI = NDVI, VIQ = VIQ, DOY = DOY, PR = PR,
                           verbose = verbose)
plot(Sample.Senescence)
```

To inspect broad measures of timing, compare histograms:

```
breaks<-seq(from=0,to=365,by=5)
hist(as.matrix(Sample.Senescence), breaks=breaks, ylim=c(0,150),
     col="ORANGE", main = "Snowmelt, Greenup, and Senescence",
     xlab="Day of Year")
hist(as.matrix(Sample.Greenup), breaks=breaks,
     col="palegreen", add=T)
```

To project the duration of the growing season, compare greenup and senescence rasters:

```
Growing.Season <- as.matrix(Sample.Senescence) - as.matrix(Sample.Greenup)
GS.raster <- raster(Growing.Season, template = Sample.Greenup)
plot(GS.raster)
```

## License 

This project is licensed under the GPL-3 license.
