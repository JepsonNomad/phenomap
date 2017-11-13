# phenomap

This takes in a series of multi-layer raster files and returns a phenology projection raster.

## Motivation

Though landscape phenology is increasingly becoming a focal point of investigations of migration timing, hitherto no R packages existed that were able to reconstruct satellite-derived phenological metrics in space.  Thus, it was my motivation to develop this package *en passant* with John 2016 in order to enable a broader group of researchers to project landscape timing measures in space.  

## Prerequisistes

This package uses methods from the `dplyr`, `phenex`, `plyr`, `raster`, `stringr`, and `rgdal` packages.  I recommend loading the `raster` package so that the product of `mapPheno()` can then be inspected and further processed as needed.

## Installation
```
devtools::install_github("JepsonNomad/phenomap")
```

## Code example
```
library(phenomap)
File_List.VI <- list.files(pattern="TeenyCrop")
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
```

To project the timing of senescence, use `phase = "senescence"`:
```
Sample.Senescence <- mapPheno(File_List = File_List.VI, PhenoFactor = PhenoFactor,
                           phase = "senescence", threshold = threshold, year = year,
                           NDVI = NDVI, VIQ = VIQ, DOY = DOY, PR = PR,
                           verbose = verbose)
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

## License 

This project is licensed under the GPL-3 license.
