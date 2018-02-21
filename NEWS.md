# phenomap 1.2.0

## Major changes

* Added `mapTrend` function.  `mapTrend` identifies long-term, pixel-wise trends in phenology.  See function documentation for more.

## Minor changes

* Added glacier filter to `mapPheno`. This filter removes all pixels where sum(Possible.Snow.Ice == 1) > 18.
