# facetsHeatmap
Utility for visualization of copy number profile distributions in a clustered, multisample-cohort using 
the output from [facets-Suite](https://github.com/mskcc/facets-suite).

## Input data

The main input data file is a **cncf** file from [FACETS](https://github.com/mskcc/facets). 
There is an example file included,
[example.cncf](https://github.com/mskcc/facetsHeatmap/blob/inst/extdata/example.cncf):

    > library(facetsHeatmap)
    > x <- system.file('extdata','example.cncf', package = 'facetsHeatmap')

sid          | chrom    | loc.start   | loc.end  | tcn.em    | lcn.em
------------ | -------- | ----------- | -------- | --------- | ---------
s_67         | 1        | 1001177     | 16255758 | 2         |  0
s_67         | 1        | 16256100    | 17359660 | 2         |  1


Only the columns shown in the above are used.

## Processing copy number profiles

The segment information is processed into mean gain/loss matrices (samples in rows and positional grids in
columns) by

    > set.seed(315)
    > z <- mat2D(cncf = x, bin.size = 10, progress.bar = TRUE)
    > names(z)
    [1] "bins"    "matrix"  "dat.tcn" "dat.lcn"

