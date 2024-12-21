# facetsHeatmap
Utility for visualization of copy number profile distributions in a clustered, multisample-cohort using 
the output from [facets-Suite](https://github.com/mskcc/facets-suite).

## Input data

The main input data file is a **cncf** file from FACETS. There is an example file included,
[example.cncf](https://github.com/mskcc/facetsHeatmap/blob/inst/extdata/example.cncf):

    > library(facetsHeatmap)
    > x <- system.file('extdata','example.cncf', package = 'facetsHeatmap')

sid          | chrom    | loc.start   | loc.end  | tcn.em    | lcn.em
------------ | -------- | ----------- | -------- | --------- | ---------
s_67         | 1        | 1001177     | 16255758 | 2         |  0
s_67         | 1        | 16256100    | 17359660 | 2         |  1


Only the columns shown in the above are used.
