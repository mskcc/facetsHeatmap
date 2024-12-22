# facetsHeatmap
Utility for visualization of copy number profile distributions in a clustered, multisample-cohort using 
the output from [facets-Suite](https://github.com/mskcc/facets-suite).

## Input data

The main input data file is a **cncf** file from [FACETS](https://github.com/mskcc/facets). 
There is an example file included,
[example.cncf](https://github.com/mskcc/facetsHeatmap/blob/inst/extdata/example.cncf):


    > library(facetsHeatmap)
    > x <- system.file('extdata','example.cncf', package = 'facetsHeatmap')
    > head(x)

sid          | chrom    | loc.start   | loc.end  | tcn.em    | lcn.em
-----------: | -------: | ----------: | -------: | --------: | --------:
s_67         | 1        | 1001177     | 16255758 | 2         |  0
s_67         | 1        | 16256100    | 17359660 | 2         |  1


Only the columns shown in the above are used.

## Processing copy number profiles

The segment information is processed into mean gain/loss matrices (samples in rows and positional grids in
columns) by

    > library(facetsHeatmap)
    > set.seed(315)
    > z <- mat2D(cncf = x, bin.size = 10, progress.bar = TRUE)
    > names(z)
    [1] "bins"    "matrix"  "dat.tcn" "dat.lcn"
    > head(z$bins)
    
id         | chromosome     | start    | end      |    fam   | floss
---------: | -------------: | -------: | -------: | -------: | -----:
1          |          1     |       1  |  9958257 |    0.37  |   0.07
2          |          1     | 9958258  | 19916514 |    0.35  |   0.05

The argument shown for **mat2D** are defaults, except **x**; **bin.size** is the size of the bin in mb. 
The output list component **bins** shows the mean fraction of gains **fam** and losses **floss** in each bin.
The component **matrix** is of dimension **n_sample x n_grids**: 

    > dim(z$matrix)
    [1] 100 287

    > head(z$matrix)
    
1         |     2   |   3  |   4  |    5  |   6 
--------- | ------- | ---- | ---- | ----: | ----: 
s_67      |  2.33   | 2.5  |   2  |    2  |   2  
s_37      |  4.00   | 4.0  |   4  |    4  |   3  
s_40      |  3.50   | 4.0  |  3.5 |    3  | 4.5

where each element is the median value of **tcn.em** from **cncf** input.

The components **dat.tcn** and **dat.lcn** are abbreviated copy number profiles (median values for each chromosome)
for **tcn** and **lcn**, respectively. These are used for clustering of samples.

## Clustering the cohort

Rudimentary wrappers for clustering the samples based on latent class analysis (LCA) are provided, calling the package
[poLCA](https://cran.r-project.org/web/packages/poLCA/index.html).

     > clusterZ(z, K = 1:5)


