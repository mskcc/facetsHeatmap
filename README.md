# facetsHeatmap
Utility for visualization of copy number profile distributions in a clustered, multisample-cohort using 
the output from [facets-Suite](https://github.com/mskcc/facets-suite).

## Installation

Install it with

    > devtools::install_github('mskcc/facetsHeatmap')
    
## Input data

The main input data file is a **cncf** file from [FACETS](https://github.com/mskcc/facets). 
There is an example file included,
[example.cncf](https://github.com/mskcc/facetsHeatmap/blob/inst/extdata/example.cncf):


    > library(facetsHeatmap)
    > x <- system.file('extdata','example_with_ploidy.cncf', package = 'facetsHeatmap')
    > head(x)

sid          | chrom    | loc.start   | loc.end  | tcn.em    | lcn.em  | ploidy
-----------: | -------: | ----------: | -------: | --------: | -------:| ------:
s_67         | 1        | 1001177     | 16255758 | 2         |  0      |    2.1
s_67         | 1        | 16256100    | 17359660 | 2         |  1      |    2.1


Only the columns shown in the above are used. The last column **ploidy** is optional, and if present, 
will be used when showing relative copy number profile.

## Processing copy number profiles

The segment information is processed into mean gain/loss matrices (samples in rows and positional grids in
columns) by

    > library(facetsHeatmap)
    > set.seed(315)
    > z <- mat2D(cncf = x, bin.size = 10, amp.tcn = 4, progress.bar = TRUE)
    > names(z)
    [1] "bins"    "matrix"  "matrix.loh" "dat.tcn" "dat.lcn" "ploidy"
    > head(z$bins)
    
id         | chromosome     | start    | end      |    fgain  |   famp |  floss  | fdel   | floh 
---------: | -------------: | -------: | -------: | --------: | ------:| ------: | -----: | ------
1          |          1     |       1  |  9958257 |    0.15   |   0.22 |   0.07  | 0      |  0.18
2          |          1     | 9958258  | 19916514 |    0.12   |   0.18 |   0.04  | 0      |  0.12

The argument shown for **mat2D** are defaults, except **cncf**; **bin.size** is the size of the bin in mb. 
The output list component **bins** shows the frequencies of samples with gains, amp (tcn > **amp.tcn**), loss, homdel, 
and LOH in each bin. The component **matrix** is of dimension **n_sample x n_grids**: 

    > dim(z$matrix)
    [1] 100 287

    > head(z$matrix)
    
id        |       1  |  2    |   3   |     4  |   5 
--------- | -------: | ----: | ----: | ----: | ----: 
s_67      |  2.33    | 2.5   |   2   |     2  |   2  
s_37      |  4.00    | 4.0   |   4   |     4  |   3  
s_40      |  3.50    | 4.0   |  3.5  |     3  | 4.5

where each element is the mean value of **tcn.em** from **cncf** input. The matrix **matrix.loh** contains the density
of LOH given by **lcn.em** = 0.

The components **dat.tcn** and **dat.lcn** are abbreviated copy number profiles (median values for each chromosome)
for **tcn** and **lcn**, respectively. These are used for clustering of samples.

The component **ploidy** stores ploidy values of all unique samples (**NA** if the column was not in input):

    > head(z$ploidy)
    s_67 s_37 s_40  s_9 s_90 s_51 
    2.1  3.3  3.0  1.6  2.7  2.3 

## Clustering the cohort

Rudimentary wrappers for clustering the samples based on latent class analysis (LCA) are provided, calling the package
[poLCA](https://cran.r-project.org/web/packages/poLCA/index.html).

     > clusterZ(z, K = 1:5)

With multiple **K** (no. of clusters) values as above, return value is a **ggplot** object of BIC vs. **K**:

<img src="/inst/extdata/cluster.png" alt="drawing" width="400" height="auto"/>

## Generate heatmap

We choose **K = 3** below, specifying it in the argument **K** for **showHeatmap**:

     > plt <- showHeatmap(z, K = 3, margin=c(0.15, 0.5, 0.5, 1.3), mar = c(1.5, 2.5, 1, 2.6))

The output is a ggplot object for the heatmap combined with reduced 1D average:

     > ggsave(plt, filename = 'heatmap.png', device = 'png')

<img src="/inst/extdata/heatmap.png" alt="" width="1000"/>

By default, it will cluster the samples with **K**, reorder them based on cluster membership (labeled 'C1', 'C2', ..., 'CK')
and by hierarchical clustering within clusters. The order of clusters are with increasing mean **tcn**. 
In the 1D profile, the line plot for LOH fractions can be turned off with **show.loh = FALSE**. To skip clustering (e.g. with
too few samples) and/or 
reordering within clusters, use the options **clusterZ = FALSE** and/or **reorder.hc = FALSE**, and the samples will be
displayed with the input order.

One can specify the maximum 'tcn' with argument **tcn.max** (default 6). If the 1D and 2D part don't align, try to fudge with
changing the arguments **margin** and **mar** (above are the defaults), each parameters for **ggplot** of 2D and base R graphics
of 1D, respectively; formats are margin = c(top, right, bottom, left) and mar = c(bottom, left, top, right).

One can use relative copy numbers for the heatmap to account for whole-genome duplication:

     > plt2 <- showHeatmap(z, K = 3, relative = TRUE, margin=c(0.15, 0.32, 0.5, 1.3), tcn.max = 4, tcn.min = -4)
     > ggsave(plt2, filename = 'heatmap2.png', device = 'png')

<img src="/inst/extdata/heatmap2.png" alt="" width="1000"/>

where for each sample, known ploidy (**z$plolidy**) is subtracted from local values 
(if unknown, mean tcn is computed from the matrix and subtracted). Heatmaps are centered at zero (white color).
     
Instead of **tcn**, LOH distribution given by fraction of **lcn.em == 0** can be displayed by specifying **useLOH = TRUE**:

     > plt3 <- showHeatmap(z, K = 3, useLOH = TRUE, margin = c(0.15, 0.3, 0.5,1.3))
     > ggsave(plt3, filename='heatmap3.png', device='png')

<img src="/inst/extdata/heatmap3.png" alt="" width="1000"/>
