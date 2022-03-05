# HDB

## Hausdorff Distance of Batches

### Installation

For the time being, clone the repo, go to the repo directory and launch R (or `setwd()` from R). 

There, use `devtools::load_all()` to load the package.

### Usage

The package contains one main function, `HDB()`, and two plotting functions (see below). 

To run `HDB()` you need a `SingleCellExperiment` object (`sce`) for which there is batch information in the `colData(sce)` slot and a dimensionality reduction is already present in `reducedDim(sce)`. Parallelization is possible through `BiocParallel`. 

```{r}
hdb <- HDB(sce,  
           group = "batch", 
           q = 4,
           verbose = TRUE,
           doplot = FALSE)
```           
           
### Short tutorial

Here is a simple demonstration on how to use `HDB`. 

We acquire a simple dataset following the `kBET` [vignette](https://github.com/theislab/kBET/blob/master/vignettes/kBET_vignette.Rmd). 

This is a small dataset with cells from 3 different individuals, with 3 different replicates for each.

```{r}
library(RCurl)
library(scran)
library(scater)

# Obtain dataset
url <- "https://raw.githubusercontent.com/jdblischak/singleCellSeq/master/data/"
umi <- read.table(text=getURL(paste0(url, "molecules.txt")))
anno <- read.table(text=getURL(paste0(url,"annotation.txt")), sep = "\t", header = TRUE)

# Get quality single cells
quality_single_cells <- scan(text=getURL(paste0(url,"quality-single-cells.txt")), 
                             what = "character")

# Create SingleCellExperiment
sce <- SingleCellExperiment(assays = list(counts = umi), colData = anno)
sce <- sce[, quality_single_cells]
sce <- sce[rowSums(counts(sce)) > 0,]
spikes <- grep('ERCC', rownames(sce))

#Normalize counts
sce <- logNormCounts(sce)

#Run PCA
sce <- runPCA(sce)
```

We see that individual `NA19098` has no cells for replicate `r2`:

```{r}
table(sce$individual, sce$replicate)

         
#          r1 r2 r3
#  NA19098 85  0 57
#  NA19101 80 70 51
#  NA19239 74 68 79
```
We plot the PCA colouring by `individual` and shaping by `replicate`:

```{r}
plotPCA(sce, colour_by = "individual", shape_by = "replicate")
```
![p1](https://user-images.githubusercontent.com/21171362/156878560-8fdb6b94-f089-467f-bddf-12d794ba0dc2.png)


The batch effect is evident at the individual level and at the replicate level. 

Moreover, the replicate 2 is entirely missing in one sample, which is an interesting use case for `HDB`.

Now we apply `HDB`:

```{r}
hdb <- HDB(sce,  
           group = c("individual", "replicate"), 
           q = 4,
           verbose = TRUE,
           doplot = FALSE)
           
#Computing null distributions for group individual
#
#Computing null models for group individual
#
#Computing pHds for group individual
#
#Computing null distributions for group replicate
#
#Computing null models for group replicate
#
#Computing pHds for group replicate
```

We plot the results:

```{r}
plotHDheatmap(hdb, "individual")
plotHDheatmap(hdb, "replicate")
```
![p2](https://user-images.githubusercontent.com/21171362/156878569-63dec937-d105-4349-a8cf-6e9b8aae5a2d.png)
![p3](https://user-images.githubusercontent.com/21171362/156878571-978155f5-71e6-4a94-a392-c9497417e017.png)

An alternative plot is the pyramid plot:

```{r}
plotHDSigmas(hdb, "individual")
plotHDSigmas(hdb, "replicate")
```
![p4](https://user-images.githubusercontent.com/21171362/156878578-7de5e1df-3db0-411d-a1a8-759916485146.png)
![p5](https://user-images.githubusercontent.com/21171362/156878581-4e9347e6-8f00-45f0-9650-f9571c401c43.png)

How do we interpret these numbers and plots? 

In HDB we introduce the concept of **asymmetric** effect, i.e. an effect in which one sample A occupies the same portion of the transcriptional space as another sample B, plus an additional part of space that B is not populating. The partial Hausdorff distance is in fact not symmetric, as it allows to measure the asymmetric coverage with high precision. 

A biological interpretation of an asymmetric effect is differential composition: if A contains a cell type that is depleted in B, and assuming integration was performed correctly, then the pHD of A to B will be higher than the pHD of B to A. 

For this reason both heatmaps and pyramid plots show the asymmetric effect specifying the direction: 
- "A to B" and "B to A"" in the barplot
- columns to rows in the heatmap

As you can see in the **replicate** level batch effect, **r1** and **r2** have an asymmetric effect: **r1** is further away from **r2** (r1 to r2) than **r2** is from **r1** (r2 to r1).

To have an intuition on why this happens, we look at the PCA plot, but we colour it by replicates and shape by individuals:

```{r}
plotPCA(sce, 
        colour_by = "replicate", 
        shape_by = "individual")
```

![p6](https://user-images.githubusercontent.com/21171362/156878586-7fcece6a-467b-4007-b8a1-1ea0f9849261.png)


We can see that replicate 2 (the orange group of cells), being absent from the group belonging to individual `NA19098`, is very distant from either groups of blue points (replicate 1), albeit closer to the group from `NA19239`. This means that the r2 group in `NA19239` (squares) may establish a NN (nearest neighbor) pair with either r1 from `NA19239` (squares) or r1 from `NA19098` (circles), whereas the r1 groups have only one choice (r2 in `NA19239`). 
