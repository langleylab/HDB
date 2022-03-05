# HDB
Hausdorff Distance of Batches


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

An alternative plot is the pyramid plot:

```{r}
plotHDSigmas(hdb, "individual")
plotHDSigmas(hdb, "replicate")
```

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

We can see that replicate 2 (the orange group of cells), being absent from the group belonging to individual `NA19098`, is very distant from either groups of blue points (replicate 1), albeit closer to the group from `NA19239`. This means that the r2 group in `NA19239` (squares) may establish a NN (nearest neighbor) pair with either r1 from `NA19239` (squares) or r1 from `NA19098` (circles), whereas the r1 groups have only one choice (r2 in `NA19239`). 
