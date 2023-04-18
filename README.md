# HDB

## Hausdorff Distance of Batches

`HDB` is a package to diagnose batch effects in single cell RNA-seq data, although in principle it can be used more in general to estimate the overlap (or lack thereof) of sets of point clouds in a high dimensional space. 

`HDB` is built on the notion of the [Hausdorff distance](https://en.wikipedia.org/wiki/Hausdorff_distance), which can be defined as the largest distance between two nearest neighbors of two sets of points. 

In this package we use the **partial** and **robust** implementation of the Hausdorff distance:
- **partial** (or directed) means we take the distance from A to B separately from the distance from B to A, as they highlight two different overlaps. 
- **robust** means we do not take the absolute largest distance - as it can be greatly influenced by outliers - but rather the **k**-th largest distance. 

The brunt of the approach in this package lies in the estimation of a null distance distribution to which the observed partial Hausdorff distance is compared, which also allows for statistical testing.

Pending benchmarks on its speed and memory usage, the advantages of using HDB are the following:

- it allows the estimation of asymmetric effects (such as additional cell types or cell states) as separate phenomena from the lack of overlap
- it allows the quantification of effects with no upper bound, i.e. it has a large dynamic range  
- it only requires one round of null model estimation for most groupings
- it is flexible enough to be used for cluster labels, and can be easily extended to a single cell estimation (still under development)

### Etymology

The name `HDB` is both an acronym for **Hausdorff Distance of Batches**, and a callback to Singapore's public housing, in which neighbors of different walks of life and communities come together to create the *kampong* spirit. 

It is, in a way, a real life example of *nearest neighbors from different backgrounds*.


## Installation

Use `remotes::install_github("langleylab/HDB")` to install the package.

## Usage

The package contains one main function, `HDB()`, and two plotting functions (see below). 

To run `HDB()` you need a `SingleCellExperiment` object (`sce`) for which there is batch information in the `colData(sce)` slot and a dimensionality reduction is already present in `reducedDim(sce)`. Parallelization is possible through `BiocParallel`. 

```{r}
hdb <- HDB(sce,  
           group = "batch", 
           q = 4,
           verbose = TRUE,
           doplot = FALSE)
```           

For all the arugments, check the manual with `?HDB`.
           
## Short tutorial

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
plotHDsigmas(hdb, "individual")
plotHDsigmas(hdb, "replicate")
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


## Model

The batch effect is estimated as follows: for every couple of samples, we estimate the robust partial Hausdorff distance (i.e. the largest distance between two nearest neighbors) in a low-dimensional embedding such as PCA.

As mentioned earlier, the partial Hausdorff distance is asymmetric, i.e. pHD(A,B) is different from pHD(B,A). 

The pHD is then compared to a "null" distribution, i.e. the pHD for to well-mixed high-dimensional point sets, in the absence of systematic differences that drive the lack of overlap. We call this the "null pHD".

The null pHD is still influenced by two key aspects of the points and their coordinates: 1) their relative proportion (e.g. A is 20% of the total of A and B, meaning B is 80%) and 2) their joint density (the number of points in the high-dimensional volume). 

pHD from a high density sample to a low density sample will always be on average lower than its opposite even when the two samples are well mixed, i.e. they occupy the same space in a well distributed manner.

To convince ourselves, we simulate two 10-dimensional Gaussian distributions, with equal mean and variance, the difference being that the X distribution is made of 5000 points and the Y distribution is made of 2000 points.

```{r}
set.seed(123)
X = matrix(rnorm(5000 * 10, mean = 0, sd = 1), ncol = 10)
Y = matrix(rnorm(2000 * 10, mean = 0, sd = 1), ncol = 10)
```

In `HDB`, we defined the internal function `.pHD` as follows:

```{r}
.pHD <- function(A, B, q, distance="Euclidean",
                 BNPARAM=KmknnParam(distance=distance)){
  nn = queryKNN(A, B, k=1, get.index=FALSE, BNPARAM=BNPARAM)
  k = length(nn$distance) - q
  return(sort(nn$distance)[k])
}
```

Meaning it will take a reference dataset `A`, a query dataset `B`, and return the robust (rank = `1 - q`) partial Hausdorff distance from B to A. 

```{r}
.pHD(X, Y, q = 3)
[1] 2.892433

.pHD(Y, X, q = 3)
[1] 3.223152
```
At this point we can generate a number of random distributions and see what the pHD looks like in either direction:

```{r}
set.seed(123)
seeds = sample(seq_len(1000), size = 100, replace = FALSE)

rdists = lapply(seeds, function(x) {
  set.seed(x)
  X = matrix(rnorm(5000 * 10, mean = 0, sd = 1), ncol = 10)
  Y = matrix(rnorm(2000 * 10, mean = 0, sd = 1), ncol = 10)
  return(list(X = X, Y = Y))
})

pHD_XY = lapply(rdists, function(l) {
  .pHD(l$X, l$Y, q = 3)
})

pHD_YX = lapply(rdists, function(l) {
  .pHD(l$Y, l$X, q = 3)
})
```

Now we can visualize the distribution of both pHDs from 100 random 10-d Gaussians:

```{r}
 boxplot(list(xy = unlist(pHD_XY), yx = unlist(pHD_YX)), las = 1)
```
![image](https://user-images.githubusercontent.com/21171362/232579463-dfd03202-d5d8-4e00-9334-54679f43a360.png)

As you can see, the Y to X pHD is systematically larger than X to Y, i.e. a higher density sample (X) will have to travel less to find a nearest neighbor in a lower density sample (Y) than vice versa. 

This calculation only takes into account the *relative proportion* of point sets, since they occupy the same volume by construction. However, some real life datasets will not occupy the same volume, as they may be "more scattered" than others even within the same embedding. 

We simulate this by looking at the pHD in a similar setting, only this time we change the variance of the 10-dimensional Gaussian:


```{r}
set.seed(123)
seeds = sample(seq_len(1000), size = 100, replace = FALSE)

rdists2 = lapply(seeds, function(x) {
  set.seed(x)
  X = matrix(rnorm(5000 * 10, mean = 0, sd = 4), ncol = 10)
  Y = matrix(rnorm(2000 * 10, mean = 0, sd = 4), ncol = 10)
  return(list(X = X, Y = Y))
})

pHD_XY2 = lapply(rdists2, function(l) {
  .pHD(l$X, l$Y, q = 3)
})

pHD_YX2 = lapply(rdists2, function(l) {
  .pHD(l$Y, l$X, q = 3)
})
```

We plot the distributions again:

```{r}
 boxplot(list(xy2 = unlist(pHD_XY2), yx2 = unlist(pHD_YX2)), las = 1)
```
![image](https://user-images.githubusercontent.com/21171362/232584334-596080a7-b632-4fc0-b659-444262bfcc56.png)

Note that by "inflating" the Gaussians we have also necessarily inflated their pHDs, as can be seen by the change in scale:

```{r}
 boxplot(list(xy = unlist(pHD_XY), yx = unlist(pHD_YX), 
              xy2 = unlist(pHD_XY2), yx2 = unlist(pHD_YX2)), 
              col = c(rep("deepskyblue", 2), 
                      rep("salmon", 2)), 
          las = 1)
```

![image](https://user-images.githubusercontent.com/21171362/232584680-8659671e-34f8-45da-b295-d0c95f9a1fd5.png)

While quite intuitive in this example, this aspect is less easy to grasp once dealing with single cell data. Even though all cells are embedded in the same lower dimensional space, accounting for their differential proportion is not enough; some of these samples may occupy smaller or larger volumes regardless of how many cells they are comprised of. For this reason, we need to account for the **density** of samples on top of their proportion. 

The **density** and **relative proportion** can be plugged into a simple linear model that is used to estimate how far from the null (i.e. well mixed samples in the absence of batch effects) our observed pHDs are, given the observed densities and relative proportions. A precise model would estimate the joint density of any pair of samples, take a vector of fixed relative proportions, and shuffle labels for every level of the proportions vector for 50/100 permutations calculating the pHD in both directions at every permutation. However, this becomes rather costly, as it scales quadratically with the number of group levels. 

By taking 3 quantiles of joint density and 3 fixed relative proportions (such as 0.2, 0.5 and 0.8) and generating 3 * 3  = 9 sets of random label shuffling permutations, we can greatly reduce the computational cost shifting the burden on the number of cells rather than the number of samples. This approach is also amenable to cell downsampling as the density calculation still holds - although we expect the effect size to be reduced as sampling becomes more sparse. The bulk of the calculation is actually taken care of by the highly optimized `BiocNeighbors` library, and parallelization of the permutations is easily achieved through the `BiocParallel` library.

Once distributions of null pHDs are obtained for all 9 combinations of relative proportions $\pi(A',B')$ and joint densities $\rho(A',B')$ for a given number of permutations, these are plugged into a linear model as follows:

$$
\mathrm{pHD}_{null}(A', B') \sim \pi(A',B') + \rho(A',B') + \epsilon
$$

where $A'$ and $B'$ indicate null distribution sampels, and $\epsilon$ are the residuals of the model. 

Then, the batch effect (BE) is calculated as follows. We plug the observed values of $\pi(A,B)$ and $\rho(A,B)$ in the model for every pHD(A,B) and determine the distance between the distance $\delta$ between the observed pHD and the predicted null pHD in the model. Then, we look at how many standard deviations (sigmas) are contained within $\delta(A,B)$.

$$
\mathrm{BE}(A,B) = \frac{\delta(A,B)}{\sigma_{\epsilon}}
$$

where $\sigma_{\epsilon}$ is the standard deviation of the residuals. 

To obtain a p-value (H0: samples are well mixed, i.e. the observed pHD comes from a null distribution) we calculate a Z-statistic for the observed $\delta$ given the null model residuals, and its relative p-value:

$$
Z(A,B) = \frac{\delta(A,B) - \mu_{\epsilon}}{\sigma_{\epsilon}}
$$

`HDB` does this calculation in both directions: A,B and B,A, so the resulting table will contain values of `dist` ($\delta$), `sigmas` and p-values in either direction. 
