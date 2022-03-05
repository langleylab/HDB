#' Robust partial Hausdorff distance from distance matrix
#'
#' Returns the  robust partial Hausdorff distance between two sets of points
#'
#' @param A a matrix with the target set of points
#' @param B a matrix with the query set of points
#' @param q numeric, the rank of the minimum distance for robust calculations
#'
#' @return a numeric with the Hausdorff distance from A to B

.pHD <- function(A, B, q){
  nn <- BiocNeighbors::queryKNN(A, B, k = 1, get.index = FALSE)
  k = length(nn$distance) - q
  return(sort(nn$distance)[k])
}


#' Hausdorff Distance of Batches
#'
#' Calculates the pairwise batch effect using the robust partial Hausdorff distance
#'   measured against the null distribution of Hausdorff distances. The null distribution
#'   is empirically calculated by permutation of the batch label at different
#'   relative proportions.
#'
#' @param sce a SingleCellExperiment object
#' @param dimred character, the name of the dimensional reduction slot
#' @param group character or vector of characters containing the name(s) of the \code{colData} column containing the batch information
#' @param distmethod character, the distance method to be used by \code{parDist}. Default is "euclidean", can be "cosine", "manhattan".
#' @param q numeric, the rank for robust pHD calculation
#' @param dims numeric, the number of dimensions in the reduced dimension to be used
#' @param nthreads numeric, the number of threads to be used for parallel computation. Default is \code{1}, which selects serial processing.
#' @param fit character, one of \code{"linear"} (default) and \code{"nonlinear"}
#' @param props vector of numerics, the relative proportions to be used for the null distribution estimation
#' @param samples numeric, the number of iterations of label shuffling for every proportion. Default is \code{50}.
#' @param doplot logical, should the results be plotted at the end? Default is \code{TRUE}
#' @param verbose logical, should the function return updates on its progress? This includes progress bars for parallel processing. Default is FALSE.
#'
#' @importFrom SummarizedExperiment colData
#' @importFrom SingleCellExperiment reducedDimNames reducedDim
#' @importFrom stats lm predict residuals sd
#' @importFrom utils combn
#'
#' @return a named list of results, one for each element of \code{group}, each list containing another list with the following items:
#'  * \code{model}: a list of linear models, one for every pairwise pHD calculation within the group
#'  * \code{results}: a data frame of pHD results, containing the pHD value and the number of standard deviations from the linear model (\code{sigmas})
#'  * \code{null}: a list of null pHD distributions, one for every relative proportion within the group
#'
#' @export

HDB <- function(sce,
                dimred = "PCA",
                group,
                distmethod = "euclidean",
                q = 3,
                dims = 20,
                nthreads = 1,
                fit = "linear",
                props = c(0.2, 0.4, 0.5, 0.6, 0.8),
                samples = 50,
                doplot = TRUE,
                verbose = FALSE) {

  # Checks
  match.arg(group, colnames(colData(sce)), several.ok = TRUE)
  match.arg(dimred, reducedDimNames(sce), several.ok = FALSE)
  if(any(props >= 1) | any(props <= 0)) stop("props must be between 0 and 1 excluded")

  for(g in group) if(any(table(colData(sce)[,g]) < q)) stop("q cannot be higher than the size of a group")

  space <- reducedDim(sce, dimred)[, seq_len(dims), drop = FALSE]

  results <- list()


  for(g in group){

    if(verbose) message("Computing null distributions for group ", g, "\n")

    combs <- data.frame(t(combn(as.character(unique(colData(sce)[,g])), m = 2)))
    colnames(combs) <- c("A", "B")
    nullist <- list()

    if(nthreads == 1) {
      par_param = BiocParallel::SerialParam()
    } else {
      par_param = BiocParallel::MulticoreParam(workers = nthreads, progressbar = verbose)
    }

    phds_1 = phds_2 = vector()

    props_actual = c(apply(combn(table(colData(sce)[,g]), m = 2), 2, function(x) x[1]/(x[1] + x[2])),
                     apply(combn(table(colData(sce)[,g]), m = 2), 2, function(x) x[2]/(x[1] + x[2])))

    if(any(props_actual < 0.2)) {

      message("Some elements of this group are highly unbalanced (relative proportion < 0.2). Smaller and larger proportions will be added to the proportion vector to estimate the null model.\n This can result in a suboptimal linear model fitting and inflated sigmas. \nYou can select \"fit = \"nonlinear\" to overcome this issue.")

      props = union(props, range(props_actual))
    }
    if(nrow(combs) <= 25) {

      nullist <- lapply(seq_len(nrow(combs)), function(j) {

        space_j_A <- space[which(colData(sce)[,g] == combs[j, "A"]), ]
        space_j_B <- space[which(colData(sce)[,g] == combs[j, "B"]), ]

        space_j <- rbind(space_j_A, space_j_B)

        nulls <- list()

        for(i in seq_len(length(props))){
          nulls[[i]] <- BiocParallel::bplapply(seq_len(samples), function(x) {
            row_prop <- floor(nrow(space_j) * props[i])
            sampled <- sample(seq_len(nrow(space_j)), row_prop, replace = FALSE)
            space_s_A <- space_j[sampled, ]
            space_s_B <- space_j[-sampled, ]
            pHD_AB <- .pHD(space_s_B, space_s_A, q = q)
            return(pHD_AB)
          }, BPPARAM = par_param)
        }

        return(nulls)

      })
    } else  {

      nullist <- BiocParallel::bplapply(seq_len(nrow(combs)), function(j) {

        space_j_A <- space[which(colData(sce)[,g] == combs[j, "A"]), ]
        space_j_B <- space[which(colData(sce)[,g] == combs[j, "B"]), ]

        space_j <- rbind(space_j_A, space_j_B)

        nulls <- list()

        for(i in seq_len(length(props))){
          nulls[[i]] <- lapply(seq_len(samples), function(x) {
            row_prop <- floor(nrow(space_j) * props[i])
            sampled <- sample(seq_len(nrow(space_j)), row_prop, replace = FALSE)
            space_s_A <- space_j[sampled, ]
            space_s_B <- space_j[-sampled, ]
            pHD_AB <- .pHD(space_s_B, space_s_A, q = q)
            return(pHD_AB)
          })
        }

        return(nulls)

      }, BPPARAM = par_param)
    }

    names(nullist) <- paste0(combs$A, "_", combs$B)

    nullist_dfs <- lapply(nullist, function(x) {
      df = data.frame("phd" = unlist(x))
      df$prop = rep(props, each = samples)
      return(df)
    })

    if(verbose) message("Computing null models for group ", g, "\n")

    if(fit == "linear") {

      nullmodels <- lapply(nullist_dfs, function(x) lm(data = x, formula = phd ~ prop))

    } else if(fit == "nonlinear") {
      nullmodels <- lapply(nullist_dfs, function(x) lm(data = x, formula = phd ~ poly(prop, 3)))
    }

    if(verbose) message("Computing pHds for group ", g, "\n")

    props <- vector()
    phds <- vector()

    comp_df <- expand.grid(unique(colData(sce)[,g]), unique(colData(sce)[,g]))
    colnames(comp_df) <- c("A", "B")

    for(i in 1:nrow(comp_df)) {
      select_A = which(colData(sce)[,g] == comp_df[i, "A"])
      select_B = which(colData(sce)[,g] == comp_df[i, "B"])
      props[i] <- length(select_A) / (length(select_A) + length(select_B))

      space_j_A <- space[select_A, ]
      space_j_B <- space[select_B, ]

      phds[i] <- .pHD(space_j_B, space_j_A, q = q)
    }

    dat_df <- data.frame(phd = phds,
                         prop = props,
                         from = as.character(comp_df$A),
                         to = as.character(comp_df$B))

    dat_df$comps <- rownames(dat_df) <- paste0(dat_df$from, "_", dat_df$to)


    dists = vector()
    sigmas = vector()

    for(i in seq_len(nrow(dat_df))) {

      p1 = dat_df[i,5]
      p2 = paste0(dat_df[i,4], "_", dat_df[i,3])

      if(p1 == p2) {
        dists[i] <- NA
        sigmas[i] <- NA
      } else {

        model.select = nullmodels[[which(names(nullmodels) %in% c(p1, p2))]]

        dists[i] <- dat_df$phd[i] - predict(newdata = data.frame(prop = dat_df$prop[i]), object = model.select)
        sigmas[i] <- dists[i]/sd(residuals(model.select))

      }
    }

    dat_df$dist <- dists
    dat_df$sigmas <- sigmas

    results[[g]] <- list(model = nullmodels, results = dat_df, null = nullist_dfs)
  }
  if(doplot) {

    plotHdSigmas(hdb = results,
                 group = g)
  }

  return(results)
}

#' Plot HDB sigmas
#'
#' Plots the sigma values for every HDB comparison in a heatmap
#'
#' @param hdb the result of \code{HDB()}
#' @param group character containing the name of the \code{HDB} results list
#' @param rng vector containing a range in the form \code{c(a, b)} where \code{a} and \code{b} contain the lower and upper bound of the scale.
#'   Used to compare different heatmaps using the same scale.
#'
#' @return a heatmap with pairwise sigma values for the HDB calculation within \code{group}.
#'
#' @importFrom ggplot2 ggplot geom_tile scale_fill_gradient2 geom_text theme_minimal theme ggtitle element_text element_blank aes_string unit
#'
#' @export


plotHdSigmas <- function(hdb, group, rng = NULL) {

  current = hdb[[group]]

  phdmat <- matrix(current$results$sigmas, nrow = length(unique(current$results$from)))
  rownames(phdmat) <- unique(as.character(current$results$from))
  colnames(phdmat) <- rownames(phdmat)

  M <- reshape2::melt(phdmat)

  M$Var1 <- as.character(M$Var1)
  M$Var2 <- as.character(M$Var2)

  if(is.null(rng)) rng <- range(M$value, na.rm = TRUE)

  p <- ggplot(M, aes_string(x = "Var2", y = "Var1", fill = "value")) +
    geom_tile(color = "black", width=0.95, height=0.95) +
    geom_text(aes_string(label = round(M$value, 2)), size = 6) +
    scale_fill_gradient2(low="white", mid="lightblue", high="steelblue", na.value = "gray",
                         midpoint = mean(0:max(rng)),
                         breaks = seq(0, max(rng), floor(max(rng)/5)),
                         limits = c(0, max(rng))) +
    theme_minimal() +
    theme(axis.text = element_text(size = 18),
          axis.title = element_blank(),
          legend.key.size = unit(2, 'cm'),
          legend.title = element_text(size = 20),
          legend.text = element_text(size=18),
          plot.title = element_text(size=22)) +
    ggtitle("HDB sigmas")

  print(p)

}
