#' Robust partial Hausdorff distance from distance matrix
#'
#' Returns the  robust partial Hausdorff distance between two sets of points
#'
#' @param A a matrix with the target set of points
#' @param B a matrix with the query set of points
#' @param q numeric, the rank of the minimum distance for robust calculations
#' @param distance string specifying the distance method, one of "Euclidean",
#'     "Cosine" or "Manhattan"
#' @param BNPARAM a \code{BiocNeighborParam} class
#'
#' @importFrom BiocNeighbors KmknnParam queryKNN
#'
#' @return a numeric with the Hausdorff distance from A to B

.pHD <- function(A, B, q, distance="Euclidean",
                 BNPARAM=KmknnParam(distance=distance)){
  nn = queryKNN(A, B, k=1, get.index=FALSE, BNPARAM=BNPARAM)
  k = length(nn$distance) - q
  return(sort(nn$distance)[k])
}

#' Hausdorff Distance of Batches (old implementation)
#'
#' Calculates the pairwise batch effect using the robust partial Hausdorff
#'    distance measured against the null distribution of Hausdorff distances.
#'    The null distribution is empirically calculated by permutation of the
#'    batch label at different relative proportions.
#'
#' @param sce a SingleCellExperiment object
#' @param dimred character, the name of the dimensional reduction slot
#' @param group character or vector of characters containing the name(s) of the
#'     \code{colData} column containing the batch information
#' @param distance string, the distance method to be used by
#'    the \code{BNPARAM}. Default is "Euclidean", can be "Cosine", "Manhattan".
#' @param q numeric, the rank for robust pHD calculation
#' @param dims numeric, the number of dimensions in the reduced dimension to
#'    be used
#' @param nthreads numeric, the number of threads to be used for parallel
#'    computation. Default is \code{1}, which selects serial processing.
#' @param fit character, one of \code{"linear"} (default) and \code{"nonlinear"}
#' @param props vector of numerics, the relative proportions to be used for the
#'    null distribution estimation
#' @param samples numeric, the number of iterations of label shuffling for every
#'    proportion. Default is \code{50}.
#' @param BNPARAM a \code{BiocNeighborParam} class. Inherits \code{distance}
#'    from the function argument. If different distances are specified, will
#'    default to the one specified in the BNPARAM.
#' @param doplot logical, should the results be plotted at the end? Default is
#'    \code{TRUE}
#' @param verbose logical, should the function return updates on its progress?
#'    This includes progress bars for parallel processing. Default is FALSE.
#'
#' @importFrom SummarizedExperiment colData
#' @importFrom SingleCellExperiment reducedDimNames reducedDim
#' @importFrom BiocNeighbors KmknnParam bndistance
#' @importFrom stats lm predict residuals sd
#' @importFrom utils combn
#'
#' @return The function returns a named list of results, one for each element of
#'    \code{group}, each list containing another list with the following items:
#'  \itemize{
#'  \item{\code{model}}: a list of linear models, one for every pairwise pHD
#'      calculation within the group
#'  \item{\code{results}}: a data frame of pHD results, containing the pHD value
#'       and the number of standard deviations from the linear
#'      model (\code{sigmas})
#'  \item{\code{null}}: a list of null pHD distributions, one for every relative
#'  proportion within the group
#' }
#' @export

HDB_old <- function(sce, dimred="PCA", group, distance="Euclidean", q=3, dims=20,
                nthreads=1, fit="linear", props=c(0.2, 0.4, 0.5, 0.6, 0.8),
                samples=50, BNPARAM=KmknnParam(distance=distance),
                doplot=TRUE, verbose=FALSE) {

  # Checks
  match.arg(group, colnames(colData(sce)), several.ok=TRUE)
  match.arg(dimred, reducedDimNames(sce), several.ok=FALSE)
  if(any(props >= 1) | any(props <= 0)) {
    stop("props must be between 0 and 1 excluded")
  }

  for(g in group) {
    if(any(table(colData(sce)[,g]) < q)) {
    stop("q cannot be higher than the size of a group")
    }
  }

  if(bndistance(BNPARAM) != distance){
    message("Conflicting distance methods found - will default to the method
            specified in BNPARAM")
    distance = bndistance(BNPARAM)
  }

  space = reducedDim(sce, dimred)[, seq_len(dims), drop = FALSE]

  results <- list()

  for(g in group){

    if(verbose) message("Computing null distributions for group ", g, "\n")

    combs = data.frame(t(combn(as.character(unique(colData(sce)[,g])), m=2)))
    colnames(combs) = c("A", "B")
    nullist = list()

    if(nthreads == 1) {
      par_param = BiocParallel::SerialParam()
    } else {
      par_param = BiocParallel::MulticoreParam(workers = nthreads,
                                               progressbar = verbose)
    }

    phds_1 = phds_2 = vector()

    props_actual = c(apply(combn(table(colData(sce)[,g]), m=2), 2, function(x)
      x[1]/(x[1] + x[2])),
                     apply(combn(table(colData(sce)[,g]), m=2), 2, function(x)
                       x[2]/(x[1] + x[2])))

    if(any(props_actual < 0.2)) {

      message("Some elements of this group are highly unbalanced
              (relative proportion < 0.2). Smaller and larger proportions will
              be added to the proportion vector to estimate the null model.\n
              This can result in a suboptimal linear model fitting and inflated
              sigmas. \nYou can select \"fit = \"nonlinear\" to overcome
              this issue.")

      props = union(props, range(props_actual))
    }
    if(nrow(combs) <= 25) {

      nullist <- lapply(seq_len(nrow(combs)), function(j) {

        space_j_A = space[which(colData(sce)[,g] == combs[j, "A"]), ]
        space_j_B = space[which(colData(sce)[,g] == combs[j, "B"]), ]

        space_j = rbind(space_j_A, space_j_B)

        nulls = list()

        for(i in seq_len(length(props))){
          nulls[[i]] = BiocParallel::bplapply(seq_len(samples), function(x) {
            row_prop = floor(nrow(space_j) * props[i])
            sampled = sample(seq_len(nrow(space_j)), row_prop, replace=FALSE)
            space_s_A = space_j[sampled, ]
            space_s_B = space_j[-sampled, ]
            pHD_AB = .pHD(A=space_s_B, B=space_s_A, q=q, BNPARAM=BNPARAM)
            return(pHD_AB)
          }, BPPARAM = par_param)
        }

        return(nulls)

      })
    } else  {

      nullist <- BiocParallel::bplapply(seq_len(nrow(combs)), function(j) {

        space_j_A <- space[which(colData(sce)[,g] == combs[j, "A"]), drop=FALSE]
        space_j_B <- space[which(colData(sce)[,g] == combs[j, "B"]), drop=FALSE]

        space_j <- rbind(space_j_A, space_j_B)

        nulls <- list()

        for(i in seq_len(length(props))){
          nulls[[i]] <- lapply(seq_len(samples), function(x) {
            row_prop <- floor(nrow(space_j) * props[i])
            sampled <- sample(seq_len(nrow(space_j)), row_prop, replace = FALSE)
            space_s_A <- space_j[sampled, ]
            space_s_B <- space_j[-sampled, ]
            pHD_AB <- .pHD(A=space_s_B, B=space_s_A, q=q, BNPARAM=BNPARAM)
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

      nullmodels <- lapply(nullist_dfs, function(x)
        lm(data=x, formula=phd ~ prop))

    } else if(fit == "nonlinear") {
      nullmodels <- lapply(nullist_dfs, function(x)
        lm(data=x, formula=phd ~ poly(prop, 3)))
    }

    if(verbose) message("Computing pHds for group ", g, "\n")

    props <- vector()
    phds <- vector()

    comp_df <- expand.grid(unique(colData(sce)[,g]), unique(colData(sce)[,g]))
    colnames(comp_df) <- c("A", "B")

    for(i in seq_len(nrow(comp_df))) {
      select_A = which(colData(sce)[,g] == comp_df[i, "A"])
      select_B = which(colData(sce)[,g] == comp_df[i, "B"])
      props[i] = length(select_A) / (length(select_A) + length(select_B))

      space_j_A = space[select_A, ]
      space_j_B = space[select_B, ]

      phds[i] = .pHD(A=space_j_B, B=space_j_A, q=q, BNPARAM=BNPARAM)
    }

    dat_df = data.frame(phd=phds,
                         prop=props,
                         from=as.character(comp_df$A),
                         to=as.character(comp_df$B))

    dat_df$comps = rownames(dat_df) = paste0(dat_df$from, "_", dat_df$to)


    dists = vector()
    sigmas = vector()

    for(i in seq_len(nrow(dat_df))) {

      p1 = dat_df[i,5]
      p2 = paste0(dat_df[i,4], "_", dat_df[i,3])

      if(p1 == p2) {
        dists[i] = NA
        sigmas[i] = NA
      } else {

        model.select = nullmodels[[which(names(nullmodels) %in% c(p1, p2))]]

        dists[i] = dat_df$phd[i] - predict(newdata=data.frame(prop=dat_df$prop[i]),
                                           object=model.select)
        sigmas[i] = dists[i]/sd(residuals(model.select))

      }
    }

    dat_df$dist = dists
    dat_df$sigmas = sigmas

    results[[g]] = list(model=nullmodels, results=dat_df, null=nullist_dfs)
  }
  if(doplot) {

    plotHDsigmas(hdb=results,
                 group=g)
  }

  return(results)
}

#' Plot HDB sigmas (heatmap)
#'
#' Plots the sigma values for every HDB comparison in a heatmap
#'
#' @param hdb the result of \code{HDB()}
#' @param group character containing the name of the \code{HDB} results list
#' @param rng vector containing a range in the form \code{c(a, b)} where
#'    \code{a} and \code{b} contain the lower and upper bound of the scale.
#'    Used to compare different heatmaps using the same scale.
#'
#' @return a heatmap with pairwise sigma values for the HDB calculation within
#'    \code{group}.
#'
#' @importFrom ggplot2 ggplot geom_tile scale_fill_gradient2 geom_text
#'    theme_minimal theme ggtitle element_text element_blank aes_string unit
#'
#' @export


plotHDheatmap <- function(hdb, group, rng = NULL) {

  current = hdb[[group]]

  phdmat = matrix(current$results$sigmas,
                  nrow=length(unique(current$results$from)))
  rownames(phdmat) = unique(as.character(current$results$from))
  colnames(phdmat) = rownames(phdmat)

  M = reshape2::melt(phdmat)

  M$Var1 = as.character(M$Var1)
  M$Var2 = as.character(M$Var2)

  if(is.null(rng)) rng <- range(M$value, na.rm = TRUE)

  p = ggplot(M, aes_string(x = "Var2", y = "Var1", fill = "value")) +
    geom_tile(color = "black", width=0.95, height=0.95) +
    geom_text(aes_string(label = round(M$value, 2)), size = 6) +
    scale_fill_gradient2(low="white", mid="lightblue",
                         high="steelblue", na.value="gray",
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


#' Plot HDB sigmas (pyramid barplot)
#'
#' Plots the sigma values for every HDB comparison in a pyramid barplot
#'
#' @param hdb the result of \code{HDB()}
#' @param group character containing the name of the \code{HDB} results list
#' @param rng vector containing a range in the form \code{c(a, b)} where
#'    \code{a} and \code{b} contain the lower and upper bound of the scale.
#'   Used to compare different barplots using the same scale.
#'
#' @return a pyramid barplot with pairwise sigma values for the HDB calculation
#'    within \code{group}.
#'
#' @importFrom ggplot2 ggplot geom_hline aes_string geom_bar theme theme_minimal
#'    ggtitle coord_flip element_text scale_y_continuous xlab
#' @importFrom stats reorder
#'
#' @export

plotHDsigmas <- function(hdb, group, rng = NULL) {

  current = hdb[[group]]

  df = current$results[!is.na(current$results$dist),]

  comparisons = apply(df[,c("from", "to")], 1, function(x)
    paste(sort(x), collapse = "_"))

  df$comparison = comparisons

  df$direction = c("A_B", "B_A")[as.numeric(sapply(names(comparisons), function(x)
    comparisons[x] == x)) + 1]

  df$sigmas[df$direction == "A_B"] = -df$sigmas[df$direction == "A_B"]

  sigma_sum = sapply(unique(df$comparison), function(x)
    sum(abs(df$sigmas[df$comparison == x])))

  df$sum = sigma_sum[df$comparison]

  df$comparison = gsub(df$comparison, pattern = "_", replacement = " and ")
  df$direction = gsub(df$direction, pattern = "_", replacement = " to ")

  if(is.null(rng)) {
    rng = rep(max(range(abs(df$sigmas), na.rm=TRUE)), 2) * c(-1, 1)
  } else {
    rng = rep(max(rng), 2) * c(-1, 1)
  }

  p = ggplot(df, aes_string(x=reorder(df$comparison, df$sum), y="sigmas",
                             fill="direction")) +
    geom_bar(stat = "identity") +
    geom_hline(yintercept=0) +
    geom_hline(yintercept=c(-3, 3), lty=2, lwd=0.5) +
    theme_minimal()  +
    scale_y_continuous(breaks=round(seq(rng[1], rng[2], 10)),
                     labels=as.character(abs(round(seq(rng[1], rng[2], 10))))) +
    xlab("Comparisons (A and B)") +
    ggtitle(paste0("HDB sigmas for group ", group)) +
    coord_flip()

  print(p)
}


#' Hausdorff Distance of Batches
#'
#' Calculates the pairwise batch effect using the robust partial Hausdorff
#'    distance measured against the null distribution of Hausdorff distances.
#'    The null distribution is empirically calculated by permutation of the
#'    batch label at different relative proportions and joint densities.
#'
#' @param sce a SingleCellExperiment object
#' @param dimred character, the name of the dimensional reduction slot
#' @param group character or vector of characters containing the name(s) of the
#'     \code{colData} column containing the batch information
#' @param distance string, the distance method to be used by
#'    the \code{BNPARAM}. Default is "Euclidean", can be "Cosine", "Manhattan".
#' @param q numeric, the rank for robust pHD calculation
#' @param dims numeric, the number of dimensions in the reduced dimension to
#'    be used
#' @param props vector of numerics, the relative proportions to be used for the
#'    null distribution estimation
#' @param samples numeric, the number of iterations of label shuffling for every
#'    proportion. Default is \code{50}.
#' @param BNPARAM a \code{BiocNeighborParam} class. Inherits \code{distance}
#'    from the function argument. If different distances are specified, will
#'    default to the one specified in the BNPARAM.
#' @param BPPARAM a \code{BiocParallel} param class. Default is \code{SerialParam()}
#'    meaning no parallelization is implemented.
#' @param doplot logical, should the results be plotted at the end? Default is
#'    \code{TRUE}
#' @param verbose logical, should the function return updates on its progress?
#'    Default is FALSE.
#'
#' @importFrom SummarizedExperiment colData
#' @importFrom SingleCellExperiment reducedDimNames reducedDim
#' @importFrom BiocNeighbors KmknnParam bndistance
#' @importFrom BiocParallel SerialParam
#' @importFrom stats lm predict residuals sd pnorm rstandard
#' @importFrom utils combn
#'
#' @return The function returns a named list of results, one for each element of
#'    \code{group}, each list containing another list with the following items:
#'  \itemize{
#'  \item{\code{model}}: an object of class \code{lm} containing the linear model
#'      used to determine the batch effect sigmas and residuals
#'  \item{\code{results}}: a data frame of pHD results, containing the pHD value
#'       and the number of standard deviations from the linear
#'      model (\code{sigmas}) and the nominal p-value
#'  \item{\code{null}}: a list of null pHD distributions, one for every joint density
#'  within the group
#' }
#' @export

HDB <- function(sce, dimred="PCA", group, distance="Euclidean", q=3, dims=20,
                 props=c(0.2, 0.5, 0.8),
                 samples=50, BNPARAM=KmknnParam(distance=distance), BPPARAM=SerialParam(),
                 doplot=TRUE, verbose=FALSE) {

  # Checks
  match.arg(group, colnames(colData(sce)), several.ok=TRUE)
  match.arg(dimred, reducedDimNames(sce), several.ok=FALSE)

  if(any(props >= 1) | any(props <= 0)) {
    stop("props must be between 0 and 1 excluded")
  }

  for(g in group) {
    if(any(table(colData(sce)[,g]) < q)) {
      stop("q cannot be higher than the size of a group")
    }
  }

  if(bndistance(BNPARAM) != distance){
    message("Conflicting distance methods found - will default to the method
            specified in BNPARAM")
    distance = bndistance(BNPARAM)
  }

  space = reducedDim(sce, dimred)[, seq_len(dims), drop = FALSE]

  results <- list()

  for(g in group){

    if(verbose) message("Computing null distributions for group ", g, "\n")

    combs = data.frame(t(combn(as.character(unique(colData(sce)[,g])), m=2)))
    colnames(combs) = c("A", "B")

    if(verbose) message("Calculating joint densities.\n")
    dens = getJointDensities(sce,
                             dimred = dimred,
                             group = g,
                             combs = combs, verbose = verbose)
    combs$dens = dens
    combs = combs[order(combs$dens),]

    if(nrow(combs) > 3) {
      dens_select = unique(c(2, round(nrow(combs)/2)+1, nrow(combs)))
    } else {
      dens_select = seq_len(nrow(combs))
    }

    combs$props_AB = apply(combn(table(colData(sce)[,g]), m=2), 2, function(x)
      x[1]/(x[1] + x[2]))
    combs$props_BA = apply(combn(table(colData(sce)[,g]), m=2), 2, function(x)
        x[2]/(x[1] + x[2]))

    null_combs = expand.grid(combs$dens, props)
    null_combs$A = rep(combs$A, table(as.numeric(null_combs[,1])))
    null_combs$B = rep(combs$B, table(as.numeric(null_combs[,1])))

    nullist = lapply(seq_len(nrow(null_combs)), function(j) {

      space_j_A = space[which(colData(sce)[,g] == null_combs[j, "A"]), ]
      space_j_B = space[which(colData(sce)[,g] == null_combs[j, "B"]), ]

      space_j = rbind(space_j_A, space_j_B)

      nulls = list()

      for(i in seq_len(length(props))){
        nulls[[i]] = BiocParallel::bplapply(seq_len(samples), function(x) {
          row_prop = floor(nrow(space_j) * props[i])
          sampled = sample(seq_len(nrow(space_j)), row_prop, replace=FALSE)
          space_s_A = space_j[sampled, ]
          space_s_B = space_j[-sampled, ]
          pHD_AB = .pHD(A=space_s_B, B=space_s_A, q=q, BNPARAM=BNPARAM)
          return(pHD_AB)
        }, BPPARAM = BPPARAM)
      }
      return(nulls)
    })

   names(nullist) <- paste0(combs$A, "_", combs$B)

    nullist_dfs <- lapply(seq_len(length(dens_select)), function(x) {
      df = data.frame("phd" = unlist(nullist[[x]]))
      df$prop = rep(props, each = samples)
      df$dens = rep(combs$dens[dens_select[x]], each = samples)
      return(df)
    })

    null_df = do.call(rbind, nullist_dfs)

    if(verbose) message("Computing null models for group ", g, "\n")

    nullmodel <- lm(data=null_df, formula = phd ~ prop + dens)

    if(verbose) message("Computing pHds for group ", g, "\n")

    props_actual <- vector()
    phds <- vector()

    comp_df <- expand.grid(unique(colData(sce)[,g]), unique(colData(sce)[,g]))
    colnames(comp_df) <- c("A", "B")

    for(i in seq_len(nrow(comp_df))) {
      select_A = which(colData(sce)[,g] == comp_df[i, "A"])
      select_B = which(colData(sce)[,g] == comp_df[i, "B"])
      props_actual[i] = length(select_A) / (length(select_A) + length(select_B))

      space_j_A = space[select_A, ]
      space_j_B = space[select_B, ]

      phds[i] = .pHD(A=space_j_B, B=space_j_A, q=q, BNPARAM=BNPARAM)
    }

    dat_df = data.frame(phd=phds,
                        prop=props_actual,
                        from=as.character(comp_df$A),
                        to=as.character(comp_df$B))

    dat_df$comps = rownames(dat_df) = paste0(dat_df$from, "_", dat_df$to)

    dat_df$dens = sapply(seq_len(nrow(dat_df)), function(x)
      combs$dens[combs$A %in% dat_df[x, c("from", "to")] & combs$B %in% dat_df[x, c("from", "to")]])
    dat_df$dens[dat_df$phd == 0] <- NA
    dat_df$dens <- unlist(dat_df$dens)

    dists = vector()
    sigmas = vector()
    pvals = vector()

    for(i in seq_len(nrow(dat_df))) {

      p1 = dat_df[i,5]
      p2 = paste0(dat_df[i,4], "_", dat_df[i,3])

      if(p1 == p2) {
        dists[i] = NA
        sigmas[i] = NA
        pvals[i] = NA
      } else {


        dists[i] = dat_df$phd[i] - predict(newdata=data.frame(prop=dat_df$prop[i], dens = dat_df$dens[i]),
                                           object=nullmodel)

        sigmas[i] = dists[i]/sd(residuals(nullmodel))

        z = (dists[i] - mean(residuals(nullmodel)))/sd(residuals(nullmodel))

        pvals[i] = pnorm(z, mean = mean(rstandard(nullmodel)),
                         sd = sd(rstandard(nullmodel)), lower.tail = FALSE)

      }
    }

    dat_df$dist = dists
    dat_df$sigmas = sigmas
    dat_df$pvalue = pvals

    results[[g]] = list(model=nullmodel, results=dat_df, null=nullist_dfs)
  }

  if(doplot) {

    plotHDsigmas(hdb=results,
                 group=g)
  }

  return(results)
}

