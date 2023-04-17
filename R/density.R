#' Get intervals
#' Find which cells lie within a given set of non-contiguous intervals
#'
#' @param points matrix, a cell x coordinates matrix of cell embeddings
#' @param d numeric, the dimension (column number of `points`)
#' @param l numeric, the size of the interval
#'
#' @return a data.frame where each cell is assigned to a specific interval number

getIntervals = function(points, d, l) {

  this_dim = sort(points[,d,drop=TRUE])

  keep = which(this_dim <= this_dim[1] + l)
  boxes = data.frame("cells" = names(this_dim[keep]), "interval_n" = 1)
  nextone = keep[length(keep)] + 1
  i = 1

  while(nextone <= length(this_dim)) {
    keep2 = which(this_dim >= this_dim[nextone] & this_dim <= this_dim[nextone] + l)
    boxes = rbind(boxes, data.frame("cells" = names(this_dim[keep2]),
                                    "interval_n" = i))
    nextone =  keep2[length(keep2)] + 1
    i = i + 1
  }
  return(boxes)
}

#' Greedy plaid covering
#' Estimate the hypercube side and regions that cover an n-dimensional point set
#'
#' @param points matrix, a cell x coordinates matrix of cell embeddings
#' @param k numeric, the target number of hypercubes/boxes to use for plaid covering
#' @param alpha numeric, the tolerance on the final number of boxes given as
#'     fraction of k
#' @param verbose logical, should the function print messages on the progress?
#'     Default is `FALSE`.
#'
#' @return a list containing two elements: `boxes`, a data.frame containing the
#'    assignment of each cell to a hypercube, and `hypercube_side` the length
#'    of the side of the hyperbcube.
#'
#' @details This algorithm was originally implemented by Hie et al.
#'     (Cell Systems 2017) in Python. Briefly, an initial guess on the hypercube
#'     side is `l` is made, and then hypercubic boxes are constructed so that they are
#'     contiguous if possible. If a box is empty it is discarded. Two non-adjacent
#'     boxes can exist at a distance that is not necessarily a multiple of `l`.
#'     The final number of boxes is counted. The algorithm then proceeds by
#'     binary search increasing or decreasing `l` according to the number of boxes.
#'     If the desired number of boxes `k` +/- `k * alpha` is reached, the algorithm
#'     stops. The final result is an assignment of each cell to a hypercubic
#'     region of a defined side `l`, which can be used to estimate the density as
#'     `(number of points) / l^(number of dimensions)`.
#'
#' @export

boxDensity = function(points, k = NULL, alpha = 0.1, verbose = FALSE) {

  if(is.null(k)) k = floor(sqrt(nrow(points)))

  # Translating and scaling
  points = t(t(points) - apply(points, 2, min))
  points = points/(max(points))

  # Initialization
  l = (diff(range(points)) * ncol(points))/k
  if(verbose) message(cat("Initial l: ", l, sep = ""))
  boxes =  lapply(seq_len(ncol(points)), function(x) getIntervals(points, d = x, l = l))

  for(i in seq_len(length(boxes))) {
    rownames(boxes[[i]]) = boxes[[i]]$cells
    boxes[[i]] = boxes[[i]][boxes[[1]]$cells,]
  }

  boxes = lapply(boxes, function(x) x[,2,drop=FALSE])
  all_intervals = do.call(cbind, boxes)
  all_intervals$box_name = as.numeric(factor(apply(all_intervals, 1, function(x) paste0(x, collapse = ""))))
  nboxes = length(unique(all_intervals$box_name))

  # Binary search
  while(nboxes < round(k*(1 - alpha)) | nboxes > round(k*(1 + alpha) )) { # tolerance
    if(nboxes < round(k*(1 - alpha))) {
      lmax = l
      lmin = l/2
      l = mean(c(lmax, lmin))
      if(verbose) message(cat("Too few boxes (", nboxes, "). Decreasing l to: ", l, sep = ""))

      boxes =  lapply(seq_len(ncol(points)), function(x) getIntervals(points, d = x, l = l))

      for(i in seq_len(length(boxes))) {
        rownames(boxes[[i]]) = boxes[[i]]$cells
        boxes[[i]] = boxes[[i]][boxes[[1]]$cells,]
      }
      boxes = lapply(boxes, function(x) x[,2,drop=FALSE])
      all_intervals = do.call(cbind, boxes)
      all_intervals$box_name = as.numeric(factor(apply(all_intervals, 1, function(x) paste0(x, collapse = ""))))
      nboxes = length(unique(all_intervals$box_name))
    } else if(nboxes > round(k*(1 + alpha) )) {
      lmax = l*2
      lmin = l
      l = mean(c(lmax, lmin))
      if(verbose) message(cat("Too many boxes (", nboxes, "). Increasing l to: ", l, sep = ""))
      boxes = lapply(seq_len(ncol(points)), function(x) getIntervals(points, d = x, l = l))

      for(i in seq_len(length(boxes))) {
        rownames(boxes[[i]]) = boxes[[i]]$cells
        boxes[[i]] = boxes[[i]][boxes[[1]]$cells,]
      }
      boxes = lapply(boxes, function(x) x[,2,drop=FALSE])
      all_intervals = do.call(cbind, boxes)
      all_intervals$box_name = as.numeric(factor(apply(all_intervals, 1, function(x)
        paste0(x, collapse = ""))))
      nboxes = length(unique(all_intervals$box_name))
    }
  }
  if(verbose) message(cat("Final number of boxes: ", nboxes,
                          "\nFinal hypercube side: ", l, sep = ""))

  return(list("boxes" = all_intervals[,"box_name", drop = FALSE],
              "hypercube_side" = l))
}

#' Estimate densities
#' Uses greedy plaid covering to estimate per-batch global densities
#'
#' @param sce a `SingleCellExperiment` object
#' @param dimred character, the name of the reduced dimension slot in `sce`
#' @param group character, name in the `colData` slot of `sce` containing information
#'     on the group/batch identity.
#' @param ndims numeric, the number of dimensions to consider in `dimred`. The
#'     first `ndims` will be considered (1:ndims). Default is NULL meaning all
#'     available dimensions will be used.
#' @param k numeric, the target number of hypercubes/boxes to use for plaid covering
#' @param alpha numeric, the tolerance on the final number of boxes given as
#'     fraction of k
#' @param verbose logical, should the function print messages on the progress?
#'     Default is `FALSE`.
#'
#' @return a vector containing the density for each level of `group`.
#'
#' @importFrom SingleCellExperiment reducedDim
#' @importFrom SummarizedExperiment colData
#'
#' @export

getDensities <- function(sce, dimred, group, ndims = NULL, k = NULL, alpha = 0.1, verbose = FALSE) {

  space = reducedDim(sce, dimred)
  if(is.null(ndims)) dims = colnames(dimred) else dims = seq_len(ndims)
  space = space[,dims]

  if(verbose) message("Plaid covering.\n")

  boxes = boxDensity(space, k = k, alpha = alpha, verbose = verbose)
  d = ncol(space)

  if(verbose) message(cat("Calculating density.\n"))

  densities = sapply(unique(colData(sce)[,group]), function(x) {
    if(verbose) message("Density for group level ", as.character(x))
    cells = rownames(colData(sce))[colData(sce)[,group] == as.character(x)]
    boxes_in = unique(boxes$boxes[cells, "box_name"])
    total_volume = length(boxes_in) * boxes$hypercube_side^d
    density = length(cells)/total_volume
    return(density)
  })

  if(verbose) message("All done.\n")

  names(densities) = unique(colData(sce)[,group])

  return(densities)
}

#' Estimate joint densities
#' Uses greedy plaid covering to estimate pairwise per-batch global densities
#'
#' @param sce a `SingleCellExperiment` object
#' @param dimred character, the name of the reduced dimension slot in `sce`
#' @param group character, name in the `colData` slot of `sce` containing information
#'     on the group/batch identity.
#' @param combs a data.frame containing pairwise combinations of group levels.
#'     default is NULL meaning combinations will be constructed for all levels.
#' @param ndims numeric, the number of dimensions to consider in `dimred`. The
#'     first `ndims` will be considered (1:ndims). Default is NULL meaning all
#'     available dimensions will be used.
#' @param k numeric, the target number of hypercubes/boxes to use for plaid covering
#' @param alpha numeric, the tolerance on the final number of boxes given as
#'     fraction of k
#' @param verbose logical, should the function print messages on the progress?
#'     Default is `FALSE`.
#'
#' @return a named vector containing the density for each pair of `group` levels.
#'
#' @importFrom SingleCellExperiment reducedDim
#' @importFrom SummarizedExperiment colData
#'
#' @export

getJointDensities <- function(sce, dimred, group, combs = NULL, ndims = NULL,
                              k = NULL, alpha = 0.1, verbose = FALSE) {

  space = reducedDim(sce, dimred)
  if(is.null(ndims)) dims = seq_len(ncol(space)) else dims = seq_len(ndims)
  space = space[,dims]

  if(is.null(combs)) {
    combs = data.frame(t(combn(as.character(unique(colData(sce)[,group])), m=2)))
    colnames(combs) = c("A", "B")
  }
  if(verbose) message(cat("Plaid covering.\n"))

  boxes = boxDensity(space, k = k, alpha = alpha, verbose = verbose)
  d = ncol(space)

  if(verbose) message("Calculating density.\n")

  densities = sapply(seq_len(nrow(combs)), function(x) {
    curr = as.character(combs[x,1:2])
    if(verbose) message("Joint density for groups ", paste0(curr, collapse = ", "))
    cells = rownames(colData(sce))[colData(sce)[,group] %in% curr]
    boxes_in = unique(boxes$boxes[cells, "box_name"])
    total_volume = length(boxes_in) * boxes$hypercube_side^d
    density = length(cells)/total_volume
    return(density)
  })

  if(verbose) message("All done.\n")

  names(densities) = apply(combs[,1:2], 1, function(x) paste0(x, collapse = "_"))

  return(densities)
}
