#' Aggregated and overdispersed sampling for geographic coordinates
#'
#' @author Bruno Vilela
#'
#' @description Generate spatially aggregated or overdispersed subsamples from
#' any given set of coordinates. Results can be used to design sampling schemes
#' for future research, for resampling proposes, and for removing spatial bias
#' from data.
#'
#' @param x \code{matrix} or \code{data.frame} indicating the coordinates
#' (first column = longitude; second column = latitude). Row names should be given.
#' @param n A positive integer number indicating the size of returned sample.
#' @param alpha Number indicating the strength of aggregation (if negative) or
#' overdispersion (if positive). When alpha = 0 sample is random.
#' There are no limits to alpha, but combinations of big numbers and big
#' distances may result in an error depending on your R configurations.
#' @param dist.func A distance function to calculate coordinates distance.
#' Default is \code{\link{rdist.earth}} from package \code{fields}.
#' @param n_start Number of initial selected sample units. Default is one starting sample unit.
#' @param return_start if \code{TRUE} the starting sample units(s) is(are) returned.
#' @param starting Character vector indicating the starting sample unit(s)
#' (= to row names). If not provided, starting sample unit(s) is(are) randomly
#' selected.
#'
#' @details The function uses the algorithm in \code{\link{run_sampler}}, but
#' accepts a two column matrix of coordinates as input.
#'
#'
#' @return The function returns a subset of the original matrix of coordinates.
#' If return_start is TRUE, a list is returned with the first
#' element being the Sampling_selection - subset matrix based on the selected
#' sampling units - and Starting_points - selected starting point(s).
#' @seealso \code{\link{run_sampler}}
#'
#' @examples
#' require(sp)
#' require(maptools)
#' data(wrld_simpl)  # World map
#' Brazil <- wrld_simpl[wrld_simpl$NAME == "Brazil", ]  # Brazil (polygon)
#' coords <- slot(spsample(Brazil, 100, "regular"), "coords")
#' rownames(coords) <- paste0("t", 1:nrow(coords))
#' ## Subsample 50%
#' ### Overdispersed
#' selection.geo <- run_sampler_geo(x = coords, n = 10, alpha = 100, starting = "t10")
#' ### Aggregated
#' selection.geo2 <- run_sampler_geo(x = coords, n = 10, alpha = -100, starting = "t10")
#' ### Random
#' selection.geo3 <- run_sampler_geo(x = coords, n = 10, alpha = 0, starting = "t10")
#'
#' ## Plot
#' par(mfrow = c(1, 3), mar = c(1, 1, 15, 1))
#' plot(Brazil, main = "Overdispersed 50% sampling (red were selected)")
#' points(selection.geo, cex = 2, pch = 19, col = "red")
#' plot(Brazil, main = "Aggregated 50% sampling (blue were selected)")
#' points(selection.geo2, cex = 2, pch = 19, col = "blue")
#' plot(Brazil, main = "Random 50% sampling (green were selected)")
#' points(selection.geo3, cex = 2, pch = 19, col = "green")


run_sampler_geo <- function (x, n, alpha, dist.func = rdist.earth,
                             n_start = 1, return_start = FALSE,
                             starting = NULL) {
  #
  if (class(x) != 'matrix' & class(x) != 'data.frame' & ncol(x) != 2) {
    stop("x must be a a matrix or data.frame with 2 columns")
  }
  if (class(dist.func) != "function") {
    stop("dist.func must be a function")
  }
  coords <- x
  x <- dist.func(coords)
  .error_control2(x, rownames(coords))

  selected <- run_sampler(x, n, alpha, n_start, return_start,
                          starting)
  if (return_start) {
    selected[[1]] <- coords[rownames(coords) %in% selected[[1]], ]
    return(selected)
  } else {
    selected <- coords[rownames(coords) %in% selected, ]
    return(selected)
  }
}


