#' Aggregated and overdispersed sampling for geographic coordinates
#'
#' @author Bruno Vilela
#'
#' @description Generate aggregated or overdispersed sampling designs
#' for any given coordinates. Results can
#' be used to design experiments/samples, for resample proposes and data bias removal.
#'
#' @param x \code{matrix} or \code{data.frame} indicating the coordinates
#' (first column = longitude; second column = latitude). Row names should be given.
#' @param n A positive integer number indicating the sample size.
#' @param alpha Number indicating the strenght of aggregation (if negative) or
#' overdispersion (if positive). When alpha = 0 sample is random.
#' @param dist.func A distance function to calculate coordinates distance.
#' Default is \code{fields::rdist.earth}.
#' @param n_start Number of initial selected points. Default is one starting point.
#' @param return_start if \code{TRUE} the starting point is returned.
#' @param starting Character vector indicating the starting point (= to row names).
#' If not provided, random starting value(s) is(are) selected.
#'
#' @details \code{run_sampler} resample \code{n} sample units with an attraction
#' or repulsive effect determined by \code{alpha} and given a distance matrix
#' (\code{x}). The algorithim begins selecting one random starting point \code{i}.
#' The following sample unit is then selected based on the a probability given
#' by the distance of \code{i} to each remaing units raised to the power of
#' \code{alpha} (pr(j | i) = dij ^ alpha). The following selections will then use
#' the average distance of the remaing units to the selected ones. The procedure
#' is repeated until the selected points reach \code{n}. Positive values of
#' \code{alpha} generate overdispersed sample designs, as sample units disntant from
#' the selected unit(s) will have a higher probability of being selected. Inverselly,
#' negative values will generate an aggregated design. Note that as \code{alpha}
#' approximate the infinity (+ or -), the sample design becomes more deterministic.
#'
#'
#' @return The function returns a vector indicating the selected rows.
#' If return_start is TRUE, a list is returned with the first element being the
#' Sampling_selection - selected sampling units - and
#' Starting_points - selected starting point(s).
#'
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


