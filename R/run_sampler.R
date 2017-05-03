#' Aggregated and overdispersed sampling
#'
#' @author Bruno Vilela
#'
#' @description Generate aggregated or overdispersed sampling designs
#' for any given distance matrix (class matrix). Results can
#' be used to design experiments/samples, for resample proposes and data bias removal
#'
#' @param x \code{matrix} indicating the distance (any unit) between sample units.
#' Row and column names should be given.
#' @param n A positive integer number indicating the sample size.
#' @param alpha Number indicating the strenght of aggregation (if negative) or
#' overdispersion (if positive). When alpha = 0 sample is random.
#'
#'
#' @details \code{run_sampler} resample \code{n} sample units with an attraction
#' or repulsive effect determined by \code{alpha} and given a disatance matrix
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
#' @values The function returns a vector indicating the selected rows.
#'
#' @examples
#' # Generate a random tree
#' require(ape)
#' set.seed(2)
#' tree <- rcoal(10)
#' # Calculate the distance
#' dist <- cophenetic(tree)
#' # Highly overdispersed 50% resample design (alpha = 100)
#' selection <- run_sampler(x = dist, n = 5, alpha = 1000)
#' # Prune tree
#' overdispersed <- drop.tip(tree, tree$tip.label[!tree$tip.label %in% selection])
#' # Plot to compare
#' par(mfrow = c(1, 2))
#' plot(tree)
#' plot(overdispersed)
#'

run_sampler <- function (x, n, alpha) {

  # Error control
  .error_control(x, n, alpha)
  # Names
  row_names_x <- rownames(x)

  # Starter
  n_row <- nrow(x)
  starter <- sample(1:n_row, 1)
  first <- row_names_x[starter]
  # Sample second
  dist_rela <- x[starter, -starter] ^ alpha
  if (any(is.infinite(dist_rela))) {
    stop("alpha is too high, infinite values generated.")
  }
  prob <- dist_rela/sum(dist_rela)
  second <- sample(names(prob), 1, prob = prob)

  selected <- character(n)
  selected[1:2] <- c(first, second)
  # loop
  for (i in 3:n) {
    positions <- row_names_x %in% selected
    dist_rela <- colMeans(x[positions, !positions]) ^ alpha
    prob <- dist_rela/sum(dist_rela)
    selected[i] <- sample(names(prob), 1, prob = prob)
  }
  return(selected)
}



# Error control function
.error_control <- function (x, n, alpha) {
  if (!is.matrix(x)) {
    stop("x must be a matrix.")
  }
  if (nrow(x) != ncol(x)) {
    stop("x must be a symmetric distance matrix.")
  }
  if (is.null(colnames(x)) | is.null(rownames(x))) {
    stop("x should have both columns and rows named.")
  }
  if (!is.numeric(n)) {
    stop("n must be a number")
  }
  if (length(n) != 1) {
    stop("n should have length 1")
  }
  if (n < 2) {
    stop("n must be higher than 1")
  }
  if (n/n != 1) {
    stop("n must be an integer")
  }
  if (!is.numeric(alpha)) {
    stop("alpha should a number")
  }
  if (length(alpha) != 1) {
    stop("alpha should have length 1")
  }
  invisible(NULL)
}
