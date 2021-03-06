#' Aggregated and overdispersed sampling
#'
#' @author Bruno Vilela
#'
#' @description Generates aggregated (underdispersed) or overdispersed samples
#' of row names from any given distance matrix.
#'
#' @param x \code{matrix} indicating the distance (any unit) between sample units.
#' Row and column names should be given.
#' @param n A positive integer number indicating the size of returned sample.
#' @param alpha Number indicating the strength of aggregation (if negative) or
#' overdispersion (if positive). When alpha = 0 sample is random.
#' There are no limits to alpha, but combinations of big numbers and big
#' distances may result in an error depending on your R configurations.
#' @param n_start Number of initial selected sample units. Default is one starting sample unit.
#' @param return_start if \code{TRUE} the starting sample unit(s) is(are) returned.
#' @param starting Character vector indicating the row name of the starting
#' sample unit. If not provided starting sample unit(s) is(are) randomly selected.
#'
#'
#' @details Given a distance matrix (\code{x}), \code{run_sampler} resamples \code{n}
#' sample units with an attraction (positive) or repulsive (negative)
#' effect determined by \code{alpha}(\eqn{\alpha}).
#' The algorithm begins selecting one random starting sample unit \code{i}.
#' The following sample unit is then selected based on the probability given
#' by the distance of \code{i} to each remaining units raised to the power of
#' \code{alpha}  (\eqn{pr(j | i) = d_{i,j} ^ \alpha}). ). Each following
#' selection will then use a joint probability. The first is calculated as the
#' average distance\code{d} of all remaining units \code{j} to all selected
#' units \code{k} (\eqn{pr1(j | k) = d_{k,j} ^ \alpha}).
#' The second as the minimum distance \code{dmin} of the remaining units to the
#' selected units (\eqn{pr2(j | k) = dmin_{k,j} ^ \alpha}).
#' The second probability prevents overdispersed samples from selecting only
#' sample units located at edges. The procedure is repeated until the selected
#' sample units reach \code{n}. Positive values of \code{alpha} generate overdispersed sample
#' designs, as sample units distant from the selected unit(s) will have a higher
#' probability of being selected. Inversely, negative values will generate an
#' aggregated design. Note that as \code{alpha} approximate the infinity (+ or -), the
#' sample design becomes more deterministic.
#'
#'
#' @return The function returns a vector of row names associated with selected
#' sampling units. If return_start is \code{TRUE}, a list is returned with the first
#' element being the Sampling_selection - selected sampling units - and
#' Starting_points - selected starting point(s).
#'
#' @seealso \code{\link{run_sampler_phy}}
#' @seealso \code{\link{run_sampler_geo}}
#'
#' @examples
#' # Phylogeny example:
#' ## Generate a random tree
#' require(ape)
#' tree <- rcoal(10)
#' ## Calculate the distance
#' dist <- cophenetic(tree)
#' ## Highly overdispersed 50% resample design (alpha = 50)
#' selection <- run_sampler(x = dist, n = 5, alpha = 50, starting = "t10")
#' ## Highly aggregated 50% resample design (alpha = -50)
#' selection2 <- run_sampler(x = dist, n = 5, alpha = -50, starting = "t10")
#' ## Random 50% resample design (alpha = 0)
#' selection3 <- run_sampler(x = dist, n = 5, alpha = 0, starting = "t10")
#' ## Plot to compare
#' par(mfrow = c(1, 3))
#' plot(tree,tip.color=ifelse(tree$tip.label %in% selection, "red","black"),
#'  main = "Overdispersed 50% sampling (red were selected)", cex = 1)
#'  axis(1)
#' plot(tree,tip.color=ifelse(tree$tip.label %in% selection2, "blue","black"),
#' main = "Aggregated 50% sampling (blue were selected)", cex = 1)
#' axis(1)
#' plot(tree,tip.color=ifelse(tree$tip.label %in% selection3, "green","black"),
#' main = "Random 50% sampling (green were selected)", cex = 1)
#' axis(1)
#'
#' # Geography example
#' require(sp)
#' require(maptools)
#' require(fields)
#' data(wrld_simpl)  # World map
#' Brazil <- wrld_simpl[wrld_simpl$NAME == "Brazil", ]  # Brazil (polygon)
#' coords <- slot(spsample(Brazil, 100, "regular"), "coords")
#' rownames(coords) <- paste0("t", 1:nrow(coords))
#' ## Calculate the geographic distance
#' dist.geo <- rdist.earth(coords)
#' ## Subsample 50%
#' ### Overdispersed
#' selection.geo <- run_sampler(x = dist.geo, n = 25, alpha = 100, starting = "t10")
#' ### Aggregated
#' selection.geo2 <- run_sampler(x = dist.geo, n = 25, alpha = -100, starting = "t10")
#' ### Random
#' selection.geo3 <- run_sampler(x = dist.geo, n = 25, alpha = 0, starting = "t10")
#'
#' ## Plot
#' par(mfrow = c(1, 3), mar = c(1, 1, 15, 1))
#' plot(Brazil, main = "Overdispersed 50% sampling (red were selected)")
#' points(coords, cex = 2, pch = 19,
#' col = ifelse(rownames(coords) %in% selection.geo, "red","gray"))
#' plot(Brazil, main = "Aggregated 50% sampling (blue were selected)")
#' points(coords, cex = 2, pch = 19,
#' col = ifelse(rownames(coords) %in% selection.geo2, "blue","gray"))
#' plot(Brazil, main = "Random 50% sampling (green were selected)")
#' points(coords, cex = 2, pch = 19,
#' col = ifelse(rownames(coords) %in% selection.geo3, "green","gray"))
#'
#' # Trait example
#' ## Fake body size
#' set.seed <- 1
#' body_size <- runif(1000)
#' # Biased sample towards large species
#' set.seed <- 1
#' body_size_bias <- sample(body_size, 500, prob = body_size)
#' par(mfrow = c(1, 3))
#' hist(body_size, main = "Species body size distribution\n(n = 1000)",  xlab = "Body size")
#' hist(body_size_bias, main = "Biased samplig towards larger species\n(n = 500)",
#' xlab = "Body size")
#' # Use sampler to reduce the bias
#' dist_bs <- as.matrix(dist(body_size_bias))
#' rownames(dist_bs) <- colnames(dist_bs) <- 1:length(body_size_bias)
#' selection.bs <- run_sampler(x = dist_bs, n = 100, alpha = 100)
#' hist(body_size_bias[as.numeric(selection.bs)],
#'  main = "Overdispersed sampling of biased information \n(n = 100)",
#'  xlab = "Body size")
#'
#'
#' # Real time simulation
#' require(raster)
#' par(mfrow = c(1, 1))
#' r <- raster(res = 25)
#' values(r) <- runif(ncell(r))
#' plot(r)
#' coords <- xyFromCell(r, 1:ncell(r))
#' rownames(coords) <- 1:ncell(r)
#' dist.geo <- as.matrix(dist(coords))
#' startingI <- c(1)
#' # Change alpha and see how it works
#' for(i in (length(startingI)+1):30) {
#' selection.geo <- run_sampler(x = dist.geo, n = i, alpha = 100,
#' starting = startingI)
#' startingI <- as.numeric(selection.geo)
#' r2 <- r
#' values(r2)[as.numeric(selection.geo)] <- 1
#' values(r2)[-as.numeric(selection.geo)] <- NA
#' plot(r2, col = "gray")
#' Sys.sleep(time = .2)
#' }

run_sampler <- function (x, n, alpha, n_start = 1, return_start = FALSE,
                         starting = NULL) {

  # Error control
  .error_control(x, n, alpha, n_start, starting, return_start)
  x <- (x - min(x)) / (max(x) - min(x))
  # Names
  row_names_x <- rownames(x)

  # Starter
  n_row <- nrow(x)
  if (any(is.null(starting))) {
    starter <- sample(1:n_row, n_start)
    first <- row_names_x[starter]
  } else {
    if (length(starting) != n_start) {
      warning("Length of starting and n_start are different,
              setting n_start to be equal to length(starting)")
    }
    n_start <- length(starting)
    if (n_start >= n) {
      stop('starting length should be smaller than n')
    }
    first <- starting
    starter <- which(row_names_x %in% starting)
  }
  selected <- character(n)

  if (n_start == 1) {
    # Sample second
    dist_rela <- x[starter, -starter] ^ alpha
    if (any(is.infinite(dist_rela))) {
      stop("alpha is too high, infinite values generated.")
    }
    prob <- dist_rela/sum(dist_rela)
    second <- sample(names(prob), 1, prob = prob)
    selected[1:2] <- c(first, second)
    start_loop <- 3
  } else {
    start_loop <- n_start + 1
    selected[1:n_start] <- first
  }
  if (start_loop <= n) {
    # loop
    for (i in start_loop:n) {
      positions <- row_names_x %in% selected
      dist_rela <- apply(x[positions, !positions], 2, min) ^ alpha
      if (any(is.infinite(dist_rela) | is.na(dist_rela))) {
        stop("alpha is too high, infinite values generated.")
      }
      prob1 <- dist_rela / sum(dist_rela)
      dist_rela2 <- apply(x[positions, !positions], 2, min) ^ alpha
      if (any(is.infinite(dist_rela2)  | is.na(dist_rela))) {
        stop("alpha is too high, infinite values generated.")
      }
      prob2 <- dist_rela2 / sum(dist_rela2)
      prob <- (prob2 + prob1) / 2
      if (any(is.na(prob))) {
        stop("alpha is too high, try smaller values.")
      }
      selected[i] <- sample(names(prob), 1, prob = prob)
    }
  }
  if (return_start) {
    return(list('Sampling_selection' = selected, "Starting_points" = first))
  } else {
    return(selected)
  }
}



# Error control function
.error_control <- function (x, n, alpha, n_start, starting, return_start) {
  if (!is.matrix(x)) {
    stop("x must be a matrix.")
  }
  if (nrow(x) != ncol(x)) {
    stop("x must be a symmetric distance matrix.")
  }
  if (is.null(colnames(x)) | is.null(rownames(x))) {
    stop("x should have both columns and rows named.")
  }
  if (any(duplicated(rownames(x)))) {
    stop("duplicated names in x.")
  }
  if (any(is.na(x))) {
    stop("x contains NA.")
  }
  if (any(is.infinite(x))) {
    stop("x contains Inf.")
  }
  if (any(x < 0)) {
    stop("x contains negative numbers.")
  }
  if (!is.numeric(n)) {
    stop("n must be a number")
  }
  if (length(n) != 1) {
    stop("n should have length 1")
  }
  if (n < 2) {
    stop("n must be 2 or high")
  }
  if (n > nrow(x)) {
    stop("n should be equal or smaller than x")
  }
  if (n != round(n)) {
    stop("n must be an integer")
  }
  if (!is.numeric(alpha)) {
    stop("alpha should a number")
  }
  if (length(alpha) != 1) {
    stop("alpha should have length 1")
  }
  if (!is.numeric(n_start)) {
    stop("n_start must be a number")
  }
  if (length(n_start) != 1) {
    stop("n_start should have length 1")
  }
  if (n_start < 1) {
    stop("n_start must be a positive value equals or bigger than 1")
  }
  if (n_start != round(n_start)) {
    stop("n_start must be an integer")
  }
  if (!is.logical(return_start)){
    stop("return_start must be logical")
  }
  if (n_start >= n) {
    stop("n_start should be smaller than n")
  }
  if (!any(is.null(starting))) {
    if (any(!starting %in% colnames(x))) {
      stop("starting was not found in x")
    }
  }
  invisible(NULL)
}
