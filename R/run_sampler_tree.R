#' Aggregated and overdispersed sampling for phylo objects
#'
#' @author Bruno Vilela
#'
#' @description Generate aggregated or overdispersed sampling designs
#' for any given phylogenetic tree (class phylo). Results can
#' be used to design experiments/samples, for resample proposes and data bias removal.
#'
#' @param x \code{phylo} object.
#' @param n A positive integer number indicating the sample size.
#' @param alpha Number indicating the strenght of aggregation (if negative) or
#' overdispersion (if positive). When alpha = 0 sample is random.
#' @param dist.func Function to calculate the phylogenetic distance.
#' The function should be able to receive a phylo object and return
#' a \code{matrix} indicating the pairwise distance
#' between tips. Row and column names should be given.
#' Default uses \code{cophenetic} from \code{ape} package.
#' Other examples include
#' @param n_start Number of initial selected tips. Default is one starting tip.
#' @param return_start if \code{TRUE} the starting tip is returned.
#' @param starting Character vector indicating the starting tips. If not provided,
#' random starting value(s) is(are) selected.
#'
#' @details \code{run_sampler} resample \code{n} sample tips with an attraction
#' or repulsive effect determined by \code{alpha} (\code{x}).
#' The algorithim begins selecting one random starting tip \code{i}.
#' The following tip is then selected based on the a probability given
#' by the distance of \code{i} to each remaing tip raised to the power of
#' \code{alpha} (pr(j | i) = dij ^ alpha). The following selections will then use
#' the average distance of the remaing tips to the selected ones. The procedure
#' is repeated until the selected tips reach \code{n}. Positive values of
#' \code{alpha} generate overdispersed sample designs, as tips within higher distant from
#' the selected tip(s) will have a higher probability of being selected. Inverselly,
#' negative values will generate an aggregated design. Note that as \code{alpha}
#' approximate the infinity (+ or -), the sample design becomes more deterministic.
#'
#'
#' @values The function returns a prunned phylogenetic tree.
#'
#' @examples
#' # Generate a random tree
#' require(ape)
#' set.seed(100)
#' tree <- rcoal(10)
#' set.seed(2)
#' # Highly overdispersed 50% resample design (alpha = 100)
#' overdispersed <- run_sampler_phy(tree, 5, alpha = 100, starting = "t10")
#' # Highly aggregated 50% resample design (alpha = -100)
#' aggregated <- run_sampler_phy(tree, 5, alpha = -100, starting = "t10")
#' # Random 50% resample design (alpha = 0)
#' random <- run_sampler_phy(tree, 5, alpha = 0, starting = "t10")
#' # Plot to compare
#' par(mfrow = c(2, 2))
#' plot(tree, main = "Full tree", cex = 1)
#' axis(1)
#' plot(overdispersed, main = "Overdispersed 50% sampling", cex = 1)
#' axis(1)
#' plot(aggregated, main = "Aggregated 50% sampling", cex = 1)
#' axis(1)
#' plot(random, main = "Random 50% sampling", cex = 1)
#' axis(1)


run_sampler_phy <- function (x, n, alpha, dist.func = cophenetic,
                             n_start = 1, return_start = FALSE,
                             starting = NULL) {
  #
  if (class(x) != 'phylo') {
    stop("x must be a phylogeny")
  }
  if (class(dist.func) != "function") {
    stop("dist.func must be a function")
  }
  tree <- x
  x <- dist.func(tree)
  if(!all(tree$tip.label %in% rownames(x))) {
    stop("dist.func provided do not result in a distance matrix
         with names corresponding to the tree tip labels.")
  }
  selected <- run_sampler(x, n, alpha, n_start, return_start,
                          starting)
  if (return_start) {
    selected[[1]] <- drop.tip(tree, tree$tip.label[!tree$tip.label %in% selected[[1]]])
    return(selected)
  } else {
    selected <- drop.tip(tree, tree$tip.label[!tree$tip.label %in% selected])
    return(selected)
  }
}


