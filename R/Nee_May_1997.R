#' Algorithms that maximize the amount of evolutionary history in a subsample
#'
#' @author Bruno Vilela
#'
#' @description Generate a subsample from a phylogenetic tree that maximizes
#' the remaining amount of evolutionary history (see Nee and May 1997).
#'
#' @param x \code{phylo} object.
#' @param n A positive integer number indicating the sample size.
#'
#'
#' @details The algorithm selects n - 1 oldest nodes, defining n clades. One species
#' from each clade is randonly selected (Nee and May 1997).
#'
#'
#' @return A prunned phylogenetic tree.
#'
#' @references Sean Nee and Robert M. May. "Extinction and the loss of evolutionary history." Science 278.5338 (1997): 692-694.
#'
#' @examples
#' require(ape)
#' set.seed(50)
#' tree <- rcoal(50)
#' k = 10
#' par(mfrow = c(1, 2))
#' plot(tree, main = paste("Original tree \n n =", Ntip(tree)))
#' plot(Nee_May_1997(tree, k), main = paste("Tree with phylogenetic history maximized
#' n =",  k,  "(Nee & May 1997)"))
#'
#' # Compare the sampler algorithm with Nee and May (1997) method:
#' PD_ob <- sum(Nee_May_1997(tree, k)$edge.length)
#' n_a = 20
#' alphas <- seq(-5, 5, length.out = n_a)
#' rep = 100
#' PD <- matrix(ncol = n_a, nrow = rep)
#' for (j in 1:rep) {
#'  for (i in 1:n_a) {
#'  PD[j, i] <- sum(run_sampler_phy(tree, k, alpha = alphas[i], starting = "t10")$edge.length)
#'  }
#' }
#' PD_max <- apply(PD, 2, max)
#' PD_min <- apply(PD, 2, min)
#' PD_mean <- apply(PD, 2, mean)
#' par(mfrow = c(1, 1))
#' plot(alphas, PD_mean, ylab = "Phylogenetic Diveristy",
#'  ylim = c(min(PD), max(c(PD, PD_ob))), type = "l", lwd = 2, col = "black")
#'  lines(alphas, PD_max, lty = 2, col = "gray20")
#'  lines(alphas, PD_min, lty = 2, col = "gray20")
#'  points(max(alphas) + max(alphas) *.05, PD_ob, col = "red", pch = 18, cex = 2)
#'  legend(x = max(alphas) - max(alphas) * .90, y = min(PD) + min(PD) * 2,
#'        legend = c("Sampler (mean)", "Sampler (min and max)", "Nee & May"),
#'        pch = c(NA, NA, 18), lty = c(1, 2, NA),
#'        col = c("black", "gray20","red"))


Nee_May_1997 <- function(x, n) {
  if (class(x) != 'phylo') {
    stop("x must be a phylogeny of class phylo")
  }
  if (!is.numeric(n)) {
    stop("n must be a positive integer")
  }
  if (length(n) > 1) {
    stop("n must be a positive integer of length one")
  }
  if (n < 2 | n > Ntip(x) | n != round(n)) {
    stop("n must be a positive integer of length one, bigger than 2,
         and smaller than the number of tips in the tree")
  }

  nodes <- 1:(Nnode(x) + Ntip(x))
  x <- node.age(x)
  ages <- c(0, x$ages)
  pos0 <- c(x$edge[(!x$edge[, 1] %in% x$edge[, 2]), 1][1],
            x$edge[, 2])
  pos <- pos0[order(pos0)]
  ages <- ages[order(pos0)]
  pick <- order(ages)[1:(n-1)]
  children <- unlist(Children(x, pick))
  choices <- children[!children %in% pick]
  choices2 <- Descendants(x, choices)
  chosen <- sapply(choices2, function(x)if(length(x) > 1){sample(x, 1)}else{x})
  x2 <- drop.tip(x, x$tip.label[-chosen])
  return(x2)
}

