#' Aggregated and overdispersed sampling
#'
#' The sampler package generates aggregated or overdispersed samples from any
#' given phylogenetic tree, spatial points defined by two coordinates, or
#' generic distance matrix. Resulting samples can remove data bias in existing
#' datasets or help design sampling schemes for new experiments.
#'
#' The package includes 4 functions:
#' \cr*\code{\link{run_sampler}}: Generates aggregated (underdispersed) or
#' overdispersed samples of values from any given distance matrix.
#' \cr*\code{\link{run_sampler_geo}}: Generates aggregated (underdispersed) or
#' overdispersed samples of points from any given set of geographic coordinates.
#' \cr*\code{\link{run_sampler_phy}}: Generates aggregated (underdispersed) or
#' overdispersed samples of tips from any given phylogenetic tree.
#' \cr*\code{\link{Nee_May_1997}}: Generates a subsample of tips from a
#' phylogenetic tree that maximize the remaining amount of evolutionary history
#'  for each node (see Nee and May 1997).

#'
#' @name sampler-package
#' @aliases sampler
#' @docType package
#'
#' @author \strong{Bruno Vilela} (maintainer), \strong{Ty Tuff}, \strong{Trevor Fristoe},
#' \strong{Edric Gavin}, \strong{Lucas Jardim}, \strong{Sara Varela}
#'  & \strong{Carlos Botero}
#'
#' @references Sean Nee and Robert M. May. "Extinction and the loss of evolutionary history." Science 278.5338 (1997): 692-694.
#'
#' @export run_sampler run_sampler_phy Nee_May_1997 run_sampler_geo
#' @import ape
#' @importFrom fields rdist.earth
#' @importFrom picante node.age
#' @importFrom phangorn Descendants Children
#' @importFrom stats cophenetic
#'
NULL

