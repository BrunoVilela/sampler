#' Aggregated and overdispersed sampling
#'
#' The sampler package is meant to Generate aggregated or overdispersed sampling
#' designs for any given phylogenetic tree, spatial points or distance matrix.
#' Results can be used to design experiments/samples, for resample proposes and
#' data bias removal.
#'
#' The package includes 4 functions:
#' \cr*\code{\link{run_sampler}}: Generate aggregated or overdispersed
#' sampling designs for any given distance matrix.
#' \cr*\code{\link{run_sampler_geo}}: Applies the run_sampler algorithm to
#'  geographic coordinates.
#' \cr*\code{\link{run_sampler_phy}}: Applies the run_sampler algorithm to
#'  phylogenetic trees.
#' \cr*\code{\link{Nee_May_1997}}: Generate a subsample from a phylogenetic
#' tree that maximizes the remaining amount of evolutionary history
#' (see Nee and May 1997).
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
#' @import ape fields sp picante phangorn maptools
#' @importFrom stats cophenetic
#'
NULL

