% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/run_sampler_geo.R
\name{run_sampler_geo}
\alias{run_sampler_geo}
\title{Aggregated and overdispersed sampling for geographic coordinates}
\usage{
run_sampler_geo(x, n, alpha, dist.func = rdist.earth, n_start = 1,
  return_start = FALSE, starting = NULL)
}
\arguments{
\item{x}{\code{matrix} or \code{data.frame} indicating the coordinates
(first column = longitude; second column = latitude). Row names should be given.}

\item{n}{A positive integer number indicating the sample size.}

\item{alpha}{Number indicating the strength of aggregation (if negative) or
overdispersion (if positive). When alpha = 0 sample is random.}

\item{dist.func}{A distance function to calculate coordinates distance.
Default is \code{\link{rdist.earth}} from package \code{fields}.}

\item{n_start}{Number of initial selected points. Default is one starting point.}

\item{return_start}{if \code{TRUE} the starting point is returned.}

\item{starting}{Character vector indicating the starting point (= to row names).
If not provided, random starting value(s) is(are) selected.}
}
\value{
The function returns a vector indicating the selected rows.
If return_start is TRUE, a list is returned with the first element being the
Sampling_selection - selected sampling units - and
Starting_points - selected starting point(s).
}
\description{
Generate aggregated or overdispersed sampling designs
for any given coordinates. Results can
be used to design experiments/samples, for resample proposes and data bias removal.
}
\details{
The function uses the algorithm in \code{\link{run_sampler}},
 but here it accepts a two column matrix of coordinates as input.
}
\examples{
require(sp)
require(maptools)
data(wrld_simpl)  # World map
Brazil <- wrld_simpl[wrld_simpl$NAME == "Brazil", ]  # Brazil (polygon)
coords <- slot(spsample(Brazil, 100, "regular"), "coords")
rownames(coords) <- paste0("t", 1:nrow(coords))
## Subsample 50\%
### Overdispersed
selection.geo <- run_sampler_geo(x = coords, n = 10, alpha = 100, starting = "t10")
### Aggregated
selection.geo2 <- run_sampler_geo(x = coords, n = 10, alpha = -100, starting = "t10")
### Random
selection.geo3 <- run_sampler_geo(x = coords, n = 10, alpha = 0, starting = "t10")

## Plot
par(mfrow = c(1, 3), mar = c(1, 1, 15, 1))
plot(Brazil, main = "Overdispersed 50\% sampling (red were selected)")
points(selection.geo, cex = 2, pch = 19, col = "red")
plot(Brazil, main = "Aggregated 50\% sampling (blue were selected)")
points(selection.geo2, cex = 2, pch = 19, col = "blue")
plot(Brazil, main = "Random 50\% sampling (green were selected)")
points(selection.geo3, cex = 2, pch = 19, col = "green")
}
\seealso{
\code{\link{run_sampler}}
}
\author{
Bruno Vilela
}
