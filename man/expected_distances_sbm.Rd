\name{expected_distances_sbm}
\alias{expected_distances_sbm}
\title{
Expected distances traversed by a Spherical Brownian Motion.
}
\description{
Given a Spherical Brownian Motion (SBM) process with constant diffusivity, compute the expected geodesic distance traversed over specific time intervals. This quantity may be used as a measure for how fast a lineage disperses across the globe over time.
}
\usage{
expected_distances_sbm(diffusivity,
                       radius,
                       deltas)
}
\arguments{
\item{diffusivity}{
Numeric, the diffusivity (aka. diffusion coefficient) of the SBM. The units of the diffusivity must be consistent with the units used for specifying the \code{radius} and time intervals (\code{deltas}); for example, if \code{radius} is in km and \code{deltas} are in years, then \code{diffusivity} must be specified in km^2/year.
}
\item{radius}{
Positive numeric, the radius of the sphere.
}
\item{deltas}{
Numeric vector, listing time intervals for which to compute the expected geodesic distances.
}
}
\details{
This function returns expected geodesic distances (i.e. accounting for spherical geometry) for a diffusion process on a sphere, with isotropic and homogeneous diffusivity. 
}

\value{
A non-negative numeric vector of the same length as \code{deltas}, specifying the expected geodesic distance for each time interval in \code{deltas}.
}


\author{Stilianos Louca}

\references{
S. Louca (2021). Phylogeographic estimation and simulation of global diffusive dispersal. Systematic Biology. 70:340-359.
}

\seealso{
\code{\link{fit_sbm_const}}
}

\examples{
# compute the expected geodistance (in km) after 100 and 1000 years
# assuming a diffusivity of 20 km^2/year
expected_distances = expected_distances_sbm(diffusivity = 20,
                                            radius      = 6371,
                                            deltas      = c(100,1000))
}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{SBM model}
