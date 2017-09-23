\name{count_clades_over_time}
\alias{count_clades_over_time}
\title{
Count number of clades over time.
}
\description{
Given a rooted phylogenetic tree whose edge lengths represent time intervals, calculate the number of clades at various time points (e.g., spanning from 0 to the maximum time of any tip). The root is interpreted as time 0, and the distance of any node or tip from the root is interpreted as time elapsed since the root. This function defines an equidistant sequence of time points, and counts how many edges "cross" each time point. Optionally, the slopes and relative slopes of the clade-counts-vs-time curve are also returned. The slopes and relative slopes are approximations for the species birth rate and the per-capita species birth rate (assuming no extinctions occurred).
}
\usage{
count_clades_over_time(tree, Ntimes=NULL, times=NULL, include_slopes=FALSE)
}
\arguments{
\item{tree}{
A rooted tree of class "phylo", where edge lengths represent time intervals (or similar). The root is assumed to be the unique node with no incoming edge.
}
\item{Ntimes}{
Integer, number of equidistant time points for which to calculade clade counts. Can also be \code{NULL}, in which case \code{times} must be provided.
}
\item{times}{
Integer vector, listing time points (in ascending order) for which to calculate clade counts. Can also be \code{NULL}, in which case \code{Ntimes} must be provided.
}
\item{include_slopes}{
Logical, specifying whether the slope and the relative slope of the returned clades-per-time-point curve should also be returned.
}
}


\details{
If \code{tree$edge.length} is missing, then every edge in the tree is assumed to be of length 1. The tree may include multi-furcations as well as mono-furcations (i.e. nodes with only one child). The tree need not be ultrametric, although in general this function only makes sense for dated trees (e.g., where edge lengths are time intervals or similar).

Either \code{Ntimes} or \code{times} must be non-\code{NULL}, but not both.
}


\value{
A list with the following elements:
\item{Ntimes}{
Integer, indicating the number of returned time points. Equal to the provided \code{Ntimes} if applicable.
}
\item{times}{
Numeric vector of size Ntimes, listing the considered time points in increasing order. If \code{times} was provided as an argument to the function, then this will be the same as provided.
}
\item{diversities}{
Integer vector of size Ntimes, listing the number of clades for each time point.
}
\item{slopes}{
Numeric vector of size Ntimes, listing the slopes (finite-difference approximation of 1st derivative) of the curve clade_counts vs time_point.
}
\item{relative_slopes}{
Numeric vector of size Ntimes, listing the relative slopes of the curve clade_counts vs time_point, i.e. \code{slopes} divided by a sliding-window average of \code{clade_counts}.
}
}


\author{Stilianos Louca}

%\references{
%}

%\seealso{
%}

\examples{
# generate a random tree
tree = generate_random_tree(list(birth_rate_intercept=1), max_tips=1000)$tree

# count clades over time
results = count_clades_over_time(tree, Ntimes=100)

# plot curve (number of clades vs time)
plot(results$times, results$diversities, type="l", xlab="time", ylab="# clades")
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{trait evolution}
