\name{count_clades_over_time}
\alias{count_clades_over_time}
\title{
Count number of clades over time.
}
\description{
Given a rooted phylogenetic tree whose edge lengths represent time intervals, calculate the number of clades at various time points (spanning from 0 to the maximum time of any tip). The root is interpreted as time 0, and the distance of any node or tip from the root is interpreted as time elapsed since the root. This function defines an equidistant sequence of time points, and counts how many edges "cross" each time point. Optionally, the slopes and relative slopes of the clade-counts-vs-time curve are also returned. The slopes and relative slopes are approximations for the species birth rate and the per-capita species birth rate (assuming no extinctions occurred).
}
\usage{
count_clades_over_time(tree, Ntimes, include_slopes=FALSE)
}
\arguments{
\item{tree}{
A rooted tree of class "phylo", where edge lengths represent time intervals (or similar). The root is assumed to be the unique node with no incoming edge.
}
\item{Ntimes}{
Integer, number of equidistant time points to consider.
}
\item{include_slopes}{
Logical, specifying whether the slope and the relative slope of the returned clades-per-time-point curve should also be returned.
}
}


\details{
If \code{tree$edge.length} is missing, then every edge in the tree is assumed to be of length 1. The tree may include multi-furcations as well as mono-furcations (i.e. nodes with only one child). The tree need not be ultrametric, although in general this function only makes sense for dated trees (e.g., where edge lengths are time intervals or similar).
}


\value{
A list with the following elements:
\item{time_points}{
Numeric vector of size Ntimes, listing the considered time points in increasing order.
}
\item{clade_counts}{
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
# generate a random tree according to the Yule model
tree = generate_random_tree(max_tips=1000, birth_rate_intercept=1)

# count clades over time
results = count_clades_over_time(tree, Ntimes=100)

# plot curve (number of clades vs time)
plot(results$time_points, results$clade_counts, type="l", xlab="time", ylab="# clades")
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ Prunning }
