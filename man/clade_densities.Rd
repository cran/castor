\name{clade_densities}
\alias{clade_densities}
\title{
Estimate the density of tips & nodes in a timetree.
}
\description{
Given a rooted timetree (i.e., a tree whose edge lengths represent time intervals), estimate the density of tips and nodes as a function of age (tips or nodes per time unit), on a discrete grid. Optionally the densities can be normalized by the local number of lineages. If the tree is full (includes all extinct & extant clades), then the normalized tip (node) density is an estimate for the per-capita extinction (speciation) rate.
}
\usage{
clade_densities(tree, 
                Nbins       = NULL,
                min_age     = NULL,
                max_age     = NULL,
                normalize   = TRUE,
                ultrametric	= FALSE)		
}
\arguments{
\item{tree}{
A rooted tree of class "phylo", where edge lengths represent time intervals (or similar).
}
\item{Nbins}{
Integer, number of equidistant age bins at which to calculate densities.
}
\item{min_age}{
Numeric, minimum age to consider. If \code{NULL}, it will be set to the minimum possible.
}
\item{max_age}{
Numeric, maximum age to consider. If \code{NULL}, it will be set to the maximum possible.
}
\item{normalize}{
Logical, whether to normalize densities by the local number of lineages (in addition to dividing by the age interval). Hence, tip (or node) densities will represent number of tips (or nodes) per time unit per lineage.
}
\item{ultrametric}{
Logical, specifying whether the input tree is guaranteed to be ultrametric, even in the presence of some numerical inaccuracies causing some tips not have exactly the same distance from the root.
}
}


\details{
This function discretizes the full considered age range (from \code{min_age} to \code{max_age}) into \code{Nbins} discrete disjoint bins, then computes the number of tips and nodes in each bin, and finally divides those numbers by the bin width. If \code{normalize==True}, the densities are further divided by the number of lineages in the middle of each age bin. For typical timetrees it is generally recommended to omit the most recent tips (i.e., extant at age 0), by setting \code{min_age} to a small non-zero value; otherwise, the first age bin will typically be associated with a high tip density, i.e., tip densities will be zero-inflated.
}


\value{
A list with the following elements:
\item{Nbins}{
Integer, indicating the number of discrete age bins.
}
\item{ages}{
Numeric vector of size Nbins, listing the centres of the age bins.
}
\item{tip_densities}{
Numeric vector of size Nbins, listing the tip densities in the corresponding age bins.
}
\item{node_densities}{
Numeric vector of size Nbins, listing the node densities in the corresponding age bins.
}
}


\author{Stilianos Louca}

%\references{
%}

\seealso{
count_lineages_through_time
}

\examples{
# generate a random full tree, including all extinct & extant tips
tree = generate_random_tree(list(birth_rate_intercept=1), 
                            max_tips=1000, coalescent=FALSE)$tree

# compute node densities, as an estimate for the speciation rate
densities = clade_densities(tree, Nbins=10, normalize=TRUE)

# plot node densities
plot(densities$ages, densities$node_densities, type="l", xlab="age", ylab="node density")
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{lineages through time}
