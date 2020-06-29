\name{geographic_acf}
\alias{geographic_acf}
\title{
Phylogenetic autocorrelation function of geographic locations.
}
\description{
Given a rooted phylogenetic tree and geographic coordinates (latitudes & longitudes) of each tip, calculate the phylogenetic autocorrelation function (ACF) of the geographic locations. The ACF is a function of phylogenetic distance x, i.e., ACF(x) is the autocorrelation between two tip locations conditioned on the tips having phylogenetic ("patristic") distance x.
}
\usage{
geographic_acf(tree, tip_latitudes, tip_longitudes, Npairs=10000, Nbins=10)
}
\arguments{
\item{tree}{
A rooted tree of class "phylo". The root is assumed to be the unique node with no incoming edge.
}
\item{tip_latitudes}{
A numeric vector of size Ntips, specifying the latitudes (decimal degrees) of the tips. Note that \code{tip_latitudes[i]} (where i is an integer index) must correspond to the i-th tip in the tree, i.e. as listed in \code{tree$tip.label}.
}
\item{tip_longitudes}{
A numeric vector of size Ntips, specifying the longitudes (decimal degrees) of the tips. Note that \code{tip_longitudes[i]} (where i is an integer index) must correspond to the i-th tip in the tree, i.e. as listed in \code{tree$tip.label}.
}
\item{Npairs}{
Total number of random tip pairs to draw. A greater number of tip pairs will improve the accuracy of the estimated ACF within each distance bin. Tip pairs are drawn randomly with replacement. If \code{Npairs<=0}, then every tip pair is included exactly once.
}
\item{Nbins}{
Number of distance bins to consider within the range of phylogenetic distances encountered between tip pairs in the tree. A greater number of bins will increase the resolution of the ACF as a function of phylogenetic distance, but will decrease the number of tip pairs falling within each bin (which reduces the accuracy of the estimated ACF).
}
}
\details{
The autocorrelation between random geographic locations is defined as the expectation of \eqn{<X,Y>}, where <> is the scalar product and \eqn{X} and \eqn{Y} are the unit vectors pointing towards the two random locations on the sphere. For comparison, for a spherical Brownian Motion model with constant diffusivity \eqn{D} and radius \eqn{r} the autocorrelation function is given by \eqn{ACF(t)=e^{-2Dt/r^2}} (see e.g. \code{\link{simulate_sbm}}).

The phylogenetic autocorrelation function (ACF) of the geographic distribution of species can give insight into the dispersal processes shaping species distributions over global scales. An ACF that decays slowly with increasing phylogenetic distance indicates a strong phylogenetic conservatism of the location and thus slow dispersal, whereas a  rapidly decaying ACF indicates weak phylogenetic conservatism and thus fast dispersal. Similarly, if the mean distance between two random tips increases with phylogenetic distance, this indicates a phylogenetic autocorrelation of species locations. Here, phylogenetic distance between tips refers to their patristic distance, i.e. the minimum cumulative edge length required to connect the two tips.

Since the phylogenetic distances between all possible tip pairs do not cover a continuoum (as there is only a finite number of tips), this function randomly draws tip pairs from the tree, maps them onto a finite set of equally-sized distance bins and then estimates the ACF for the centroid of each distance bin based on tip pairs in that bin. In practice, as a next step one would usually plot the estimated ACF (returned vector \code{autocorrelations}) over the centroids of the distance bins (returned vector \code{distances}).

The tree may include multi-furcations (i.e. nodes with more than 2 children) as well as mono-furcations (i.e. nodes with only one child). If \code{tree$edge.length} is missing, then every edge is assumed to have length 1. The input tree must be rooted at some node for technical reasons (see function \code{\link{root_at_node}}), but the choice of the root node does not influence the result.
}
\value{
A list with the following elements:
\item{success}{
Logical, indicating whether the calculation was successful. If \code{FALSE}, an additional element \code{error} (character) is returned that provides a brief description of the error that occurred; in that case all other return values may be undefined.
}
\item{distances}{
Numeric vector of size Nbins, storing the centroid phylogenetic distance of each distance bin in increasing order. The first and last distance bin approximately span the full range of phylogenetic distances encountered between any two random tips in the tree.
}
\item{autocorrelations}{
Numeric vector of size Nbins, storing the estimated Pearson autocorrelation of the trait for each distance bin.
}
\item{mean_geodistances}{
Numeric vector of size Nbins, storing the mean geographic distance between tip pairs in each distance bin.
}
\item{Npairs_per_distance}{
Integer vector of size Nbins, storing the number of random tip pairs associated with each distance bin.
}
}


\author{Stilianos Louca}

%\references{
%J. R. Zaneveld and R. L. V. Thurber (2014). Hidden state prediction: A modification of classic ancestral state reconstruction algorithms helps unravel complex symbioses. Frontiers in Microbiology. 5:431.
%}

\seealso{
\code{\link{get_trait_depth}}, 
\code{\link{get_trait_acf}}
}

\examples{
# generate a random tree
tree = generate_random_tree(list(birth_rate_intercept=1),max_tips=1000)$tree

# simulate spherical Brownian Motion on the tree
simul = simulate_sbm(tree, radius=1, diffusivity=0.1)
tip_latitudes  = simul$tip_latitudes
tip_longitudes = simul$tip_longitudes

# calculate geographical autocorrelation function
ACF = geographic_acf(tree, tip_latitudes, tip_longitudes, Nbins=10)

# plot ACF (autocorrelation vs phylogenetic distance)
plot(ACF$distances, ACF$autocorrelations, type="l", xlab="distance", ylab="ACF")
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{trait evolution}
