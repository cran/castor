\name{correlate_phylo_geodistances}
\alias{correlate_phylo_geodistances}
\title{
Correlations between phylogenetic & geographic distances.
}
\description{
Given a rooted phylogenetic tree and geographic coordinates (latitudes & longitudes) of each tip, examine the correlation between pairwise phylogenetic and geographic distances of tips. The statistical significance is computed by randomly permuting the tip coordinates, which is essentially a phylogenetic version of the Mantel test and accounts for shared evolutionary history between tips.
}
\usage{
correlate_phylo_geodistances(tree,
                             tip_latitudes,
                             tip_longitudes,
                             correlation_method,
                             max_phylodistance  = Inf,
                             Npermutations      = 1000,
                             alternative        = "right",
                             radius             = 1)
}
\arguments{
\item{tree}{
A rooted tree of class "phylo".
}
\item{tip_latitudes}{
Numeric vector of size Ntips, specifying the latitudes (decimal degrees) of the tree's tips. By convention, positive latitudes correspond to the northern hemisphere. Note that \code{tip_latitudes[i]} must correspond to the i-th tip in the tree, i.e. as listed in \code{tree$tip.label}.
}
\item{tip_longitudes}{
Numeric vector of size Ntips, specifying the longitudes (decimal degrees) of the tree's tips. By convention, positive longitudes correspond to the eastern hemisphere.
}
\item{correlation_method}{
Character, one of "\code{pearson}", "\code{spearman}" or "\code{kendall}".
}
\item{max_phylodistance}{
Numeric, maximum phylodistance between tips to consider, in the same units as the tree's edge lengths. If \code{Inf}, all tip pairs will be considered.
}
\item{Npermutations}{
Integer, number of random permutations to consider for estimating the statistical significance (P value). If 0, the significance will not be computed. A larger number improves accuracy but at the cost of increased computing time.
}
\item{alternative}{
Character, one of "\code{two_sided}", "\code{right}" or "\code{left}", specifying which part of the null model's distribution to use as P-value.
}
\item{radius}{
Optional numeric, radius to assume for the sphere. If 1, then all geodistances are measured in multiples of the sphere radius. This does not affect the correlation or P-value, but it affects the returned geodistances. Note that Earth's average radius is about 6371 km.
}
}
\details{
To compute the statistical significance (P value) of the observed correlation \code{C}, this function repeatedly randomly permutes the tip coordinates, each time recomputing the corresponding "random" correlation, and then examines the distribution of the random correlations. If \code{alternative="right"}, the P value is set to the fraction of random correlations equal to or greater than \code{C}.
}
\value{
A named list with the following elements:
\item{correlation}{
Numeric between -1 and 1, the correlation between phylodistances and geodistances.
}
\item{Npairs}{
Integer, the number of tip pairs considered.
}
\item{Pvalue}{
Numeric between 0 and 1, estimated statistical significance of the correlation. Only returned if \code{Npermutations>0}.
}
\item{mean_random_correlation}{
Numeric between -1 and 1, the mean correlation obtained across all random permutations. Only returned if \code{Npermutations>0}.
}
\item{phylodistances}{
Numeric vector of length \code{Npairs}, listing the pairwise phylodistances between tips, used to compute the correlation.
}
\item{geodistances}{
Numeric vector of length \code{Npairs}, listing the pairwise geodistances between tips, used to compute the correlation.
}
}


\author{Stilianos Louca}

%\references{
%J. R. Zaneveld and R. L. V. Thurber (2014). Hidden state prediction: A modification of classic ancestral state reconstruction algorithms helps unravel complex symbioses. Frontiers in Microbiology. 5:431.
%}

\seealso{
\code{\link{geographic_acf}}
}

\examples{
# Generate a random tree
Ntips = 50
tree = generate_random_tree(list(birth_rate_intercept=1),max_tips=Ntips)$tree

# simulate spherical Brownian motion (a dispersal model) on the tree
simul = simulate_sbm(tree, radius=6371, diffusivity=50)

# Analyze correlations between geodistances & phylodistances
coranal = correlate_phylo_geodistances(tree               = tree,
                                       tip_latitudes      = simul$tip_latitudes,
                                       tip_longitudes     = simul$tip_longitudes,
                                       correlation_method = "spearman",
                                       Npermutations      = 100,
                                       max_phylodistance  = 100,
                                       radius             = 6371)
print(coranal$correlation)
print(coranal$Pvalue)
plot(coranal$phylodistances, coranal$geodistances,
     xlab="phylodistance", ylab="geodistance", type="p")
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{phylogeography}
\keyword{dispersal}
