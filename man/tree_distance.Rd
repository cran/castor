\name{tree_distance}
\alias{tree_distance}
\title{
Calculate the distance between two trees.
}
\description{
Given two rooted phylogenetic trees with identical tips, calculate their difference using a distance metric.
}
\usage{
tree_distance(treeA, treeB, tipsA2B=NULL, metric="RF", normalized=FALSE)
}
\arguments{
\item{treeA}{
A rooted tree of class "phylo". The root is assumed to be the unique node with no incoming edge.
}
\item{treeB}{
A rooted tree of class "phylo", with the same number of tips as \code{treeA}. The root is assumed to be the unique node with no incoming edge.
}
\item{tipsA2B}{
Optional integer vector of size Ntips, mapping \code{treeA} tip indices to \code{treeB} tip indices (i.e. \code{tipsA2B[a]} is the tip index in \code{treeB} corresponding to tip index \code{a} in \code{treeA}). The mapping must be one-to-one. If left unspecified, it is determined by matching tip labels between the two trees (this assumes that the same tip labels are used in both trees).
}
\item{metric}{
Character, specifying the distance metric to be used. Currently only the Robinson-Foulds ("RF") metric is implemented, which is the number of clusters (sets of tips descending from a node) in either of the trees but not shared by both trees (Robinson and Foulds, 1981; Day, 1985). The Robinson-Foulds metric does not take into account branch lengths.
}
\item{normalized}{
Logical, specifying whether the calculated distance should be normalized to be between 0 and 1. For the Robinson-Foulds metric, the distance will be normalized by dividing it by the total number of nodes in the two trees.
}
}
\details{
If the trees differ in theis tips, they must be pruned down to their common set of tips. If tips have different labels in the two trees, but are nevertheless equivalent, the mapping between the two trees must be provided using \code{tipsA2B}. The trees may include multi-furcations as well as mono-furcations (i.e. nodes with only one child).

Note that under some Robinson-Foulds variants the trees can be unrooted; in this present implementation trees must be rooted and the placement of the root influences the distance, following the definition by Day (1985).

The asymptotic average time complexity of this function is O(Nedges*log(Nedges)*log(log(Nedges))) for a balanced bifurcating tree.
}

\value{
A single non-negative number, representing the distance between the two trees.
}


\author{Stilianos Louca}

\references{
Robinson, D. R., Foulds, L. R. (1981). Comparison of phylogenetic trees. Mathematical Biosciences. 53: 131-147.

Day, W. H. E. (1985). Optimal algorithms for comparing trees with labeled leaves. Journal of Classification. 2:7-28.
}

\seealso{
\code{\link{congruent_divergence_times}}
}

\examples{
# generate a random tree
Ntips = 1000
treeA = generate_random_tree(list(birth_rate_intercept=1),
                            max_tips=Ntips)$tree
                            
# create a second tree with slightly different topology
treeB = treeA
shuffled_tips = sample.int(Ntips, size=Ntips/10, replace=FALSE)
treeB$tip.label[shuffled_tips] = treeB$tip.label[sample(shuffled_tips)]

# calculate Robinson-Foulds distance between trees
distance = tree_distance(treeA, treeB, metric="RF")
}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{tree comparison}