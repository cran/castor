\name{get_all_distances_to_tip}
\alias{get_all_distances_to_tip}
\title{
Get distances of all tips/nodes to a focal tip.
}
\description{
Given a tree and a focal tip, calculate the phylogenetic ("patristic") distances between the focal tip and all other tips & nodes in the tree.
}
\usage{
get_all_distances_to_tip(tree, focal_tip)
}
\arguments{
\item{tree}{
A rooted tree of class "phylo".
}
\item{focal_tip}{
Either a character, specifying the name of the focal tip, or an integer between 1 and Ntips, specifying the focal tip's index.
}
}
\details{
The "patristic distance" between two tips and/or nodes is the shortest cumulative branch length that must be traversed along the tree in order to reach one tip/node from the other. If \code{tree$edge.length} is missing, then each edge is assumed to be of length 1.

The tree may include multi-furcations as well as mono-furcations (i.e. nodes with only one child). The input tree must be rooted at some node for technical reasons (see function \code{\link{root_at_node}}), but the choice of the root node does not influence the result.
}
\value{
A numeric vector of length Ntips+Nnodes, specifying the distances of all tips (entries 1,..,Ntips) and all nodes (entries Ntips+1,..,Ntips+Nnodes) to the focal tip.
}


\author{Stilianos Louca}

%\references{
%}

\seealso{
\code{\link{get_all_pairwise_distances}},
\code{\link{get_pairwise_distances}}
}

\examples{
# generate a random tree
Ntips = 100
tree  = generate_random_tree(list(birth_rate_intercept=1),
                             max_tips = Ntips,
                             tip_basename="tip.")$tree

# calculate all distances to a focal tip
distances = get_all_distances_to_tip(tree, "tip.39")
print(distances)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{phylogenetic distance}
