\name{find_root}
\alias{find_root}
\title{
Find the root of a tree.
}
\description{
Find the root of a phylogenetic tree. The root is defined as the unique node with no parent.
}
\usage{
find_root(tree)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
\item{tree}{
A tree of class "phylo".
}
}


\details{
By convention, the root of a "phylo" tree is typically the first node (i.e. with index \code{length(tree$tip.label)+1}), however this is not always guaranteed. This function finds the root of a tree by searching for the node with no parent. If no root exists, NA is returned. If multiple roots exist (this should generally not happen), only one of them will be returned.
The asymptotic time complexity of this function is O(Nedges), where Nedges is the number of edges in the tree.
}

\value{
Index of the tree's root, as listed in \code{tree$edge}. An integer ranging from Ntips+1 to Ntips+Nnodes, where Ntips and Nnodes is the number of tips and nodes in the tree, respectively. Typically (but not always), the root will be Ntips+1.
}


\author{Stilianos Louca}

%\references{
%}

\seealso{
\code{\link{root_at_node}}
}

\examples{
# generate a random tree
Ntips = 100
tree = generate_random_tree(list(birth_rate_intercept=1),Ntips)$tree

# reroot the tree at the 20-th node
new_root_node = 20
tree = root_at_node(tree, new_root_node, update_indices=FALSE)

# find new root index and compare with expectation
cat(sprintf("New root is \%d, expected at \%d\n",find_root(tree),new_root_node+Ntips))
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ Pruning }
