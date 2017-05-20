\name{collapse_tree_at_resolution}
\alias{collapse_tree_at_resolution}
\title{
Collapse nodes of a tree at a phylogenetic resolution.
}
\description{
Given a rooted tree and a phylogenetic resolution threshold, collapse all nodes whose distance to all descending tips does not exceed the threshold, into new tips. This function can be used to obtain a "coarser" version of the tree, or to collapse identical tips.
}
\usage{
collapse_tree_at_resolution(tree, resolution=0, 
                            by_edge_count=FALSE)
}
\arguments{
\item{tree}{
A rooted tree of class "phylo". The root is assumed to be the unique node with no incoming edge.
}
\item{resolution}{
Numeric, specifying the phylogenetic resolution at which to collapse the tree. This is the maximum distance a descending tip can have from a node, such that the node is collapsed into a new tip. If set to 0 (default), then only nodes whose descending tips are identical to the node will be collapsed.
}
\item{by_edge_count}{
Logical. Instead of considering edge lengths, consider edge counts as phylogenetic distance between nodes and tips. This is the same as if all edges had length equal to 1.
}
}


\details{
The input tree may include multi-furcations (i.e. nodes with more than 2 children) as well as mono-furcations (i.e. nodes with only one child).
}


\value{
A list with the following elements:
\item{collapsed_tree}{
A new rooted tree of class "phylo", containing the collapsed tree. Tip and node labels are inheritted from the original tree. Tip labels for collapsed nodes will be the node labels from the original tree.
}
\item{new2old_clade}{
Integer vector of length equal to the number of tips+nodes in the collapsed tree, with values in 1,..,Ntips+Nnodes, mapping tip/node indices of the collapsed tree to tip/node indices in the original tree. In particular,

\code{c(tree$tip.label,tree$node.label)[new2old_clade]}

will be equal to:

\code{c(collapsed_tree$tip.label,collapsed_tree$node.label)}.
}
\item{new2old_edge}{
Integer vector of length equal to the number of edges in the collapsed tree, with values in 1,..,Nedges, mapping edge indices of the collapsed tree to edge indices in the original tree. In particular, \code{tree$edge.length[new2old_edge]} will be equal to \code{collapsed_tree$edge.length} (if edge lengths are available).
}
}

\author{Stilianos Louca}

%\references{
%}

%\seealso{
%}

\examples{
# generate a random tree
tree = generate_random_tree(max_tips=1000, max_time=100, birth_rate_intercept=1);

# print number of nodes
cat(sprintf("Simulated tree has \%d nodes\n",tree$Nnode))

# collapse any nodes with tip-distances < 20
collapsed = collapse_tree_at_resolution(tree, resolution=20)$collapsed_tree

# print number of nodes
cat(sprintf("Collapsed tree has \%d nodes\n",collapsed$Nnode))
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{prunning}
\keyword{subtree}
